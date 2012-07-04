function [iter,v_vecs,y_vecs,mass_infeas,mass_action_infeas] ...
         = SolverFpIterationHomPDCO(Y,Ak, ...
                   mass_infeas_stop_tol,max_iter,IterStep,x_ini)

% Solves the fixed point problem
%   v*(r,r_0),v*_0(r,r_0) = (r,r_0),          (v* = v^\star)
%   v*_0 = r_0
% where v* is given by
%   min  v' D(log(v)-1)
%   st   YDv =  YA^Tr;

% Parameters
% Y            Stoichiometric Matrix
% At           Graph laplacian
% mass_action_infeas_stop_tol
%              Maximum infeasibility allowed for the linear constraint
% max_iter     Maximum number of iterations for the fixed point

% Returned Values:
% iter         Number of iterations.
% v_vecs       Minimizer at every iteration.
% y_vecs       Dual variable of the linear constraint at every iteration.
% mass_infeas  Vector with the norm of the residual of the linear
%              constraints at every iterate.
% mass_action_infeas
%              Norm of the vector Y^T\lambda - \log(v)

% Read the problem size
m = size(Y,1);
n = size(Ak,1);
if (size(Y,2)~=n)
   error('The stoichiometry matrix and the graph laplacian are of inconsistent sizes')
end

% Check if the initial point is defined.
if nargin<9
    x_ini = ones(n,1);
else
    if length(x_ini)~=n
        error('The initial point must be of size n\n');
    end
end

% Extract the adjacency matrix and the diagonal matrix.
d     = -diag(Ak);
At    = Ak+diag(d);
% Form the products used to define the constraints
YD    = Y*diag(d);
YAt   = Y*At;

% Calculate the constraints for the first iteration.
mu    = YAt*x_ini(1:end);


% Given the value of x_ini we can easily construct
% a feasible point.
v_feas = At*(x_ini)./d;


% Evaluate the norm of the initial infeasibility.
InInf  = norm(YD*v_feas - mu,inf);
fprintf('Infeasibility of initial point')
fprintf('\n||YDv_0-mu||_2: %d\n', InInf)

% Set up the pdco arguments
pdObj = @(x)phiP(x,d);   % Objective function returns f,g,diag(H)
pdMat = YD;          % Linear constraint
bl    = zeros(n,1);        % Lower bound on variables
bu    = Inf(n,1);          % Upper bounds on variables
d1    = 1e-4;                % Primal regularization
d2    = 1e-4;                % Dual regularization
options         = pdcoSet;   % PDCO options
options.MaxIter = 100;       % Maximum iterations in internal loop
options.Print   = 0;         % Option to suppress output
options.Method 	= 1;         % Use cholesky to solve the Newton systems
options.mu0     = 1;         % 0 let's PDCO set the initial mu

y      = ones(m,1);          % Initial dual of the equality constraint
z      = zeros(n,1);       % Initial dual of bounds	
% xsize = max(max(abs(YD)));% Estimate of the largest x at solution
xsize  = 1;
zsize  = 1;                  % Estimate of the largest z at solution
                             % (really the largest y)

iter    = 1;
hdr_str = '%4s %10s %10s %10s %10s %10s %10s %10s %10s %6s %6s %10s\n';
out_str = '%4d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %6d %6d %10.2e\n';
fprintf('\n')
fprintf(hdr_str,'iter','InitialInf','min(v*)','max(v*)',   ...
        'min(z)','max(z)','Dual Inf ','phi(v)','||YAk*v||', ...
        'inform','PDitns','time')

% Matrices to store the iterates
mass_infeas        = mass_infeas_stop_tol + 1;
mass_action_infeas = [];
objectives         = [];
mu_vecs            = [];
v_vecs             = [];
y_vecs             = [];

% Variables to calculate the initial dual infeasibility.
rho        = pdMat'*ones(m,1);
[f,g,h]    = pdObj(v_feas);
DualInfeas = norm(pdMat'*y + z - g,inf);

% Start the main iteration loop
while mass_infeas(iter) > mass_infeas_stop_tol && iter <= max_iter
    t1 = tic;
	
    % Call pdco with the options above, and the given parameter mu
    x0 = v_feas;
    y0 = y;
    z0 = z;
    [v,y,z,inform,PDitns,CGitns,time] = ...
        pdco(pdObj,pdMat,mu,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize);

    % Store the objective value
    objectives = [objectives, pdObj(v)];
    % Set v equal to the solution at each iteration and store them
    v_vecs     = [v_vecs v];
    y_vecs     = [y_vecs y];

    mass_infeas = [mass_infeas, norm(Y*Ak*v,inf)];	
    mass_action_infeas = [mass_action_infeas, norm(Y'*y - log(v),inf)];
    t2 = toc(t1);
    fprintf(out_str,iter,InInf,min(v),max(v),min(z),max(z), ...
            DualInfeas,objectives(iter),mass_infeas(iter),inform,PDitns,t2);

    if inform>0
       break     % pdco failed
    end

    % Define the new linear constraint.
    mu      = (1-IterStep)*mu + IterStep*(YAt*v);
    % Store the new constraint
    mu_vecs  = [mu_vecs mu];

    % Find a new initial feasible point, use a convex combination of the 
    % past initial feasible point and one calculated from the new minimizer.
	
    % Calculate a new point from the new minimizer. This point would 
    % be feasible if IterStep = 1.
	
    v_feasN  = At*(v)./d;
    
    % Make the convex combination of the last initial feasible point and the 
    % point v_feasN
    v_feas   = (1-IterStep)*v + IterStep*v_feasN;
    InInf    = norm(YD*v_feas-mu,inf);

    % increase iteration count
    iter = iter + 1;
end

mass_infeas = mass_infeas(2:iter);
iter        = iter - 1;



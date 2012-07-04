%This function executes the fixed point iteration for a provided model and stops at convergence
%Use 
function [iter,v_vecs,lmda_vecs,mass_infeas,mass_action_infeas,mu_vecs,objectives]= FpIterationPDCO(Y,Ak,mass_infeas_stop_tol,max_iter,x0)

m = size(Y,1);
n = size(Ak,1);

%Check if v0 is assigned
if(nargin==4)
		x0 = ones(n,1);
end
d = -diag(Ak)
A = Ak'-diag(d);
YAt = Y*A';
YD = Y*diag(d);

%inital parameter of the iterations
mu = YD*x0; 		

% solve with pdco
% pdco solves the problem 
%
%    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
%
% input parameters and options:
%
pdMat           = YD;            % constraint matrix, A
b               = mu;               % RHS of constraint
bl              = zeros(n,1);     % lower bound on x
bu              = Inf(n,1);       % upper bound on x
d1              = 1e-4;             % primal regularization 
d2              = 1e-4;             % dual regularization 
options         = pdcoSet;          % options for pdco
options.MaxIter = 500;              % increase the max iterations
options.Print   = 0;                % i hate when it prints shit out.
y0              = ones(m,1);        % init lagrange mult of linear constraint
z0              = zeros(n,1);     % init lagrange mult of box constraints
xsize           = max(max(abs(YD))); % estimates of the biggest x at solution
zsize           = 1;                % estimates of the biggest z at solution

% initial value of v and pdco iteration
v         = x0;
lmda      = y0;
z         = z0;
iter      = 1;
%v_vecs    = v;
%lmda_vecs = lmda;
mu_vecs   = mu;
gm        = 1;
fact      = 1.5;
vtol      = 1e-2;

% header for output 
hdr_str = '%4s %10s %10s %10s %10s %10s %10s\n';
out_str = '%2d %20d %10d %13.3d %10d %10d %10g\n';
fprintf(hdr_str,'iter','Initial Infeasibility','phi(v)','||YAk*v||','solved','inform','time');

% store the difference in norm between YAt*v and mu
diffnorm = [];
mass_infeas = mass_infeas_stop_tol + 1;
mass_action_infeas = [];
objectives = [];
mu_vecs = [];
v_vecs = [];
lmda_vecs = [];
% start timing
tic
% iterate mu until mu = Y*A'*v or we reach the maximum number of iterations
%while(( mass_infeas(iter) > mass_infeas_stop_tol || mass_action_infeas(iter) > mass_action_infeas_stop_tol) && iter < max_iter)
while(( mass_infeas(iter) > mass_infeas_stop_tol ) && iter < max_iter)

	%calculate the norm of the initial infeasibility
	InInf = norm(YD*v-mu);
	% call pdco with the options above, and the given parameter mu
	[v,lmda,z,inform,PDitns,CGitns,time] = ...
	    pdco(@(x)phiP(x,d),pdMat,mu,bl,bu,d1,d2,options,v,lmda,z,xsize,zsize);

	%store the objective value
	objectives =[objectives, phiP(v,d)];
	% set v equal to the solution at each iteration and store them
	v_vecs    = [v_vecs v];
	lmda_vecs = [lmda_vecs lmda];

	% output some shiz
	if inform < 1
		slvd = 1;
	else
		slvd = 0;
	end
	%calculate the result of the mapping
	mun = YAt*v;
	
	%diffnorm = [diffnorm, norm(mun - mu,inf)];
	mass_infeas = [mass_infeas, norm(Y*Ak*v,inf)];	
	mass_action_infeas = [mass_action_infeas, norm(Y'*lmda - log(v),inf)];
	fprintf(out_str,iter,InInf,objectives(iter),mass_infeas(iter),slvd,inform,toc);

	% iterate mu by setting mu = (mu + Y*A'*v)/2
	% store the new parameter vector
	mu      = (mu + mun)/2;
	mu_vecs  = [mu_vecs mu];

	% define a new v that's feasible for the next iteration
	v = 0.5*(v+(1./(d).*(A'*v)));

	% increase iteration
	iter = iter + 1;

end
mass_infeas = mass_infeas(2:iter);
iter = iter -1;

end

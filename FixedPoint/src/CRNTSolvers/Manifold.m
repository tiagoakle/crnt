%This script tries to solve for the log concentrations 
%of a steady state by using a gradient descent type method.


%Set the number of complexes to n and of species to m
%{
m = 5;
n = 50;
r = 3;

cd ../Generator
Y = YGenerator(m,n,r);
A = AkGenerator(n,0.2,1);
cd ../CRNTSolvers
%}

%Requires Y,Ak.
[m,n] = size(Y);

%XXX
%{
%Form an inhomogeneous term
gamm = randn(n,1);
iota = Y*A*gamm;
%}

%XXX
d        = - diag(Ak);
At       = A + diag(d);
eta      = abs(randn(n,1));
s        = 0.1*ones(n,1);

etaPlus  = eta + s;
etaMinus = eta + diag(1./d)*At*s;
iota     = Y*(At*etaPlus-d.*etaMinus);


eta   = ones(m,1); %initial log concentrations
rho   = exp(Y'*eta);

alpha  = 0.1;
lambda = 1e-8;
assert(size(A,1) == n);
iter = 0;

norms = [];
done = false;
while ~done
    %Find the tangent space to the manifold of complex fluxes.
    R   = diag(rho);
    rhs = iota-Y*A*rho; 
    nu  = (Y*A*R*Y'+lambda*eye(m))\(rhs);
    %Backtracking search
   % alpha = 1;
   % bt  = true;
   % b_trials = 0; 
   % while ~bt
   %     etan = eta + alpha*nu;
   %     
   % end
    eta = eta + alpha*nu; %maybe smaller step
    rho = exp(Y'*eta);
    nor = norm(Y*A*rho-iota);
    fprintf('\n Norm at iteration %i : %d, norm nu %d, norm eta %d, norm rho %d ', iter,nor,norm(nu),norm(eta),norm(rho));
    iter = iter + 1;
    norms = [norms,nor];
    if(iter > 1000) 
        done = true;
    end
end
plot(log(norms));

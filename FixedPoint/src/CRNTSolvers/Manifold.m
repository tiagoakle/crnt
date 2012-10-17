%This script tries to solve for the log concentrations 
%of a steady state by using a gradient descent type method.


%Set the number of complexes to n and of species to m

m = 500;
n = 5000;
r = 30;

cd ../Generator
Y = YGenerator(m,n,r);
A = AkGenerator(n,0.3,3);
cd ../CRNTSolvers


%Requires Y,Ak.
[m,n] = size(Y);


%Form an inhomogeneous term
gamm = randn(n,1);
iota = Y*A*gamm;


d        = - diag(A);
At       = A + diag(d);
eta      = abs(randn(n,1));
s        = 0.1*ones(n,1);

etaPlus  = eta + s;
etaMinus = eta + diag(1./d)*At*s;
iota     = Y*(At*etaPlus-d.*etaMinus);
n_iota   = norm(iota);

eta   = ones(m,1); %initial log concentrations
rho   = exp(Y'*eta);

alpha  = 1;
delta = 1e-8;
%delta = 0;
assert(size(A,1) == n);
iter = 0;

norms = [];
done = false;
while ~done
    %Find the tangent space to the manifold of complex fluxes.
    RHO   = diag(rho);
    rhs = iota-Y*A*rho; 
    %nu  = (Y*A*RHO*Y'+delta*eye(m))\(rhs);
    [Q,R,E] = qr(Y*A*RHO*Y');
    R(m,m) = 1; 
    nu   = R\(Q'*rhs);
    nu   = E*nu;
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
    if(nor/n_iota<1e-10) 
        done = true;
    end
end
fprintf('\n');
plot(log10(norms));

clear all
addpath ./Generator
addpath ./Generator/GrTheory
addpath ./PDCO
addpath ./CRNTSolvers

%Generate a random single linkage class and try 
%to converge to the fixed point of the mapping
m = 500;
n = 50;
r = 10;

%Set the seed to be able to reproduce the experiment
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',2));
%Generate a random set of n complexes on n species, each with at most 3 species
Y = YGenerator(m,n,r);
%Generate the strongly connected graph with a single linkage class`
Ak = AkGenerator(n,0.2,3);

%Extract the matrices to generate an inhomogeneous term in the
%range of YAk. 
d  = diag(Ak);
At = Ak - diag(d);


% We will have to find etaPlus and etaMinus anyway so it makes
% sense to generate b = YAtetaPlus - YDetaMinus, and generate etaPlus and etaMinus 
% strictly positive and random.

etaPlus  = zeros(n,1);
etaMinus = zeros(n,1);


%Call the solver
[iter,v_vecs,lmda_vecs,mass_infeas,mass_action_infeas]=...
		SolverFpIterationHomPDCO(Y,Ak,1.e-5,200,0.5);

fprintf('Final Mass Balance infeasibility %d\n Final Mass Action Infeasibility %d\n',mass_infeas(iter),mass_action_infeas(iter))

% plot the convergence of mu to [R F]*v
h = semilogy(1:iter, mass_infeas(1:iter),'b.-');
title('Mass balance infeasibility vs iteration count for multiple linkage network')
xlabel('Iteration')
ylabel('log ||YA_k v||_\infty')
saveas(h,'InfeasibilityVsIterationMultiple.pdf','pdf')




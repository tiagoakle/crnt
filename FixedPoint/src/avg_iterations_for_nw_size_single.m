clear all
addpath ./Generator
addpath ./Generator/GrTheory
addpath ./PDCO
addpath ./CRNTSolvers

%Calculates average iteration count vs network size
%Set the 50 network sizes to sove for
Sizes = [100:200:5000];
%Sizes = 100;
%The fraction of species is maintained constant
Sp_fraction = 0.1;
%The number of maximum active species per complex is maintained constant
r = 3;

mass_infeas_stop_tol = 1e-6;

SamplesPerNetworkSize = 20;

LinkageClasses = 1;

iter_step = 0.5;

max_iter = 500;

%Matrix to store all iteration counts per sample
Results = zeros(SamplesPerNetworkSize,length(Sizes));

%The random seed will increase by one in every generated network.
random_Seed = 0;


for j = 1:SamplesPerNetworkSize
	tic
	for k = 1:length(Sizes)
		%Set the seeds sequentially so the experiment can be reproduced
		RandStream.setDefaultStream(RandStream('mt19937ar','seed',random_Seed));
		random_Seed = random_Seed + 1;
		%Generate a random network of the given size	
		Y = YGenerator(floor(Sizes(k)*Sp_fraction),Sizes(k),r);
		Ak = AkGenerator(Sizes(k),0.2,LinkageClasses);
		%Call the fixed point iteration
		[iter,v,lmda,mi,mai]=SolverFpIterationHomPDCO(Y,Ak,mass_infeas_stop_tol,max_iter,iter_step);
		Results(j,k)= iter;
	end
	%Store a progress indicator
	fid = fopen('../Results/ProgressSingleLinkage.txt','a+');
	fprintf(fid,'Networks sampled %i, time up until last write %.2f seconds \n',k,toc);
	fclose(fid);
	%Save Results Up until now
	save './Results/NetworkSizeVsItCountSingleLink.mat' Results
end


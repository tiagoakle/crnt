%This function generates a matrix Y which given a single linkage network, would produce a 
%mass conserving network 1^TYA_k = 0. For a weight vector e = 1.
%Use Y = YGenerator(m,n,f), m and n are the dimension of the matrix, f is the maximum number of species
%per complex 
function [Y]=YGenerator(m,n,f)
	if(f>=m)
		fprintf('Error: Maximum number of species per complx can not be larger than number of species\n');
		return
	end
	Y = sparse(m,n);
	for j = 1:n
		%Choose number of non zero entries
		k = randi(f);
		%Pick values for the entries
		c = abs(rand(k));
		%pick the sparsity pattern
		p = randperm(m);
		for i=1:k
			Y(p(i),j) = c(i);
		end
		%Normalize the matrix
		Y(:,j) = Y(:,j)/norm(Y(:,j),1);
	  end	
end


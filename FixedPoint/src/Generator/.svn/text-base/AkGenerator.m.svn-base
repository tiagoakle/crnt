%This function generates a matrix Ak which has a backbone structure for a
%single linkage class. The rate-constants are chosen arbitrarily from
%Use Ak = AkGenerator(n,p), n is the number of complexes and p is the 
%probability in the Erdos-Renyi random graph generator. This also uses 
%grDecOrd.m from the GrTheory toolbox available at http://webscripts.softpedia.com/script/Scientific-Engineering-Ruby/Mathematics/matlab-grTheory-33957.html
function Ak = AkGenerator(n,p,l)
    addpath .\GrTheory\
    addpath ./GrTheory/
    classSize = randi([3,ceil(2*n/l)],l,1);
    deficit = n - sum(classSize);
    classSize = classSize + ceil(deficit/l);
    [y i] = max(classSize); 
    classSize(i) = classSize(i) - (sum(classSize)-n);  
    Ak = []; Gt=[];
    for t = 1:l    
        Dec = eye(2);
        GzeroCount = 1;
        while(size(Dec,2) > 1 || GzeroCount ~=0)
            G = rand(classSize(t),classSize(t)) < p;
            G = G - diag(diag(G));
            G = G + eye(classSize(t));
            [Ci Cj] = find(G>0);
            if (size([Ci Cj],2) ==2)
            	Dec = grDecOrd([Ci Cj]);
           	Dz = G*ones(classSize(t),1);
            	GzeroCount = classSize(t) - nnz(Dz);
            end           
        end
        A = 10*sprand(G);
        Da = A*ones(classSize(t),1);
        classCount = size(Dec,2);
        zeroCount = classSize(t) - nnz(Da);
        A = A - diag(Da);
        Aknew = A' - diag(A*ones(classSize(t),1));
        Ak = blkdiag(Ak,Aknew);
    end
%    figure(1)
%    imagesc(sign(Ak))
%    drawnow; pause(0.1)
end

function [ K ] = assembleK( k,B )
%Function to assemble the stiffness matrix.
nrows = size(B,1);
ncols = size(B,2);

sizeofK = max(B(:));
K = zeros(sizeofK);
%Begin Assembly of K using the connectivity matrix
for P=1:sizeofK
    for Q=1:sizeofK
        if(P~=Q)
            for M=1:nrows
                for N=1:ncols
                    for j=1:ncols
                        if(B(M,j)==P && B(M,N)==Q)
                            K(P,Q) = k(j,N,M);
                            K(Q,P) = k(N,j,M);                       
                        end
                    end
                end
            end
        end
        if(P==Q)
            K(P,P)=0;
            for M=1:nrows
                for N=1:ncols
                    if(P==B(M,N))
                        K(P,P) = K(P,P)+k(N,N,M);
                    end
                end
            end
        end
    end
end

end
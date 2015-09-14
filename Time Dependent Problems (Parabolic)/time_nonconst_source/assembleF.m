function [ F11 ] = assembleF(F,B,n)
%Function to assemble source vector.

%Assembled Matrix F11
nrows = size(B,1);
ncols = size(B,2);

sizeofF11 = n+1;
F11 = sym(zeros(sizeofF11,1));
%Begin Assembly using the connectivity matrix
for i=1:sizeofF11
    F11(i)=0;
    for M=1:nrows
        for N=1:ncols
            if(i==B(M,N))
                F11(i) = F11(i)+F(M,N);
            end
        end
    end
end


end


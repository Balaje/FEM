function [ y ] = solvesys( K,M,F,t_s,~,del_t,alpha,y,n )
%Function to solve time march scheme.

for z=1:n
    t_s = ones(n,1);
    for t=1:(1/del_t)-1
        t_s1 = K \ (M*t_s + del_t*(alpha*F + (1-alpha)*F) );
        y(t,z) = t_s1(z,1);
        t_s = t_s1;
    end
end

end
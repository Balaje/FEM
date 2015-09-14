function [ y ] = solvesys1( K,M,F,u,del_t,n,y,a1,h,b1 )
%Function to solve time march scheme.

F = matlabFunction(F);
for z=1:n
    t_s = zeros(n,1);
    l=1;
    for j=a1+h:h:b1-h
        t_s(l) = u(j);
        l = l+1;
    end
    dd=0;
    for t=0:(1/del_t)-1
        t_s1 = K \ (M*t_s + feval(F,t));
        y(t+1,z) = t_s1(z,1);
        t_s = t_s1;
        dd = dd+del_t;
    end
end

end
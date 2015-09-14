function [ sol ] = exactsol(n,alpha,px,qx,fx,del_t)
%Function to compute exact solution
    function [c,f,s] = pdefun(x,t,u,DuDx)
        c=qx/px;
        f=DuDx;
        s=fx/px;
    end

    function u0 = icfun(x)
        u0 = 1;
    end

    function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
        pl = ul;
        ql = 0;
        pr = 0;
        qr = 1;
    end

    x = linspace(0,1,n);
    t = linspace(0,1,1/del_t);
    
    m=0;
    sol = pdepe(m, @pdefun, @icfun, @bcfun, x, t);
end


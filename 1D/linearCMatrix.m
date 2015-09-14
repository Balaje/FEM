%3. Linear interpolation functions using connectivity matrix

clear;
clc;
syms x;

%Input functions p,q,r,f;
p = @(x) x^2;
q = @(x) x;
r = @(x) 1;
f = @(x) x;

%Domain.
a1 = 1;
b1 = exp(1);

%Boundary conditions (nan) if not specified.
ua = exp(1);
ub = 1;
Qa = nan;
Qb = nan;

%Exact Solution.
syms y(x);
Dy = diff(y);
a(x) = dsolve( p(x)*diff(y,x,x)+q(x)*diff(y,x)+r(x)*y-f(x)==0, y(1)==exp(1), y(exp(1))==1);

%for l=3:10
%Number of elements.
n = 4;
%--Algorithm Begins
h = (b1-a1)/n; %Interval size

% Connectivity Matrix for linear elements B
B = zeros(n,2);
num = 1;
for i=1:n
    for j=1:2
        B(i,j) = num;
        if(j==1)
            num = num+1;
        end
    end
end

psi1 = @(x) 1-x/h; %Interpolation functions
psi2 = @(x) x/h;

%CONSTRUCTING ELEMENT EQUATIONS
%stiffness matrix k
P = p(x+a1);
Q = q(x+a1);
R = r(x+a1);

k11 = int( (Q/P)*psi1(x)*diff(psi1(x),x) + (R/P)*psi1(x)*psi1(x) - diff(psi1(x),x)*diff(psi1(x),x),x,0,h);
k12 = int( (Q/P)*psi1(x)*diff(psi2(x),x) + (R/P)*psi1(x)*psi2(x) - diff(psi1(x),x)*diff(psi2(x),x),x,0,h);
k21 = int( (Q/P)*psi2(x)*diff(psi1(x),x) + (R/P)*psi2(x)*psi1(x) - diff(psi2(x),x)*diff(psi1(x),x),x,0,h);
k22 = int( (Q/P)*psi2(x)*diff(psi2(x),x) + (R/P)*psi2(x)*psi2(x) - diff(psi2(x),x)*diff(psi2(x),x),x,0,h);

k = [k11,k12;k21,k22];
j=2;
for i=a1+h:h:b1-h
    P = p(x+i);
    Q = q(x+i);
    R = r(x+i);
    
    k11 = int( (Q/P)*psi1(x)*diff(psi1(x),x) + (R/P)*psi1(x)*psi1(x) - diff(psi1(x),x)*diff(psi1(x),x),x,0,h);
    k12 = int( (Q/P)*psi1(x)*diff(psi2(x),x) + (R/P)*psi1(x)*psi2(x) - diff(psi1(x),x)*diff(psi2(x),x),x,0,h);
    k21 = int( (Q/P)*psi2(x)*diff(psi1(x),x) + (R/P)*psi2(x)*psi1(x) - diff(psi2(x),x)*diff(psi1(x),x),x,0,h);
    k22 = int( (Q/P)*psi2(x)*diff(psi2(x),x) + (R/P)*psi2(x)*psi2(x) - diff(psi2(x),x)*diff(psi2(x),x),x,0,h);

    k1 = [k11,k12;k21,k22];
    k(:,:,j) = k1;
    j=j+1;
end

%Source Vector
F = zeros(n,2);
j=1;
for i=a1:h:b1-h
    F1 = f(x+i)/p(x+i);
    F(j,1) = int(psi1(x)*F1,0,h);
    j=j+1;
end
j=1;
for i=a1:h:b1-h
    F1 = f(x+i)/p(x+i);
    F(j,2) = int(psi2(x)*F1,0,h);
    j=j+1;
end


%Secondary boundary conditions vector.
Qx = zeros(n+1,1);
if(isnan(Qa))
    Qa=0;
else
    Qx(1)=Qa;
end
if(isnan(Qb))
    Qb=0;
else
    Qx(n+1)=-Qb;
end

%ASSEMBLING THE FINITE ELEMENT MODEL.   
%Assembled K matrix:
nrows = size(B,1);
ncols = size(B,2);

sizeofK = max(B(:));
K = zeros(sizeofK);
%Begin Assembly using the connectivity matrix
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

%Assembled Matrix F11
sizeofF11 = n+1;
F11 = zeros(sizeofF11,1);
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

%Get the rhs matrix.
rhs = Qx+F11;

%Adjusting the assembly according to the boundary conditions
N = n;
if(~isnan(ua))
    %If the node value is 0, remove the row.
    if(ua==0)
        K(:,1) = [];
        K(1,:) = [];
        rhs(1,:) = [];
        N=N-1;
    end
    %If the node value is given, set the value in the system.
    if(ua~=0)
        K(1,1)=1;
        for i=2:N
            K(1,i) = 0;
        end
        rhs(1) = ua;
    end
end
if(~isnan(ub))
    %If the node value is 0, remove the row.
    if(ub==0)
        K(:,N+1) = [];
        K(N+1,:) = [];
        rhs(N+1,:) = [];
        N=N-1;
    end
    %If the node value is given, set the value in the system.
    if(ub~=0)
        K(N+1,N+1)=1;
        for i=1:N
           K(N+1,i) = 0;
        end
        rhs(N+1) = ub;
    end
end

%Solve the system.
res = linsolve(K,rhs);
%Construct a new matrix that includes all the node values
RES = zeros(numel(res)+2,1);
if(ub==0 && ua==0)
    RES(1) = ua;
    RES(numel(RES))=0;
    for j=1:numel(res)
        RES(j+1) = res(j);
    end
end
if(ua==0 && ub~=0)
    for j=1:numel(res)
        RES(j+1) = res(j);
    end
    RES(numel(RES))=[];
end
if(ua~=0 && ub==0)
    for j=1:numel(res)
        RES(j) = res(j);
    end
    RES(numel(RES)-1)=[];
end
if(ua~=0 && ub~=0)
    RES = zeros(numel(res),1);
    for j=1:numel(res)
        RES(j) = res(j);
    end
end

display(K);
display(rhs);
fprintf('\nThe solution of the differential equation: ');
display(vpa(RES,10));

%Compare exact and approximate solutions
p1 = a1:h:b1;
r1 = a(p1);
head = strcat('Uniform mesh for ',num2str(n),' elements');

plot(p1,r1,'--',p1,RES,'r-o'),legend('Exact','Approx'),title(head);
xlabel('x');
ylabel('y(x)');
grid on
fprintf('\n The exact solution of the differential equation: ');
display(vpa(r1)');
%pause(2);
%end
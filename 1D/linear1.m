%1. Linear interpolation functions

clear;
clc;
syms x;

%Input functions p,q,r,f;
p = @(x) 1;
q = @(x) 0;
r = @(x) 1;
f = @(x) 0;

%Domain.
a1 = 0;
b1 = 2*pi;

%Boundary conditions (nan) if not specified.
ua = 0; %Dirichlet
ub = nan;
Qa = nan; %Neumann type
Qb = 1;

%Exact Solution.
syms y(x);
Dy = diff(y);
a(x) = dsolve( p(x)*diff(y,x,x)+q(x)*diff(y,x)+r(x)*y-f(x)==0, y(0)==0, Dy(2*pi)==1);

%for l=1:6
%Number of elements.
n = 15;

%Algorithm Begins.
h = (b1-a1)/n; %Interval size

psi1 = @(x) 1-x/h; %Interpolation functions
psi2 = @(x) x/h;

%CONSTRUCTING ELEMENT EQUATIONS
k11 = zeros(1,n);
k12 = zeros(1,n);
k21 = zeros(1,n);
k22 = zeros(1,n);
%Stiffness matrix
j=1;
for i=a1:h:b1-h
    P = p(x+i);
    Q = q(x+i);
    R = r(x+i);
    k11(j) = int( (Q/P)*psi1(x)*diff(psi1(x),x) + (R/P)*psi1(x)*psi1(x) - diff(psi1(x),x)*diff(psi1(x),x),x,0,h);
    k12(j) = int( (Q/P)*psi1(x)*diff(psi2(x),x) + (R/P)*psi1(x)*psi2(x) - diff(psi1(x),x)*diff(psi2(x),x),x,0,h);
    k21(j) = int( (Q/P)*psi2(x)*diff(psi1(x),x) + (R/P)*psi2(x)*psi1(x) - diff(psi2(x),x)*diff(psi1(x),x),x,0,h);
    k22(j) = int( (Q/P)*psi2(x)*diff(psi2(x),x) + (R/P)*psi2(x)*psi2(x) - diff(psi2(x),x)*diff(psi2(x),x),x,0,h);
    j=j+1;
end
%source vector
f1 = zeros(1,n);
f2 = zeros(1,n);
j=1;
for i=a1:h:b1-h
    F1 = f(x+i)/p(x+i);
    f1(j) = int(psi1(x)*F1,0,h);
    f2(j) = int(psi2(x)*F1,0,h);
    j=j+1;
end
%Secondary boundary conditions vector.
Qx = zeros(n+1,1);
qa = 0;
qb = 0;
if(isnan(Qa))
    qa=0;
else
    Qx(1)=Qa;
end
if(isnan(Qb))
    qb=0;
else
    Qx(n+1)=-Qb;
end

%ASSEMBLING THE FINITE ELEMENT MODEL.
%Matrix F
F = zeros(n+1,1);
F(1) = f1(1);
F(n+1) = f2(n);
for i=2:n
    F(i) = f2(i-1)+f1(i);
end

%Assembled K matrix;
m1 = diag(k21,-1);
m3 = diag(k12,1);
mid = zeros(1,n+1);
mid(1) = k11(1);
mid(n+1) = k22(n);
for i=2:n
    mid(i) = k22(i-1)+k11(i);
end
m2 = diag(mid,0);
K = m1+m2+m3;

%Get the rhs matrix.
rhs = Qx+F;

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
display(vpa(r1,7)');
%pause(3)
%end
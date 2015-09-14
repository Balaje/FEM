%Program to solve an ordinary differential equation using finite element
%method

%Takes a differential equation of the form -( d/dx (a*du/dx) )+cu-f = 0 
%subject to the boundary conditions in the domain [a1,b1]

%Instructions
%1. Input f(x), a(x) and c(x) and the no of subintervals.
%Program assumes uniform mesh.

clear;
clc;
syms x;
%Input the parameters describing the differential equation
fprintf('\tEnter in terms of function handles: Eg. @(x) x^2:\n');
a = input('Enter function a(x): ');
c = input('Enter function c(x): ');
f = input('Enter function f(x): ');

%Input Domain
b1 = 0;
a1 = 0;
while b1<=a1
    dom = input('\nEnter domain in the form of [a,b]: ');
    a1 = dom(1);
    b1 = dom(2);
    
    if b1<a1
        fprintf('\nInvalid Domain.');
    end
end

%The boundary conditions

fprintf('\nEnter the boundary conditions: \n If not specified enter "nan"\n');
u0 = input('Enter the first condition u(a1)= ');
Q0 = input('Enter the second condition [a*du/dx at x=a1] = ');
u1 = input('Enter the third condition u(b1)= ');
Q1 = input('Enter the fourth condition [a*du/dx at x=b1] = ');

%Input number of subintervals.
%With this calculate the approximation functions psi1 and psi2.
n = input('\nEnter the number of subintervals: ');
h = (b1-a1)/n;
psi1 = 1-x/h;
psi2 = x/h;

%calculating the Stiffness matrix k|2x2 k = [k11,k12; k21,k22]
k11 = zeros(1,n);
k12 = zeros(1,n);
k22 = zeros(1,n);
%Calculating the elements
j=1;
for i=a1:h:b1-h
    A = a(x+i);
    C = c(x+i);
    k12(j) = int(A*diff(psi1,x)*diff(psi2,x)+C*psi1*psi2,0,h);
    k11(j) = int(A*diff(psi1,x)*diff(psi1,x)+C*psi1*psi1,0,h);
    k22(j) = int(A*diff(psi2)*diff(psi2)+C*psi2*psi2,0,h);
    j=j+1;
end
%Calculating the source vector
f1 = zeros(1,n);
f2 = zeros(1,n);
j=1;
for i=a1:h:b1-h
    F1 = f(x+i);
    f1(j) = int(psi1*F1,0,h);
    f2(j) = int(psi2*F1,0,h);
    j=j+1;
end

%Assemble the elements by using the continuity of the elements
%[K][U]={F}+{Q}
%Find the matrix {F}
F = zeros(n+1,1);
F(1) = f1(1);
F(n+1) = f2(n);
for i=2:n
    F(i) = f2(i-1)+f1(i);
end
%Find the matrix [K] --> Tridiagonal Matrix.
m1 = diag(k12,-1);
m3 = diag(k12,1);

mid = zeros(1,n+1);
mid(1) = k11(1);
mid(n+1) = k22(n);
for i=2:n
    mid(i) = k22(i-1)+k11(i);
end
m2 = diag(mid,0);

K = m1+m2+m3;
%Find the matrix [Q] --> Secondary boundary conditions.
Q = zeros(n+1,1);
if(isnan(Q0))
    Q0=0;
else
    Q(1)=-Q0;
end
if(isnan(Q1))
    Q1=0;
else
    Q(n+1)=Q1;
end

%Get the rhs matrix
rhs = Q+F;
N=n;
%If the boundary condition is given..
if(~isnan(u0))
    %If the node value is 0, remove the row.
    if(u0==0)
        K(:,1) = [];
        K(1,:) = [];
        rhs(1,:) = [];
        N=N-1;
    end
    %If the node value is given, set the value in the system.
    if(u0~=0)
        K(1,1)=1;
        for i=2:N
            K(1,i) = 0;
        end
        rhs(1) = u0;
    end
    
end

if(~isnan(u1))
    %If the node value is 0, remove the row.
    if(u1==0)
        K(:,N+1) = [];
        K(N+1,:) = [];
        rhs(N+1,:) = [];
        N=N-1;
    end
    %If the node value is given, set the value in the system.
    if(u1~=0)
        K(N+1,N+1)=1;
        for i=1:N
           K(N+1,i) = 0;
        end
        rhs(N+1) = u1;
    end
end

%Solve the system.
res = linsolve(K,rhs);

RES = zeros(numel(res)+2,1);
if(u1==0 && u0==0)
    RES(1) = u0;
    RES(numel(RES))=0;
    for j=1:numel(res)
        RES(j+1) = res(j);
    end
end

if(u0==0 && u1~=0)
    for j=1:numel(res)
        RES(j+1) = res(j);
    end
    RES(numel(RES))=[];
end

if(u0~=0 && u1==0)
    for j=1:numel(res)
        RES(j) = res(j);
    end
    RES(numel(RES)-1)=[];
end

if(u0~=0 && u1~=0)
    RES = zeros(numel(res),1);
    for j=1:numel(res)
        RES(j) = res(j);
    end
end

fprintf('\nThe solution of the differential equation: ');
display(vpa(RES,10));

if(~isnan(u0) && ~isnan(u1))
    syms y(x);
    exact = zeros(N,1);
    error = zeros(N,1);
    
    bc = strcat('y(',num2str(a1),')==',num2str(u0),',y(',num2str(b1),')==',num2str(u1));
    y(x) = dsolve(-a(x)*diff(y,x,x) - diff(a(x),x)*diff(y,x)+c(x)*y-f(x)==0,bc);
    fprintf('\nExact solution corresponding to the Dirichlet boundary condition:');
    display(y(x));
    
    j=1;
    for i=a1:h:b1
        exact(j) = vpa(y(i),4);
        error(j) = abs(RES(j) - exact(j));
        j=j+1;
    end
    fprintf('\nThe exact solution at the nodes:\n');
    display(vpa(exact,10));
    fprintf('\nError in computing the solution using fem:\n');
    display(error);
    
    p = a1:h:b1;
    q = a1:0.001:b1;
    r = y(q);
    subplot(1,3,1);
    plot(q,r,'m'),legend('Exact');
    subplot(1,3,2);
    plot(p,RES,'b'),legend('Approx');
    subplot(1,3,3);
    plot(p,error,'r');
end

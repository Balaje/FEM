% Time dependent IBVP
% Non homogeneous Parabolic equations with constant coefficients and source
% term
% Homogeneous neumann boundary conditions
% Non homogeneous constant initial conditon

% Input coefficients a,b,c0,c1 of the differential equation 
% of the form in J.N.Reddy pg.316
% (Heat equation with constant coefficients and constant source term subjected
% to homogeneous neumann boundary condtions)

% solvesys() to solve the system
clc
clear
syms x;
syms t;

%Input functions
a = 1;
b = 0;
c0 = 0;
c1 = 1;
f = 0;

%Domain
a1 = 0;
b1 = 1;

%Initial and boundary conditions
u0_t = 0;
uL_t = 0; % Boundary conditions
ux_0 = 1; % Initial conditions

%Number of elements (discretize x)
n = 10;
h = (b1-a1)/n;

%Enter time step size (1/del_t must be an integer)
del_t = 0.05;

%Connectivity Matrix
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

psi1 = @(x) 1-x/h;
psi2 = @(x) x/h;

k11 = int( (a*diff(psi1(x),x)*diff(psi1(x),x) + b*diff(psi1(x),x,2)*diff(psi1(x),x,2) + c0*psi1(x)*psi1(x)) ,x,0,h );
k12 = int( (a*diff(psi1(x),x)*diff(psi2(x),x) + b*diff(psi1(x),x,2)*diff(psi2(x),x,2) + c0*psi1(x)*psi2(x)) ,x,0,h );
k21 = k12;
k22 = int( (a*diff(psi2(x),x)*diff(psi2(x),x) + b*diff(psi2(x),x,2)*diff(psi2(x),x,2) + c0*psi2(x)*psi2(x)) ,x,0,h );

m11 = int( (c1*psi1(x)*psi1(x)) ,x,0,h );
m12 = int( (c1*psi1(x)*psi2(x)) ,x,0,h );
m21 = m12;
m22 = int( (c1*psi2(x)*psi2(x)) ,x,0,h );    

Ke = [k11,k12;k21,k22];
Me = [m11,m12;m21,m22];
j=2;
for i=a1+h:h:b1-h
    k11 = int( (a*diff(psi1(x),x)*diff(psi1(x),x) + b*diff(psi1(x),x,2)*diff(psi1(x),x,2) + c0*psi1(x)*psi1(x)) ,x,0,h );
    k12 = int( (a*diff(psi1(x),x)*diff(psi2(x),x) + b*diff(psi1(x),x,2)*diff(psi2(x),x,2) + c0*psi1(x)*psi2(x)) ,x,0,h );
    k21 = k12;
    k22 = int( (a*diff(psi2(x),x)*diff(psi2(x),x) + b*diff(psi2(x),x,2)*diff(psi2(x),x,2) + c0*psi2(x)*psi2(x)) ,x,0,h );    
    
    m11 = int( (c1*psi1(x)*psi1(x)) ,x,0,h );
    m12 = int( (c1*psi1(x)*psi2(x)) ,x,0,h );
    m21 = m12;
    m22 = int( (c1*psi2(x)*psi2(x)) ,x,0,h );    

    k1 = [k11,k12;k21,k22];
    m1 = [m11,m12;m21,m22];
    Ke(:,:,j) = k1;
    Me(:,:,j) = m1;
    j=j+1;
end

%Obtain source element matrix
f1 = sym(zeros(n,2));
j=1;
for i=a1:h:b1-h
    f1(j,1) = int(psi1(x)*f,x,0,h);
    j=j+1;
end
j=1;
for i=a1:h:b1-h
    f1(j,2) = int(psi2(x)*f,x,0,h);
    j=j+1;
end

%Vary the method alpha = 0.5 --> Crank Nicolson Method
alpha = 0.5;

LEFT = Me+del_t*alpha*Ke;
RIGHT = Me-del_t*(1-alpha)*Ke;

K = assembleK(LEFT,B);
M = assembleK(RIGHT,B);
F11 = assembleF(f1,B,n);

%Applying the homogeneous dirichlet boundary condition
N = n;
K(:,1) = [];
M(:,1) = [];
K(1,:) = [];
M(1,:) = [];
F11(1,:) = [];
N=N-1;

%{
K(:,N+1) = [];
M(:,N+1) = [];
K(N+1,:) = [];
M(N+1,:) = [];
F11(N+1,:) = [];
N=N-1;
%}
%Time march system.
%function outputs values of u(1,t);
y = zeros(ceil((1/del_t)-1),n);
t_s = ones(n,1);
t_s1 = zeros(n,1);
z1 = solvesys(K,M,F11,t_s,t_s1,del_t,alpha,y,n);

%sol = exactsol(n+1,alpha,a,c1,f);
%Plot the data corresponding to the second node. 0<t<10
x1 = 0:del_t:1-del_t;
A = zeros(1,n);
for i=1:n
    A(i) = ux_0;
end
z1 = [A;z1];
sol = exactsol(n,alpha,a,c1,f,del_t);
%{
for j=1:n-1
    plot(x1,y1(:,j)','r-o',x1,sol(:,j+1),'m-.o'),legend('Approximate solution','Exact Solution');
    hold on
    grid on
end
%}
%Plot for u(1,t);
plot(x1,z1(:,n-1)','r-o',x1,sol(:,n),'m-.o'),legend('Approximate solution','Exact solution');
title('Plot for u(1,t)');
xlabel('Time t');
ylabel('u(1,t)');
hold on
grid on;

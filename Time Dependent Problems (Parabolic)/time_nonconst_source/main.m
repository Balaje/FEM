% Time dependent IBVP
% Non homogeneous Parabolic equations with constant coefficients and non-constant source
% term
% Homogeneous Dirichlet boundary conditions
% Non homogeneous non constant initial conditon

% Input coefficients a,b,c0,c1 of the differential equation 
% of the form in J.N.Reddy pg.316
% (Heat equation with constant coefficients and non-constant source term subjected
% to homogeneous dirichelt boundary condtions)

% solvesys1() to solve the system
clc
clear
syms x;
syms t;

%Input functions
a = 1;
b = 0;
c0 = 0;
c1 = 1;
f = @(x,t) t*sin(x);

%Domain
a1 = 0;
b1 = pi;

%Initial and boundary conditions
u0_t = 0;
uL_t = 0; % Boundary conditions (Dirichlet)
ux_0 = @(x) sin(x); % Initial conditions

%Number of elements (discretize x)
n = 10;
h = (b1-a1)/n;

%Enter time step size t=[0,1] (1/del_t must be an integer)
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
    f1(j,1) = int(psi1(x)*f(x+i,t),x,0,h);
    j=j+1;
end
j=1;
for i=a1:h:b1-h
    f1(j,2) = int(psi2(x)*f(x+i,t),x,0,h);
    j=j+1;
end

%Vary the method alpha = 0.5 --> Crank Nicolson Method
alpha = 0.5;

LEFT = Me+del_t*alpha*Ke;
RIGHT = Me-del_t*(1-alpha)*Ke;

%Time march scheme
syms s;
RIGHT1 = del_t*( alpha*feval(matlabFunction(f1),(s+1)*del_t) + (1-alpha)*feval(matlabFunction(f1),s*del_t));

K = assembleK(LEFT,B);
M = assembleK(RIGHT,B);
F = assembleF(RIGHT1,B,n);

%Applying the Dirichlet boundary conditions
N = n;
K(:,1) = [];
M(:,1) = [];
K(1,:) = [];
M(1,:) = [];
F(1,:) = [];
N=N-1;

K(:,N+1) = [];
M(:,N+1) = [];
K(N+1,:) = [];
M(N+1,:) = [];
F(N+1,:) = [];
N=N-1;

%Compute the solutions
y = sym(zeros((1/del_t)-1,n-1));
sol = solvesys1(K,M,F,ux_0,del_t,n-1,y,a1,h,b1);

%Comparing it with the exact solution (IF KNOWN);
x1 = 0:del_t:1-del_t;
syms x2;
%Plot values for u(1,t)
g = @(x2) (2*exp(-x2)+x2-1)*sin(b1-h);
plot(x1,sol(:,n-1)','b-x',x1,g(x1),'r-+'),legend('Approximate solution','Exact Solution');
title('Plot values for u(1,t)');
xlabel('Time t');
ylabel('u(1,t)');
hold on
grid on
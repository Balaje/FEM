% Code to implement Example 8.3.1, J.N.Reddy
% u_xx + u_yy = -f0 in a square region (-A,A) x (-A,A)
% u = 0 on the boundary

% NOTE: Comment out matrices not used 
% Code: Running for T2 Mesh.

clear;
clc;
syms x;
syms y;

%Input matrices,
a11 = 1;
a12 = 0;
a21 = 0;
a22 = 1;
a00 = 0;
f = 1;

%Input no of elements,
%n=4;
%n=9;
n=16; %Run for 16 elements

%Input connectivity matrix. From this get the number of nodes(Global)
p = 3; %Should be equal to 3 in this example(Triangular element).
q = n; %Should be equal to number of elements.

%B = [1,2,3; 5,3,2; 2,4,5; 3,5,6];
%B = [1,2,3; 5,3,2; 2,4,5; 8,5,4; 4,7,8; 3,5,6; 9,6,5; 5,8,9; 6,9,10];
B = [1,2,3; 5,3,2; 2,4,5; 8,5,4; 4,7,8; 12,8,7; 7,11,12; 3,5,6; 9,6,5; 5,8,9; 13,9,8; 8,12,13; 6,9,10; 14,10,9; 9,13,14; 10,14,15];

A = [1,0,0; 1,1/sqrt(n),0; 1,1/sqrt(n),1/sqrt(n)];

PSI = [1,x,y]*inv(A);
%Element matrix K
k = zeros(3,3);
for i=1:3
    for j=1:3
        F = @(x,y) diff(PSI(i),x).*(a11*diff(PSI(j),x)+a12*diff(PSI(j),y))+diff(PSI(i),y)*(a21*diff(PSI(j),x)+a22*diff(PSI(j),y)+a00*PSI(i)*PSI(j));
        k(i,j)=int(int(F(x,y),y,0,x),x,0,1/sqrt(n));
    end
end
%Element matrix F
F1 = zeros(3,1);
for i=1:3
    F1(i) = int(int(PSI(i)*f,y,0,x),x,0,1/sqrt(n));
end

%Assembled K matrix:
nrows = size(B,1);
ncols = size(B,2);

sizeofK = max(B(:));
K = zeros(sizeofK);
%Begin Assembly using the connectivity matrix
for P=1:sizeofK
    for Q=1:sizeofK
        K(P,Q) = 0;
        K(Q,P) = 0;
        for M=1:nrows
            for N=1:ncols
                for j=1:ncols
                    if((B(M,j)==P && B(M,N)==Q))
                        K(P,Q) = K(P,Q)+k(j,N);
                        K(Q,P) = K(Q,P)+k(N,j);
                    end
                end
            end
        end
        if(P==Q)
            K(P,P)=0;
            for M=1:nrows
                for N=1:ncols
                    if(P==B(M,N))
                        K(P,P) = K(P,P)+k(N,N);
                    end
                end
            end
        end
    end
end

sizeofF11 = max(B(:));
F11 = zeros(sizeofF11,1);
%Begin Assembly using the connectivity matrix
for i=1:sizeofF11
    F11(i)=0;
    for M=1:nrows
        for N=1:ncols
            if(i==B(M,N))
                F11(i) = F11(i)+F1(N);
            end
        end
    end
end

%Given u4,u5,u6=0 [remove corresponding rows for finer meshes]
for i=sizeofK:-1:sizeofK-(sqrt(n))
    K(i,:)=[];
    K(:,i)=[];
    F11(i,:)=[];
end

res=linsolve(K,F11);
fprintf('Solution at nodes (in order). Values at the right hand side of the triangle is zero:\n');
disp(res);

%Plot the mesh
x1 = 0:1/sqrt(n):1;
%{
z = @(x) x;
y1 = z(x1);
plot(x1,y1,'o--');
hold on
%}
for i=0:1/sqrt(n):1
    z=@(x) x-i;
    y1 = z(x1);
    plot(x1,y1,'o--');
    hold on;
    ylim([0,1]);
end

sum = 0;
j=1;
for i=0:1/sqrt(n):1
    x1=0:1/sqrt(n):sum;
    x2=i*ones(size(x1));
    x3=sum:1/sqrt(n):1;
    x4=i*ones(size(x3));
    
    plot(x2,x1,'o-r');
    plot(x3,x4,'o-r');
    set(gca,'XTick',[])
    
    sum = sum + 1/sqrt(n);
    hold on
end

k=1;
for i=0:1/sqrt(n):1
    for j=0:1/sqrt(n):i
        str = [' ',' ',' ',num2str(k)];
        text(i,j,str);
        k=k+1;
    end
end
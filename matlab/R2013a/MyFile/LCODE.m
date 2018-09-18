function dx=LCODE(t,x)

%Irreducible Matrix
%n=500;
%L=rand(n)*10;
%for i=1:n
%L(i,i)=-(sum(L(i,:))-L(i,i));
%end
%L=round(L);
L=[   -18     3     1     7     7
     1   -11     3     1     4
     2     6   -18     1     8
     2     0     6   -14     5
    10     2     1     5   -18];

L_dim=length(L);
c=100;
t
E=zeros(L_dim,L_dim);
i0=maxlocation(L);
E(i0,i0)=1;

dx_length=L_dim*3+4;
lh=dx_length;
dx=zeros(dx_length,1);

% s(t)
s = x(1:3);
dx(1:3)= f(s);
 ex=zeros(3,1);
    
% x_i(t)
for i=1:L_dim
      
    xi = x(3*i+(1:3));
    ic=0;
    %p=diag(rand(1,3));
    for k=1:L_dim
    ex(1)=x(3*k+1)-x(1);
    ex(2)=x(3*k+2)-x(2);
    ex(3)=x(3*k+3)-x(3);
        ic=ic+14*ex'*ex;
    end
    
    ux=f(xi)-c*E(i,i)*x(lh)*(h(xi)-h(s));  
    
    for j=1:L_dim        
        xj=h(x(3*j+(1:3)));        
        ux=ux+c*L(i,j)*x(lh)*xj;
    end
    
    dx(3*i+(1:3))=ux;
    dx(lh)=ic;
end
end

function  fx=f(x)
%Lorenz's oscillators
fx = zeros(3,1);
fx(1)= 10*(x(2) - x(1));
fx(2)= 28*x(1)- x(1)*x(3) - x(2);
fx(3)= x(1)*x(2) -8/3*x(3);
end

function  hx=h(x)
hx=x;
end
function w=maxlocation(A)
m=length(A);
B=zeros(1,m);
for j=1:m
    B(j)=sum(A(:,j))-A(j,j);
end
b=B(1);n=1;
for i=1:m
    for j=i+1:m
        if B(j)>=B(i)&&B(j)>=b
            n=j;
            b=B(j);
        else
            q=1;
        end
    end
end
w=n;
end
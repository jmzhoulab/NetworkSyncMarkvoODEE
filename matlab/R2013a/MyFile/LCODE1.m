function dx=LCODE1(t,x)

%Irreducible Matrix
n=500;
L=rand(n)*10;
for i=1:n
L(i,i)=-(sum(L(i,:))-L(i,i));
end
L=round(L);
L_dim=length(L);
c=100;
t
E=zeros(L_dim,L_dim);
i0=maxlocation(L);
E(i0,i0)=1;

dx_length=L_dim*3+3;
lh=dx_length;
dx=zeros(dx_length,1);

% s(t)
s = x(1:3);
dx(1:3)= f(s);

% x_i(t)
for i=1:L_dim
      
    xi = x(3*i+(1:3));
    
    ux=f(xi)-c*E(i,i)*(h(xi)-h(s));  
    
    for j=1:L_dim        
        xj=h(x(3*j+(1:3)));        
        ux=ux+c*L(i,j)*xj;
    end
    
    dx(3*i+(1:3))=ux;
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
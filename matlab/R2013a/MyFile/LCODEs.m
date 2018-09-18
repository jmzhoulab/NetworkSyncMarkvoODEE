function dx=LCODEs(t,x)

%Irreducible Matrix
n=5;
t
L=rand(n)*10;
for i=1:n
L(i,i)=-(sum(L(i,:))-L(i,i));
end
L=round(L);
L_dim=length(L);

E=zeros(L_dim,L_dim);
E(1,1)=1;
E(2,2)=1;
Gamma = eye(3);
C=5;

dx_length=L_dim*3+3+1;
dx=zeros(dx_length,1);

% s(t)
s = x(1:3);
dx(1:3)= f(s);
ct = 0;
% x_i(t)
for i=1:L_dim
      
    xi = x(3*i+(1:3));    

    ux=f(xi)-x(dx_length)*E(i,i)*Gamma*(h(xi)-h(s));  
    
    for j=1:L_dim        
        xj=h(x(3*j+(1:3)));        
        ux=ux+x(dx_length)*L(i,j)*Gamma*xj;
    end    
    
    ct = ct + 0.01 * (xi - s)'*(xi - s);
	
    dx(3*i+(1:3))=ux;
end

dx(dx_length) = ct;
end

function  fx=f(x)
%Lorenz's oscillators
fx = zeros(3,1);
fx(1)= 10*(x(2) - x(1));
fx(2)= 28*x(1)- x(1)*x(3) - x(2);
fx(3)= x(1)*x(2) -8/3*x(3);
end

function  hx=h(x)
hx=5*x+sin(x);
end
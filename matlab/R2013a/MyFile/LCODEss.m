function dx=LCODEss(t,x)

%Irreducible Matrix
L=[ -3 1 1 1 0;...
    1 -3 1 1 0;...
    1 1 -3 0 1;...
    1 1 0 -3 1;...
    0 0 1 1 -2];
L_dim = 5;

E=zeros(L_dim,L_dim);
E(1,1)=1;
E(2,2)=1;
Gamma = eye(3);
C=5;

dx_length=L_dim*3+3;
dx=zeros(dx_length,1);

% s(t)
s = x(1:3);
dx(1:3)= f(s);

% x_i(t)
for i=1:L_dim
      
    xi = x(3*i+(1:3));    

    ux=f(xi)-C*E(i,i)*Gamma*(h(xi)-h(s));  
    
    for j=1:L_dim        
        xj=h(x(3*j+(1:3)));        
        ux=ux+C*L(i,j)*Gamma*xj;
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
hx=5*x+sin(x);
end
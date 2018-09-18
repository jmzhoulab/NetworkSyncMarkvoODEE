function dx=wf(t,x)
dx=zeros(2,1);
dx(1)=x(2);
dx(2)=-2*x(2)-x(1);
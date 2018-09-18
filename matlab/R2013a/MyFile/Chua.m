function dx= Chua(t,x)
dx=zeros(3,1);
p=10; q=14.87;
dx(1)=p*(x(2)-g(x(1)));
dx(2)=x(1)-x(2)+x(3);
dx(3)=-q*x(2);
end

function gx=g(x)
gx=0.32*x-0.295*(abs(x(1)+1)-abs(x(1)-1));
end
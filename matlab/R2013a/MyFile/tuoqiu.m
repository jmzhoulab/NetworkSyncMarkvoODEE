function tuoqiu=tuoqiu(a,b,c)
s=linspace(-pi/2,pi/2,100);
t=linspace(0,2*pi,100);
[s,t]=meshgrid(s,t);
x=a.*cos(s).*cos(t);         
y=b.*cos(s).*sin(t);   
z=c.*sin(s);
mesh(x,y,z);
xlabel('x');
ylabel('y');
zlabel('z');
title('Õ÷«Ú√Ê æ“‚Õº');
shading  flat;
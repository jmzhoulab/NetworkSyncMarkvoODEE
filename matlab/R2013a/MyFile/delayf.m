function fx=delayf(x,y)
    p=9.78;
    a=-1.4325;
    b=-0.7831;
    c=0.2;
    eta=9.53;
    v=0.5;
    e=0.1636;
    fx = zeros(1,3);
    fx(1)= p*(-(1+b)*x(1)+x(2)-1/2*(a-b)*(abs(x(1)+1)-abs(x(1)-1)));   
    fx(2)= x(1)-x(2)+x(3);
    fx(3)= -eta*x(2)-e*x(3)-eta*c*sin(v*y(1));
    %{
    %f(x)µÄÑÅ¿É±È¾ØÕó
    J1=[-p*(1+a)  p   0;
        1        -1   1;
        0        -eta 0];
    J2=[-p*(1+b)  p   0;
        1        -1   1;
        0        -eta 0];
    alpha1=max([norm(J1) norm(J2)]);
    alpha2=abs(eta*c*v)
    %}
end
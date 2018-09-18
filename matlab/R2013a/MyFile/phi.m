function a=phi(t,x,y,u0,v0)
    l1=3.0608;
    l2=3.4438;
    beta=3.25;
    eta1=(norm(x-y)+l1*norm(u0-v0))/l1;
    eta2=(norm(x-y)+l2*norm(u0-v0))/l2;
    a=beta*sqrt(eta1^2*(exp(l1*t)-1)^2+eta2^2*(exp(l2*t)-1)^2);
end
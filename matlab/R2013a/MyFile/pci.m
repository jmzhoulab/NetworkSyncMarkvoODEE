function a=pci(t,x,y,u0,v0)
    v=1;
    l1=3.0608;
    a1=norm(x-y)/(-v*(2*l1^0.5+v));
    a=(norm(u0-v0)*exp(-2*l1^0.5-v)*t-a1*(exp(-(2*l1^0.5-v)*t)-1))^0.5;
end
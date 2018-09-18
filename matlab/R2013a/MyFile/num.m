node_num=10;
dx_num=node_num*3+4;

%init=-30 + 60.*rand(dx_num,1);
init=-30+60.*rand(dx_num,1);

%opts = odeset('RelTol',1e-5,'AbsTol',1e-8);


[t,y]=ode45('LCODE',[0:0.01:5],init);

et=0;
for i=1:node_num
    for j=1:3
        et=et+(y(:,j)-y(:,3*i+j)).^2;
    end
end

figure;
plot(t,(et/node_num).^(1/2));
xlabel('t')
ylabel('E(t)')

figure;
plot(t,y);
xlabel('t')
ylabel('x_{ij}(t)')
node_num=5;
dx_num=node_num*3+3+1;

init=-30 + 60.*rand(dx_num,1);
%init=[12;7;16;18;-17;24;24;26;-14;-1;-15;-8;2;-7;19;-11;10;9;0.001];

%opts = odeset('RelTol',1e-1,'AbsTol',1e-2);


[t,y]=ode45('LCODEs',[0:0.01:4],init);

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

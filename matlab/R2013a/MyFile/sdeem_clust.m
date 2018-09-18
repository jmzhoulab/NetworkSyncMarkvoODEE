%   dX(t) = f(t,X(t))dt + g(t,X(t))dW(t)
%   f(t,X(t))   dirft function
%   g(t,X(t))   diffuse function
%   Euler-Maruyama (EM) discrete method
%
%   X(0) = x(0)
%   X(t+1)=X(t)+f(t,X(t))dt+g(t,X(t))\sqrt(dt)(dW(t+1)-d(W(t)))
%
%  USAGE: sdeem([-1:0],[0:10])
%
function sol = sdeem_clust

    
    inittspan = -1:0;
    tspan = 0:20;    
    
    global NODENUM NODEDIM DIM dt A B D T
    NODENUM=18;
    CLUSTERNUM=3;
    NODEDIM=2;
    
    dt =0.001;
    T = 4;
    
    L_size=[6,6,6];
%    L = -[3 -1 -1 -1 0 0 0 0.1 -0.1 0 0 0 0 0.1 -0.1 0 0 0;-1 4 -1 -1 -1 0 0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 -0.1 -0.1 0.1 -0.1 0.1;-1 -1 5 -1 -1 -1 -0.1 -0.1 0.1 0.1 0 0 -0.1 -0.1 0.2 0.1 -0.1 0;-1 -1 -1 5 -1 -1 0 0.1 0.1 -0.1 -0.1 0 0 0.1 0.1 -0.2 0 0;0 -1 -1 -1 4 -1 0 -0.1 0 -0.1 0.1 0.1 0 -0.1 -0.1 0 0.1 0.1;0 0 -1 -1 -1 3 0 0.1 0 0 0.1 -0.2 0 0.1 0 0 0.1 -0.2;0 0.1 -0.1 0 0 0 5 -1 -1 -1 -1 -1 -0.1 0.1 -0.1 0.1 -0.1 0.1;0.1 -0.1 -0.1 0.1 -0.1 0.1 -1 4 -1 -1 -1 0 0.1 0 -0.1 0 0 0;-0.1 -0.1 0.1 0.1 0 0 -1 -1 3 -1 0 0 -0.1 -0.1 0.1 0.1 -0.1 0.1;0 0.1 0.1 -0.1 -0.1 0 -1 -1 -1 5 -1 -1 0.1 0 0.1 -0.2 0 0;0 -0.1 0 -0.1 0.1 0.1 -1 -1 0 -1 4 -1 -0.1 0 -0.1 0 0.2 0;0 0.1 0 0 0.1 -0.2 -1 0 0 -1 -1 3 0.1 0 0.1 0 0 -0.2;0 0.1 -0.1 0 0 0 -0.1 0.1 -0.1 0.1 -0.1 0.1 5 -1 -1 -1 -1 -1;0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 0 -0.1 0 0 0 -1 4 -1 -1 -1 0;-0.1 -0.1 0.2 0.1 -0.1 0 -0.1 -0.1 0.1 0.1 -0.1 0.1 -1 -1 4 -1 -1 0;0 0.1 0.1 -0.2 0 0 0.1 0 0.1 -0.2 0 0 -1 -1 -1 4 0 -1;0 -0.1 -0.1 0 0.1 0.1 -0.1 0 -0.1 0 0.2 0 -1 -1 -1 0 4 -1;0 0.1 0 0 0.1 -0.2 0.1 0 0.1 0 0 -0.2 -1 0 0 -1 -1 3];
L1=[-150	50	50	50	0	0	0	-0.1	0.1	0	0	0	0	-0.1	0.1	0	0	0;
50	-200	50	50	50	0	-0.1	0.1	0.1	-0.1	0.1	-0.1	-0.1	0.1	0.1	-0.1	0.1	-0.1;
50	50	-250	50	50	50	0.1	0.1	-0.1	-0.1	0	0	0.1	0.1	-0.2	-0.1	0.1	0;
50	50	50	-250	50	50	0	-0.1	-0.1	0.1	0.1	0	0	-0.1	-0.1	0.2	0	0;
0	50	50	50	-200	50	0	0.1	0	0.1	-0.1	-0.1	0	0.1	0.1	0	-0.1	-0.1;
0	0	50	50	50	-150	0	-0.1	0	0	-0.1	0.2	0	-0.1	0	0	-0.1	0.2;
0	-0.1	0.1	0	0	0	-125	25	25	25	25	25	0.1	-0.1	0.1	-0.1	0.1	-0.1;
-0.1	0.1	0.1	-0.1	0.1	-0.1	25	-100	25	25	25	0	-0.1	0	0.1	0	0	0;
0.1	0.1	-0.1	-0.1	0	0	25	25	-75	25	0	0	0.1	0.1	-0.1	-0.1	0.1	-0.1;
0	-0.1	-0.1	0.1	0.1	0	25	25	25	-125	25	25	-0.1	0	-0.1	0.2	0	0;
0	0.1	0	0.1	-0.1	-0.1	25	25	0	25	-100	25	0.1	0	0.1	0	-0.2	0;
0	-0.1	0	0	-0.1	0.2	25	0	0	25	25	-75	-0.1	0	-0.1	0	0	0.2;
0	-0.1	0.1	0	0	0	0.1	-0.1	0.1	-0.1	0.1	-0.1	-145	29	29	29	29	29;
-0.1	0.1	0.1	-0.1	0.1	-0.1	-0.1	0	0.1	0	0	0	29	-116	29	29	29	0;
0.1	0.1	-0.2	-0.1	0.1	0	0.1	0.1	-0.1	-0.1	0.1	-0.1	29	29	-116	29	29	0;
0	-0.1	-0.1	0.2	0	0	-0.1	0	-0.1	0.2	0	0	29	29	29	-116	0	29;
0	0.1	0.1	0	-0.1	-0.1	0.1	0	0.1	0	-0.2	0	29	29	29	0	-116	29;
0	-0.1	0	0	-0.1	0.2	-0.1	0	-0.1	0	0	0.2	29	0	0	29	29	-87

];

L2=[-0.3	0.1	0.1	0.1	0	0	0	-0.1	0.1	0	0	0	0	-0.1	0.1	0	0	0;
0.1	-0.4	0.1	0.1	0.1	0	-0.1	0.1	0.1	-0.1	0.1	-0.1	-0.1	0.1	0.1	-0.1	0.1	-0.1;
0.1	0.1	-0.5	0.1	0.1	0.1	0.1	0.1	-0.1	-0.1	0	0	0.1	0.1	-0.2	-0.1	0.1	0;
0.1	0.1	0.1	-0.5	0.1	0.1	0	-0.1	-0.1	0.1	0.1	0	0	-0.1	-0.1	0.2	0	0;
0	0.1	0.1	0.1	-0.4	0.1	0	0.1	0	0.1	-0.1	-0.1	0	0.1	0.1	0	-0.1	-0.1;
0	0	0.1	0.1	0.1	-0.3	0	-0.1	0	0	-0.1	0.2	0	-0.1	0	0	-0.1	0.2;
0	-0.1	0.1	0	0	0	-0.5	0.1	0.1	0.1	0.1	0.1	0.1	-0.1	0.1	-0.1	0.1	-0.1;
-0.1	0.1	0.1	-0.1	0.1	-0.1	0.1	-0.4	0.1	0.1	0.1	0	-0.1	0	0.1	0	0	0;
0.1	0.1	-0.1	-0.1	0	0	0.1	0.1	-0.3	0.1	0	0	0.1	0.1	-0.1	-0.1	0.1	-0.1;
0	-0.1	-0.1	0.1	0.1	0	0.1	0.1	0.1	-0.5	0.1	0.1	-0.1	0	-0.1	0.2	0	0;
0	0.1	0	0.1	-0.1	-0.1	0.1	0.1	0	0.1	-0.4	0.1	0.1	0	0.1	0	-0.2	0;
0	-0.1	0	0	-0.1	0.2	0.1	0	0	0.1	0.1	-0.3	-0.1	0	-0.1	0	0	0.2;
0	-0.1	0.1	0	0	0	0.1	-0.1	0.1	-0.1	0.1	-0.1	-0.5	0.1	0.1	0.1	0.1	0.1;
-0.1	0.1	0.1	-0.1	0.1	-0.1	-0.1	0	0.1	0	0	0	0.1	-0.4	0.1	0.1	0.1	0;
0.1	0.1	-0.2	-0.1	0.1	0	0.1	0.1	-0.1	-0.1	0.1	-0.1	0.1	0.1	-0.4	0.1	0.1	0;
0	-0.1	-0.1	0.2	0	0	-0.1	0	-0.1	0.2	0	0	0.1	0.1	0.1	-0.4	0	0.1;
0	0.1	0.1	0	-0.1	-0.1	0.1	0	0.1	0	-0.2	0	0.1	0.1	0.1	0	-0.4	0.1;
0	-0.1	0	0	-0.1	0.2	-0.1	0	-0.1	0	0	0.2	0.1	0	0	0.1	0.1	-0.3
];
    A=L1;
    B=L2;
    D = zeros(NODENUM, NODENUM);
    D(1,1)=500;
    D(7,7)=500;
    D(13,13)=500;
        
    
    % the dim of X(t)
    
    DIM =(NODENUM+CLUSTERNUM)*NODEDIM;
    
    t0=tspan(1);                %0
    tfinal=tspan(end);          %10
    tbegin=inittspan(1);        %-1
    
    tspan=t0:dt:tfinal;
    inittspan=tbegin:dt:t0;
    
    ntspan=numel(tspan);
    ninittspan=numel(inittspan);
    
    % the Wiener process increments (ith simulation)
    NORMRAND = [zeros(DIM, 1), randn(DIM, ntspan)];
    dW = sqrt(dt)*NORMRAND;
    
    x=[sdehistory(inittspan), zeros(DIM, ntspan)];
    t=[inittspan,tspan];
    dw=[zeros(DIM, ninittspan), dW];
    
   
    
    %% Euler-Maruyama (EM) discrete method
    for ti=ninittspan:(ninittspan+ntspan-1)        
        x(:,ti+1)=x(:,ti)+sdedrift(ti,t(ti),x(:,ti),x)*dt+sdediffuse(ti,t(ti),x(:,ti),x)*dw(:,ti);
    end
    %% solution  
    sol.t=tspan';
    sol.x=x(:,ninittspan+1:ninittspan+ntspan);
    sol.w=cumsum(dw(:,ninittspan+1:ninittspan+ntspan),2);
    %%
    figure;
    plot(sol.t,sol.x);
    xlabel('t');
    ylabel('x_{ij}(t)');
%   
    figure;
    for i=1:2:41
        plot(sol.t,sol.x(i,:),'b');hold on;
    end
    xlabel('Time');
    ylabel('x_{i1}(t)');  
    
   
    figure;
    for i=2:2:42
        plot(sol.t,sol.x(i,:),'r');hold on;
    end
    xlabel('Time');
    ylabel('x_{i2}(t)');
%     
    ei1=0;
    for i=1:6
        ei1=ei1+(sol.x(2*i+5,:)-sol.x(1,:)).^2+(sol.x(2*i+17,:)-sol.x(3,:)).^2+(sol.x(2*i+29,:)-sol.x(5,:)).^2;
    end
    figure;
    plot(sol.t,ei1.^(1/2));
    xlabel('Time');
    ylabel('e_{i1}(t)');
    
    ei2=0;
    for i=1:6
      ei2=ei2+(sol.x(2*i+6,:)-sol.x(2,:)).^2+(sol.x(2*i+18,:)-sol.x(4,:)).^2+(sol.x(2*i+30,:)-sol.x(6,:)).^2;
    end
    figure;
    plot(sol.t,ei2.^(1/2));
    xlabel('Time');
    ylabel('e_{i2}(t)');
    
    e12=0;
    for i=1:6
        e12=e12+(sol.x(2*i+5,:)-sol.x(2*i+17,:)).^2+(sol.x(2*i+6,:)-sol.x(2*i+18,:)).^2;
    end
    figure;
    plot(sol.t,e12.^(1/2));
    xlabel('Time');
    ylabel('e_{12}(t)');
    
    e13=0;
    for i=1:6
        e13=e13+(sol.x(2*i+5,:)-sol.x(2*i+29,:)).^2+(sol.x(2*i+6,:)-sol.x(2*i+30,:)).^2;
    end
    figure;
    plot(sol.t,e13.^(1/2));
    xlabel('Time');
    ylabel('e_{13}(t)');
    
    e23=0;
    for i=1:6
        e23=e23+(sol.x(2*i+17,:)-sol.x(2*i+29,:)).^2+(sol.x(2*i+18,:)-sol.x(2*i+30,:)).^2;
    end
    figure;
    plot(sol.t,e23.^(1/2));
    xlabel('Time');
    ylabel('e_{23}(t)');
    
end

function fx = sdedrift(ti,t,x,X)
    %%   drift function F(t,X(t))
    global DIM NODENUM NODEDIM A B D T;
    fx=zeros(DIM,1);
    delay=sdedelay(t);
    x_lag_2=X(:,ti-delay(2));              % x(t-\tau_2(t))
    %
   
    
    % s(t)
    s1_arr= 1:NODEDIM;
    s2_arr = NODEDIM +(1:NODEDIM);
    s3_arr = 2*NODEDIM +(1:NODEDIM);
    fx(s1_arr)=f1(ti,t,x(s1_arr),X(s1_arr,:));
    fx(s2_arr)=f2(ti,t,x(s2_arr),X(s2_arr,:));
    fx(s3_arr)=f3(ti,t,x(s3_arr),X(s3_arr,:));
    
    s1=x(s1_arr);
    s2=x(s2_arr);
    s3=x(s3_arr);
    
    
    % coupling
    
    for i=1:NODENUM
        Xi=X(4+NODEDIM*i+1:4+NODEDIM*i+NODEDIM,:);
        xi=x(4+NODEDIM*i+1:4+NODEDIM*i+NODEDIM);
        xitau=x_lag_2(4+NODEDIM*i+1:4+NODEDIM*i+NODEDIM);
        
        if 1<i|i==1&&6>i|i==6
            ux=f1(ti,t,xi,Xi);
        elseif 7<i|i==7&&12>i|i==12
            ux=f2(ti,t,xi,Xi);
        elseif 13<i|i==13&&18>i|i==18
            ux=f3(ti,t,xi,Xi);
        end
        
        sita=0.75;
        if mod(t,  T) < sita*T|mod(t,  T) == sita*T && i == 1
            ux  = ux  - D(i,i)*(xi-s1);    
        elseif mod(t,  T) < sita*T |mod(t,  T) == sita*T && i == 7
            ux  = ux  - D(i,i)*(xi-s2);
        elseif mod(t,  T) < sita*T |mod(t,  T) == sita*T && i == 13
            ux  = ux  - D(i,i)*(xi-s3); 
        end
        
        
        for j=1:NODENUM
            if j ~= i
                xj=x(4+NODEDIM*j+1:4+NODEDIM*j+NODEDIM);
                xjtau=x_lag_2(4+NODEDIM*j+1:4+NODEDIM*j+NODEDIM);
                ux=ux+A(i,j)*(g(xj)-g(xi))+B(i,j)*(h(xjtau)-h(xitau));
            end
        end
        fx(4+NODEDIM*i+1: 4+NODEDIM*i+NODEDIM)=ux;
    end
end

function  fx =f1(i,t,x,X)
    delay=sdedelay(t);
   
    x_lag_1=X(:,i-delay(1));              % x(t-\tau_1(t))
    C=[-1 0;0 -1];
    A1=[2 -0.1;-5 4.5];
    B1=[-1.5 -0.1;-0.2 -4];

    fx=C*x+A1*tanh(x)+B1*tanh(x_lag_1);
end

function fx=f2(i,t,x,X)
    delay=sdedelay(t);
    x_lag_1=X(:,i-delay(1));              % x(t-\tau_1(t))
    C=[-1 0;0 -1];
    A2=[0.1 -0.1;-0.4 0.4];
    B2=[-0.5 -0.1;-0.1 -0.4];
    fx=C*x+A2*tanh(x)+B2*tanh(x_lag_1);
end

function fx=f3(i,t,x,X)
    delay=sdedelay(t);
    x_lag_1=X(:,i-delay(1));              % x(t-\tau_1(t))
    C=[-1 0;0 -1];
    A3=3*[0.2 -0.1;-0.3 0.4];
    B3=[-0.5 -0.2;-0.1 -0.3];
    fx=C*x+A3*tanh(x)+B3*tanh(x_lag_1);
end
    
    
function gx1=g(x)
    T=[1 0;0 0];
    gx1=T*x+0.005*T*sin(x);
end

function hx=h(x)
    T=[1 0;0 0];
    hx=T*0.1*sin(x);
end


function gx = sdediffuse(i,t,x,X)
    %   diffuse function G(t,X(t))
    global DIM;
    gx=zeros(DIM,DIM);
    
    
    delay=sdedelay(t);
    %     x_lag_1=X(i-delay(1),:);              % x(t-\tau_1(t))
    %     x_lag_2=X(i-delay(2),:);              % x(t-\tau_2(t))
    
    
    for n=7:18
        if mod(n,2)==1
            gx(n,n)=0.1*(x(n)-x(1));
        elseif mod(n,2)==0
            gx(n,n)=0.1*(x(n)-x(2));
        end
    end
    
    for n=19:30
        if mod(n,2)==1
            gx(n,n)=0.1*(x(n)-x(3));
        elseif mod(n,2)==0
            gx(n,n)=0.1*(x(n)-x(4));
        end
    end
    
    for n=31:42
        if mod(n,2)==1
            gx(n,n)=0.1*(x(n)-x(5));
        elseif mod(n,2)==0
            gx(n,n)=0.1*(x(n)-x(6));
        end
    end
end

function xinit = sdehistory(t)
    global DIM
    nt=numel(t);
    xinit=zeros(DIM,nt);
%     i=1;
%     xinit(i,:)=log(t+2);
%     i=2;
%     xinit(i,:)=sin(t);
%     i=3;
%     xinit(i,:)=cosh(t);
%     i=4;
%     xinit(i,:)=tan(t);
%     i=5;
%     xinit(i,:)=tanh(t);
%     i=6;
%     xinit(i,:)=sinh(t);
    for i=1:6
        xinit(i,:)=i*log(2*t+4);
    end
    for i=7:42
        xinit(i,:)=0.1*i*cos(t);
    end
%     for i=19:30
%         xinit(i,:)=exp(t);
%     end
%     for i=31:42
%         xinit(i,:)=i*t.^2;
%     end
end

function delay=sdedelay(t)
    global dt
    ndelay=2;
    delay=zeros(ndelay);
    delay(1)=0.01*exp(2*t)/(2+exp(2*t));
    delay(2)=0.01*exp(t)/(1+exp(t));
    delay=int16(delay/dt);      % the multiple of dt
end
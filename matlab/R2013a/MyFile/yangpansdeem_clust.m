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
    
    global NODENUM NODEDIM DIM dt A B hi T L u p ki
    NODENUM=18;
    CLUSTERNUM=3;
    NODEDIM=2;
    
    dt =0.1;
    T = 4;
   
    L_size=[6,6,6];
L1 = -[3 -1 -1  -1 0   0    0 0.1 -0.1 0 0 0 0 0.1 -0.1 0 0 0;
      -1  4 -1  -1 -1  0    0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 -0.1 -0.1 0.1 -0.1 0.1;
      -1 -1  5  -1 -1 -1     -0.1 -0.1 0.1 0.1 0 0 -0.1 -0.1 0.2 0.1 -0.1 0;
      -1 -1 -1   5 -1 -1      0 0.1 0.1 -0.1 -0.1 0 0 0.1 0.1 -0.2 0 0;
       0 -1 -1  -1  4 -1      0 -0.1 0 -0.1 0.1 0.1 0 -0.1 -0.1 0 0.1 0.1;
       0  0 -1  -1 -1  3      0 0.1 0 0 0.1 -0.2 0 0.1 0 0 0.1 -0.2;
    0 0.1 -0.1 0 0 0            5  -1 -1 -1 -1 -1    -0.1 0.1 -0.1 0.1 -0.1 0.1;
   0.1 -0.1 -0.1 0.1 -0.1 0.1   -1  4 -1 -1 -1  0     0.1 0 -0.1 0 0 0;
    -0.1 -0.1 0.1 0.1 0 0       -1 -1  3 -1  0  0      -0.1 -0.1 0.1 0.1 -0.1 0.1;
    0 0.1 0.1 -0.1 -0.1 0       -1 -1 -1  5 -1 -1       0.1 0 0.1 -0.2 0 0;
    0 -0.1 0 -0.1 0.1 0.1       -1 -1  0 -1  4 -1     -0.1 0 -0.1 0 0.2 0;
    0 0.1 0 0 0.1 -0.2          -1 0   0 -1 -1  3     0.1 0 0.1 0 0 -0.2;
    0 0.1 -0.1 0 0 0 -0.1 0.1 -0.1 0.1 -0.1 0.1          5 -1 -1 -1 -1 -1;
    0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 0 -0.1 0 0 0         -1  4 -1 -1 -1  0;
    -0.1 -0.1 0.2 0.1 -0.1 0 -0.1 -0.1 0.1 0.1 -0.1 0.1 -1 -1  4 -1 -1  0;
    0 0.1 0.1 -0.2 0 0 0.1 0 0.1 -0.2 0 0               -1 -1 -1  4  0 -1;
    0 -0.1 -0.1 0 0.1 0.1 -0.1 0 -0.1 0 0.2 0           -1 -1 -1  0  4 -1;
    0 0.1 0 0 0.1 -0.2 0.1 0 0.1 0 0 -0.2               -1  0  0 -1 -1  3];
L2 = -[3 -1 -1 -1 0 0 0 0.1 -0.1 0 0 0 0 0.1 -0.1 0 0 0;
    -1 4 -1 -1 -1 0 0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 -0.1 -0.1 0.1 -0.1 0.1;
    -1 -1 5 -1 -1 -1 -0.1 -0.1 0.1 0.1 0 0 -0.1 -0.1 0.2 0.1 -0.1 0;
    -1 -1 -1 5 -1 -1 0 0.1 0.1 -0.1 -0.1 0 0 0.1 0.1 -0.2 0 0;
    0 -1 -1 -1 4 -1 0 -0.1 0 -0.1 0.1 0.1 0 -0.1 -0.1 0 0.1 0.1;
    0 0 -1 -1 -1 3 0 0.1 0 0 0.1 -0.2 0 0.1 0 0 0.1 -0.2;
    0 0.1 -0.1 0 0 0 5 -1 -1 -1 -1 -1 -0.1 0.1 -0.1 0.1 -0.1 0.1;
    0.1 -0.1 -0.1 0.1 -0.1 0.1 -1 4 -1 -1 -1 0 0.1 0 -0.1 0 0 0;
    -0.1 -0.1 0.1 0.1 0 0 -1 -1 3 -1 0 0 -0.1 -0.1 0.1 0.1 -0.1 0.1;
    0 0.1 0.1 -0.1 -0.1 0 -1 -1 -1 5 -1 -1 0.1 0 0.1 -0.2 0 0;
    0 -0.1 0 -0.1 0.1 0.1 -1 -1 0 -1 4 -1 -0.1 0 -0.1 0 0.2 0;
    0 0.1 0 0 0.1 -0.2 -1 0 0 -1 -1 3 0.1 0 0.1 0 0 -0.2;
    0 0.1 -0.1 0 0 0 -0.1 0.1 -0.1 0.1 -0.1 0.1 5 -1 -1 -1 -1 -1;
    0.1 -0.1 -0.1 0.1 -0.1 0.1 0.1 0 -0.1 0 0 0 -1 4 -1 -1 -1 0;
    -0.1 -0.1 0.2 0.1 -0.1 0 -0.1 -0.1 0.1 0.1 -0.1 0.1 -1 -1 4 -1 -1 0;
    0 0.1 0.1 -0.2 0 0 0.1 0 0.1 -0.2 0 0 -1 -1 -1 4 0 -1;
    0 -0.1 -0.1 0 0.1 0.1 -0.1 0 -0.1 0 0.2 0 -1 -1 -1 0 4 -1;
    0 0.1 0 0 0.1 -0.2 0.1 0 0.1 0 0 -0.2 -1 0 0 -1 -1 3];

L3=[-0.3	0.1	0.1	0.1	0	0	0	-0.1	0.1	0	0	0	0	-0.1	0.1	0	0	0;
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
    B=zeros(NODENUM, NODENUM);
    hi =binornd(1,0.9,1,NODENUM);
    G=eye(NODEDIM);
    L1=[10*L1(1:6,1:6),L1(1:6,7:12),L1(1:6,13:18);L1(7:12,1:6),10*L1(7:12,7:12),L1(7:12,13:18);L1(13:18,1:6),L1(13:18,7:12),10*L1(13:18,13:18)];
    L=cat(3,L1,L1,L1);
    p=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];   
    u=1;
    ki=200;
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
    ci=zeros(NODENUM, ntspan);
    cij=zeros(NODENUM,NODENUM,ntspan);
    
    %% Euler-Maruyama (EM) discrete method
      for ti=ninittspan:(ninittspan+ntspan-1)  
        x(:,ti+1)=x(:,ti)+sdedrift(ti,t(ti),ci(:,ti),cij(:,:,ti),x(:,ti),x)*dt+sdediffuse(ti,t(ti),x(:,ti),x)*dw(:,ti);
        ci(:,ti+1)=ci(:,ti)+cidf(x(:,ti));
        cij(:,:,ti+1)=cij(:,:,ti)+cijdf(x(:,ti));
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
%   plot
    figure;
    for i=1:2:41
        plot(sol.t,sol.x(i,:),'b');hold on;
    end
    xlabel('Time');
    ylabel('x_{i1}(t)');  
    for i=1:NODENUM
        plot(sol.t,ci(i,ninittspan:ninittspan+ntspan-1),'b');hold on;
    end
    xlabel('Time');
    ylabel('c_{i}(t)');  
    for i=1:NODENUM
        for j=1:NODENUM
            plot(sol.t,squeeze(cij(i,j,ninittspan:ninittspan+ntspan-1)),'b');hold on;
        end
    end
    xlabel('Time');
    ylabel('c_{ij}(t)'); 
    figure;
    for i=2:2:42
        plot(sol.t,sol.dW(i,:),'r');hold on;
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

function fx = sdedrift(ti,t,cci,ccij,x,X)
    %%   drift function F(t,X(t))
    global DIM NODENUM NODEDIM A B hi L u p;
    v=Markov(u,p);
    A=L(:,:,v);
    A=ccij.*A;
    pho=binornd(1,0.8,1,NODENUM);
    fx=zeros(DIM,1);
    delay=sdedelay(t);
    x_lag_2=X(:,ti-delay(2));              % x(t-\tau_2(t))
    
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
            ux=f1(ti,t,xi,Xi)- hi(i)*pho(i)*cci(i)*(xi-s1);
        elseif 7<i|i==7&&12>i|i==12
            ux=f2(ti,t,xi,Xi)-hi(i)*pho(i)*cci(i)*(xi-s2);
        elseif 13<i|i==13&&18>i|i==18
            ux=f3(ti,t,xi,Xi)- hi(i)*pho(i)*cci(i)*(xi-s3); 
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
function ci=cidf(x)
        global NODENUM NODEDIM dt ki hi;
        ci=zeros(NODENUM,1);
        s1_arr= 1:NODEDIM;
        s2_arr = NODEDIM +(1:NODEDIM);
        s3_arr = 2*NODEDIM +(1:NODEDIM);
        s1=x(s1_arr);
        s2=x(s2_arr);
        s3=x(s3_arr);
        for i=1:NODENUM
            xi=x(4+NODEDIM*i+1:4+NODEDIM*i+NODEDIM);
            if 1<i|i==1&&6>i|i==6
                ci(i)=hi(i)*ki*(xi-s1)'*(xi-s1)*dt;
            elseif 7<i|i==7&&12>i|i==12
                ci(i)=hi(i)*ki*(xi-s2)'*(xi-s2)*dt;
            elseif 13<i|i==13&&18>i|i==18
                ci(i)=hi(i)*ki*(xi-s3)'*(xi-s3)*dt;
            end
        end
end
function cij=cijdf(x)
        global NODENUM NODEDIM A dt;
        cij=zeros(NODENUM,NODENUM);
        hij=1;
        for i=1:NODENUM
             xi=x(4+NODEDIM*i+1:4+NODEDIM*i+NODEDIM);
            for j=1:NODENUM
                if abs(j-i)>=6  %判断是否同簇
                    hij=0;
                end
                xj=x(4+NODEDIM*j+1:4+NODEDIM*j+NODEDIM);
                cij(i,j)=hij*A(i,j)*(xi-xj)'*(xi-xj)*dt;
            end
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
function v=Markov(u,p)
    %v=[0];
    n=length(p);
    s=1:n;
    pp=[p(u,:);s];
    random=rand(1);  %产生[0，1]上均匀分布的随机数，用于刻画转移的随机性
    min=0;max=pp(1,1);
    for i=1:n-1
        if min<random&&random<=max
            v=pp(2,1);
            break;
        end
        min=min+pp(1,i);max=max+pp(1,i+1);
        if min<random&&random<=max
            v=pp(2,i+1);
            break;
        end
    end
    %v=[v,random,min,max];
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
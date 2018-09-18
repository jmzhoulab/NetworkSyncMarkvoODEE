function sol =CSNL
    
    opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
    node_num=10;
    cluster_num=2;
    L_size=[5,5];
    L1=[-2 1 1 0 0 0 0 0 0 0;
        1 -2 0 1 0 0 0 -0.1 0 0.1;
        1 0 -2 0 1 0 0 0 0 0;
        0 1 0 -2 1 0 0 0 0 0;
        0 0 1 1 -2 0 0 0 0 0;
        0 0 0 0 0 -2 1 0 0 1;
        0 0 0 0 0 1 -2 1 0 0;
        0 0 0 0 0 0 1 -2 1 0;
        -0.1 0 0 0 0.1 0 0 1 -2 1;
        0 0 0 0 0 1 0 0 1 -2];
    L2=[-2 0 1 0 0 0 1 0 0 0;
        0 -2 0 1 0 0 0 -0.1 0 0.1;
        1 0 -2 0 1 0 0 0 0 0;
        0 1 0 -2 1 0 0 0 0 0;
        0 0 1 1 -2 0 0 0 0 0;
        1 0 0 0 0 -3 1 0 0 1;
        0 0 0 0 0 1 -2 1 0 0;
        0 0 0 0 0 0 1 -2 1 0;
        -0.1 0 0 0 0.1 0 0 1 -2 1;
        0 0 0 0 0 1 0 0 1 -2];
    L3=[-3 1 1 0 0 0 0 0 0 1;
        1 -2 0 1 0 0 0 -0.1 0 0.1;
        1 0 -2 0 1 0 0 0 0 0;
        0 1 0 -2 1 0 0 0 0 0;
        0 0 1 1 -2 0 0 0 0 0;
        0 0 0 0 0 -2 1 0 0 1;
        0 0 0 0 0 1 -2 0 0 1;
        0 0 0 0 0 0 0 -1 1 0;
        -0.1 0 0 0 0.1 0 0 1 -2 1;
        1 0 0 0 0 1 1 0 1 -4];
    L=L2;
    LL=cat(3,L1,L2,L3);
    a=0.05;
    P=diag([5,1,1]);
    
    cij_num = 24;
    cpin_num = 2;
    hi=1;
    M = L;
    M(1:11:100) = 0;
    M = abs(M);
    CIJ = find(M ~= 0);
    E=zeros(node_num,node_num);
      E(1:11:end) = 1;
    pr=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];   
    u=1;
    
    dx_num=(cluster_num+node_num)*3 + cij_num + cpin_num; % (10+2)*3 + 24 + 2
    init = rand(dx_num,1);
    init(1:36) = (init(1:36)-1)*30;
    init(37:62) = init(37:62);
    [t,x]=ode45(@LCODEs,[0:0.001:1.5],init,opts,LL,u,pr,hi,CIJ,E,P,L_size,cluster_num, cij_num , dx_num,a);
    sol.t = t;
    sol.x = x;
    
    figure;    
    plot(t,x(:,1:36));    
    xlabel('t');
    ylabel('x_i^j(t)');
    figure;    
    plot(t,x(:,1:3:36));    
    xlabel('t');
    ylabel('x_i^1(t)');
    figure;    
    plot(t,x(:,2:3:36));    
    xlabel('t');
    ylabel('x_i^2(t)');
    figure;    
    plot(t,x(:,3:3:36));    
    xlabel('t');
    ylabel('x_i^3(t)');
    
    figure;
    
    plot(t,x(:,37:60));
    
    xlabel('t');
    ylabel('c_{ij}(t)');
    
    et=0;
    beg_i=cluster_num;
    for clu=1:cluster_num
        end_i=beg_i+L_size(clu)-1;
        for i=beg_i:end_i
            for j=1:3
                et=et+(x(:,3*i+j)-x(:,3*(clu-1)+j)).^2;
            end
        end
        beg_i=end_i+1;
    end
    figure;
    plot(t,(et/node_num).^(1/2));
    xlabel('t');
    ylabel('E(t)');
    
    figure;
    plot(t,x(:,dx_num-1:dx_num));
    xlabel('t');
    ylabel('c_i(t)')
%     
end


%% LCODEs
function dx=LCODEs(t,x,LL,u,pr,hi,CIJ,E,P,L_size,cluster_num,cij_num , dx_length,a)
    L=LL(:,:,u);
    L_dim=size(L,1);
    node_num=L_dim;
%     disp(t);
    
    dx=zeros(dx_length,1);
    
    %s_j(t), The number of s_j(t) is cluster_num (2)
    for clu=1:cluster_num
        xi(1)=x(3*(clu-1)+1);
        xi(2)=x(3*(clu-1)+2);
        xi(3)=x(3*(clu-1)+3);
        fxi=f(xi,clu);
        dx(3*(clu-1)+1)=fxi(1);
        dx(3*(clu-1)+2)=fxi(2);
        dx(3*(clu-1)+3)=fxi(3);
    end
    
    %c_ij(t)
    for dxi = 1 : cij_num
        j = mod(CIJ(dxi),node_num);
        i = ceil(CIJ(dxi)/10);
        if j == 0
            j = 10;
        end
        
        dx((cluster_num+node_num)*3 + dxi) = a*L(i,j)*(x(3*(i+cluster_num - 1)+(1:3))-x(3*(j+cluster_num - 1)+(1:3)))'*P*(x(3*(i+cluster_num - 1)+(1:3))-x(3*(j+cluster_num - 1)+(1:3)));
        
    end
    
    %c_i(t) 
    
    
    dx(dx_length-1) = a*(x(3*(cluster_num)+(1:3))-x((1:3)))'*P*(x(3*(cluster_num)+(1:3))-x((1:3)));
    dx(dx_length) = a*(x(3*(5+cluster_num)+(1:3))-x(3+(1:3)))'*P*(x(3*(5+cluster_num)+(1:3))-x(3+(1:3)));
    
    
    %dx(dx_length)=0;
    
    %CODE
    beg_i=cluster_num;
    for clu=1:cluster_num
        end_i=beg_i+L_size(clu)-1;
        for i=beg_i:end_i
            
            xi=zeros(3,1);
            xi(1)=x(3*i+1);
            xi(2)=x(3*i+2);
            xi(3)=x(3*i+3);
            
            fxi=f(xi,clu);
            
            ux=zeros(3,1);
            
            ux(1)=fxi(1)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(1))-h(x(3*(clu-1)+1)))+0.2*randn();
            ux(2)=fxi(2)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(2))-h(x(3*(clu-1)+2)))+0.2*randn();
            ux(3)=fxi(3)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(3))-h(x(3*(clu-1)+3)))+0.2*randn();
            
            for j=1:node_num
                
                xj=zeros(3,1);
                xj(1)=x(3*(j+cluster_num-1)+1);
                xj(2)=x(3*(j+cluster_num-1)+2);
                xj(3)=x(3*(j+cluster_num-1)+3);
                dcij_ind = find(CIJ==(i-cluster_num+1)+j*10)+36;
                if isempty( dcij_ind) ~= 1
                    
                    ux=ux+ binornd(1,hi)*L(i-cluster_num+1,j)*x(dcij_ind)*(h(xj)-h(xi))+0.2*randn(3,1);
                end
            end
            dx(3*i+1)=ux(1);
            dx(3*i+2)=ux(2);
            dx(3*i+3)=ux(3);
        end
        beg_i=end_i+1;
        u=Markov(u,pr);
        L=LL(:,:,u);
    end
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
%f(x)
%{
function  fx=f(x,X,delay,i)
    x_lag_1=X(:,delay); 
    if i == 1
        C=[-1 1,0;0 1,-1; 0,0,1];
        A=[2 0.2,-0.1;-5 0.8,4.5;-2,0.3, 1];
        B=[-1.5,-1, -0.1;-0.2 0.3,-4;2,-3.5,1];
    elseif i == 2
        C=[-1 1,0;0 1,-1; 0,0,1];
        A=[0.1 0.2 -0.1;-0.4 0.2 0.4;1, -2.3, 0];
        B=[-0.5 -0.3 -0.1;-0.1 -0.2 -0.4;-0.2 -0.1 0.3];
    elseif i == 3
        C=[-1 1,0;0 1,-1; 0,0,1];
        A=3*[0.2 1.2 -0.1;-0.3 -0.2 0.4; 1 -1 0.3];
        B=[-0.5 -0.1 -0.2;-0.1 -0.3 -0.3;-0.2 0.1 -0.4];
    end
     fx=C*x+A*tanh(x)+B*tanh(x_lag_1);
end
%}
function  fx=f(x,i)
    if i == 1
        fx(1)=-x(1)+x(2)+0.1*tanh(x(1))+0.2*tanh(x(2))+0.4*tanh(x(3));
        fx(2)=x(2)-x(3)-0.4*tanh(x(1))+0.2*tanh(x(2))+0.4*tanh(x(3));
        fx(3)=x(3)+1*tanh(x(1))-2.3*tanh(x(2));
        %C=[-1 1,0;0 1,-1; 0,0,1];
        %A=[2 0.2,-0.1;-5 0.8,4.5;-2,0.3, 1];
        %B=[-1.5,-1, -0.1;-0.2 0.3,-4;2,-3.5,1];
    elseif i == 2
        fx(1)=-x(1)+x(2)+0.1*tanh(x(1))+0.2*tanh(x(2))-0.1*tanh(x(3));
        fx(2)=x(2)-x(3)-0.4*tanh(x(1))+0.2*tanh(x(2))+0.4*tanh(x(3));
        fx(3)=x(3)+1*tanh(x(1))-2.3*tanh(x(2));
        %C=[-1 1,0;0 1,-1; 0,0,1];
        %A=[0.1 0.2 -0.1;-0.4 0.2 0.4;1, -2.3, 0];
        %B=[-0.5 -0.3 -0.1;-0.1 -0.2 -0.4;-0.2 -0.1 0.3];
    elseif i == 3
        fx(1)=-x(1)+x(2)+0.6*tanh(x(1))+3.6*tanh(x(2))-0.3*tanh(x(3));
        fx(2)=x(2)-x(3)-0.9*tanh(x(1))-0.6*tanh(x(2))+1.2*tanh(x(3));
        fx(3)=x(3)+3*tanh(x(1))-3*tanh(x(2))+0.9*tanh(x(3));
        %C=[-1 1,0;0 1,-1; 0,0,1];
        %A=3*[0.2 1.2 -0.1;-0.3 -0.2 0.4; 1 -1 0.3];
        %B=[-0.5 -0.1 -0.2;-0.1 -0.3 -0.3;-0.2 0.1 -0.4];
    end
     %fx=(C*x'+A*tanh(x'));
end
%{
function  fx=f(x,i)
    if i == 1
        a=0;
    elseif i == 2
        a=0.001;
    elseif i == 3
        a=0.002;
    end
    fx(1)=(25*a+10)*(x(2)- x(1));
    fx(2)=(28-35*a)*x(1) - x(1)*x(3)+(29*a-1)*x(2);
    fx(3)=x(1)*x(2) - (8+a)/3*x(3);
end
%}
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
%h(x)
function hx=h(x)
    hx=5*x+sin(x);
end
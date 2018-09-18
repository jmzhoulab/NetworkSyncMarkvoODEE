function sol =CSNL
    
    opts = odeset('RelTol',1e-4,'AbsTol',1e-9);
    node_num=10;
    cluster_num=2;
    L_size=[5,5];
    
    L=[-2 1 1 0 0 0 0 0 0 0;1 -2 0 1 0 0 0 -0.1 0 0.1;1 0 -2 0 1 0 0 0 0 0;0 1 0 -2 1 0 0 0 0 0;0 0 1 1 -2 0 0 0 0 0;0 0 0 0 0 -2 1 0 0 1;0 0 0 0 0 1 -2 1 0 0;0 0 0 0 0 0 1 -2 1 0;-0.1 0 0 0 0.1 0 0 1 -2 1;0 0 0 0 0 1 0 0 1 -2];
    
    a=0.02;

    P=diag([5,1,1]);
    
    cij_num = 24;
    cpin_num = 2;
    
    M = L;
    M(1:11:100) = 0;
    M = abs(M);
    CIJ = find(M ~= 0);   
    E=zeros(node_num,node_num);
    E(1,1) =1;
    E(6,6) = 1;
    
    
    dx_num=(cluster_num+node_num)*3 + cij_num + cpin_num; % (10+2)*3 + 24 + 2
    init = rand(dx_num,1);
    [t,x]=ode45(@LCODEs,[0:10],init,opts,L,CIJ,E,P,L_size,cluster_num, cij_num , dx_num,a);
    sol.t = t;
    sol.x = x;
    
    figure;
    
   plot(t,x(:,1:36));
    
    xlabel('t');
    ylabel('x_i^j(t)');
    
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
    
end
function dx=LCODEs(t,x,L,CIJ,E,P,L_size,cluster_num,cij_num , dx_length,a)
    L_dim=size(L,1);
    node_num=L_dim;
    %disp(t);
    
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
    
    % c_i(t)
   
    
      
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
            
            ux(1)=fxi(1)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(1))-h(x(3*(clu-1)+1)));
            ux(2)=fxi(2)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(2))-h(x(3*(clu-1)+2)));
            ux(3)=fxi(3)-x(dx_length + clu-2)*E(i-cluster_num+1,i-cluster_num+1)*(h(xi(3))-h(x(3*(clu-1)+3)));
            
            for j=1:node_num
                
                xj=zeros(3,1);
                xj(1)=h(x(3*(j+cluster_num-1)+1));
                xj(2)=h(x(3*(j+cluster_num-1)+2));
                xj(3)=h(x(3*(j+cluster_num-1)+3));
                dcij_ind = find(CIJ==(i-cluster_num+1)+j*10)+36;
                if isempty( dcij_ind) ~= 1
                    ux=ux+L(i-cluster_num+1,j)*x(dcij_ind)*(xj-xi);
                end
            end
            dx(3*i+1)=ux(1);
            dx(3*i+2)=ux(2);
            dx(3*i+3)=ux(3);
        end
        beg_i=end_i+1;
    end
end
    
    %f(x)
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
    %h(x)
function hx=h(x)
       hx=5*x+sin(x);    
end
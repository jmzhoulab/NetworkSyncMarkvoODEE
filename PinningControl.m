function [ output_args ] = PinningControl()
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    global nodenum nodedim
    tic;
    nodenum=100;      %节点数
    nodedim=3;      %节点维数
    kapa=0.4;
    alpha=14.8;
    delta=0.8;
    dt=0.01;
    tn=3;
    t=0:dt:tn;       %迭代时间序列
    inititspan=2*(2*rand(nodenum,nodedim)-0.5);  
    Gama=eye(nodedim);      %内耦合结构矩阵
    s=[1,-1,2];         %目标函数的初值
    timedim=length(t);      %迭代次数
    x=zeros(nodenum,nodedim,timedim);       %生成一个三维的矩阵，用于存放数据，第一维度表示节点，第二维度表示节点状态，第三维度表示时间
    for i=1:timedim     %把时间t放在数组x的第三维中
        x(:,:,i)=t(i);
    end
    x(:,:,1)=inititspan;      %赋予初值
    Sigmma = zeros(nodenum,nodedim,timedim);
    for i=1:nodenum
        for j=1:nodenum
            Sigmma(i,j,1)=norm(x(i,:,1)-x(j,:,1));
        end
    end
    e=zeros(nodenum,nodedim,timedim);   %存放误差
    q=zeros(nodenum,timedim);
    e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %初始误差
    for i=1:nodenum
        q(i,1)=kapa*norm(e(i,:,1));
    end
    for k=1:timedim-1       %时间方向上的循环
        for i=1:nodenum
            cp=zeros(1,nodedim);
            for j=1:nodenum        %计算系统中的求和项，即耦合项的和 
                cp=cp+Sigmma(i,j,k)*g(x(j,:,k))*Gama; 
                Sigmma(i,j,k)=Sigmma(i,j,k)+alpha*norm(x(i,:,1)-x(j,:,1));
            end
            x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))+cp-normrnd(delta,1)*q(i,k)*(x(i,:,k)-s(k,:)))*dt;   %迭代方程
            q(i,k)=q(i,k)+kapa*norm(e(i,:,k));
        end
        s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %目标函数的迭代
        e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算

        if mod(k, 100) == 0
            runtime=toc;
            display(strcat('计算时间：',num2str(runtime),'秒'));
        end
    end
    figplot(t, s', three2two(x), three2two(e))
    runtime=toc;
    display(strcat('总计算时间：',num2str(runtime),'秒'));
end

function gx=g(x)
     gx=zeros(1,3);
     gx(1)=2*x(1)+0.2*sin(x(1));
     gx(2)=2*x(2)+1;
     gx(3)=3*x(3)+0.5*cos(2*x(3));
end

function fx=f(x)
    p=9.78;
    q=14.97;
    m0=-1.31;
    m1=-0.75;
    fx = x;
    fx(1)= p*(-x(1)+x(2)-m1*x(1)+1/2*(m0-m1)*(abs(x(1)+1)-abs(x(1)-1)));   
    fx(2)= x(1)-x(2)+x(3);
    fx(3)= -q*x(2);
end

function m=three2two(x)
    dim=size(x);
    m=zeros(dim(1)*dim(2), dim(3));
    for k=1:dim(3)
        tmp=x(:,:,k)';
        m(:,k)=tmp(:);
    end
end

function figplot(t, s, X, E)
    global nodenum nodedim
    color='rgb';
    figure;     %节点各状态轨迹图
    for j=1:nodedim
        plot(t,s(j,:),color(j));
        hold on;
    end
    for j=1:nodedim
       plot(t,X(j:nodedim:nodenum*nodedim,:),color(j));
       hold on;
    end
    se=zeros(nodedim,length(t));
    figure;     %节点误差轨迹
    for j=1:nodedim
       plot(t,sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)),color(j));
       hold on;
    end
    for j=1:nodedim
        se(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2));
    end
    talle=sum(se);
    figure;
    plot(t,talle);
end


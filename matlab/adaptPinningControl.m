function [kapa,alpha,delta,pinnum,t, s, X, E, q, Sigmma]=adaptPinningControl(kapa,alpha,delta,pinnum,L)
%自适应牵制控制数值模拟
%输入拉普拉斯矩阵
    global nodenum nodedim timedim t s X E q Sigmma;
    tic;
    dt=0.01;
    tn=10;
    t=0:dt:tn;       %迭代时间序列
    nodenum=length(L);    nodedim=3;    timedim=length(t);
    dl = diag(L);
    pd = sort(dl, 'descend');
    D=zeros(nodenum,nodenum);    pin=find(dl>(pd(pinnum)-1));    diagD=diag(D);  diagD(pin)=ones(length(pin),1); D=diag(diagD)+D-diag(diag(D));
    Gama=eye(nodedim);      %内耦合结构矩阵
    q=zeros(nodenum, timedim);
    A=zeros(nodenum, nodenum);
    Sigmma = zeros(sum(sum(triu(L,1))), timedim);
    X=zeros(nodenum*nodedim,timedim);    
    X(:,1)=4*(rand(nodenum*nodedim,1)-0.5);     %定义状态矩阵并进行初始化
    %X(:,1)=kron(ones(nodenum,1),[-6,1,3]');
    E=X;  s=zeros(nodedim,timedim);  
    s(:,1)=4*(rand(nodedim,1)-0.5);
    %s(:,1)=[1,-1,2]';
    E(:,1)=X(:,1)-kron(ones(nodenum,1),s(:,1));

    for k=1:timedim-1
        St=kron(ones(nodenum,1),s(:,k));
        X(:,k+1)=X(:,k)+(F(X(:,k))+kron(A,Gama)*H(X(:,k))-kron(diag(binornd(1,delta,nodenum,1).*q(:,k)).*D,Gama)*(X(:,k)-St))*dt;
        s(:,k+1)=s(:,k)+f(s(:,k))*dt;   %目标函数的迭代
        E(:,k+1)=X(:,k+1)-kron(ones(nodenum,1),s(:,k+1));       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
        font=1;
        for i=1:nodenum
            A(i,i)=0;
            xi=X((i-1)*nodedim+1:i*nodedim,k);
            for j=i+1:nodenum
                if L(i,j)~=0
                    xj=X((j-1)*nodedim+1:j*nodedim,k);
                    A(i,j)=A(i,j)+alpha*norm(xi-xj);
                    A(j,i)=A(i,j);
                    Sigmma(font,k+1)=A(i,j);
                    font=font+1;
                end
            end
            q(i,k+1)=q(i,k)+kapa*norm(xi-s(:,k));
        end
        A=A-diag(sum(A));
        if mod(k, 100) == 0
            display(strcat('计算时间：',num2str(toc),'秒'));
        end
    end
    figplot();
%     name=sprintf('adaptPinningControl_k%sa%sd%sp%s.mat',num2str(kapa),num2str(alpha),num2str(delta),num2str(pinnum));
%     save(name, 'L','X','E','t','s','q','Sigmma')
    display(strcat('总计算时间：',num2str(toc),'秒'));
end

function F=F(X)
global nodenum nodedim;
F=zeros(length(X),1);
for i=1:nodenum
    F((i-1)*nodedim+1:i*nodedim)=f(X((i-1)*nodedim+1:i*nodedim));
end
end
function H=H(X)
    global nodenum nodedim;
    sizexk=size(X);
    H=zeros(sizexk(1),sizexk(2));
    for j=1:sizexk(2)
        for i=1:nodenum
            H((i-1)*nodedim+1:i*nodedim,j)=g(X((i-1)*nodedim+1:i*nodedim,j));
        end
    end
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
function gx=g(x)
 gx=zeros(1,3);
 gx(1)=2*x(1)+0.2*sin(x(1));
 gx(2)=2*x(2)+1;
 gx(3)=3*x(3)+0.5*cos(2*x(3));
end

function figplot()
    global nodenum nodedim t s X E q Sigmma;
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
%     figure;     %节点误差轨迹
%     for j=1:nodedim
%        plot(t,sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)),color(j));
%        hold on;
%     end
    for j=1:nodedim
        se(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2));
    end
    talle=sum(se);
    figure;
    for i=1:nodenum
        plot(t,q(i,:),'r');
        hold on;
    end
    figure;
    for i=1:length(Sigmma(:,1))
        plot(t,Sigmma(i,:),'r');
        hold on;
    end
    figure;
    plot(t,talle);
end
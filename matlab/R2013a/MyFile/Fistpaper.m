function Fistpaper(L1,L2,L3,tn,dt)
    t=0:dt:tn;       %迭代时间序列
    nodenum=length(L1);    nodedim=3;    timedim=length(t);
    tanh=0.5;     %控制强度
    P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %马氏链概率转移矩阵
    LL=cat(3,L1,L2,L3);  Gama=eye(nodedim);      %外耦合结构矩阵集 %内耦合结构矩阵   
    u=[1];     L=L1;    %马氏链初始状态     %马氏链处于初始状态下的外耦合矩阵
    D=zeros(nodenum,nodenum);    pin=round(100*rand(1));    D(pin,pin)=1;
    c=0.8;    a1=0.3649;    a=0.02;    b=0.2;    p=0.8;      %随机发生耦合的概率
    trinum=zeros(nodenum,1);    %记录各个节点的激发次数
    X=zeros(nodenum*nodedim,timedim);     X(:,1)=2*(2*rand(nodenum*nodedim,1)-0.5);
    s=zeros(nodedim,timedim);  s(:,1)=[1,-1,2]';    S=kron(ones(nodenum,1),s(:,1));       %目标函数的初值
    Xk=kron(ones(1,nodenum),X(:,1));     Sk=kron(ones(1,nodenum),S(:,1)); 
    HXk=H(Xk);    HSk=H(Sk);
    LHXk=diagbrock(kron(L,Gama)*HXk,nodedim,1);    DHXSk=diagbrock(kron(D,Gama)*(HXk-HSk),nodedim,1);
    for k=1:timedim-1
        sita_t=binornd(1,p);         %伯努利随机变量描述随机耦合
        c_t=2*c*rand(1);     %随机耦合强度为[0,2c]上的均匀分布
        X(:,k+1)=X(:,k)+(F(X(:,k))-sita_t*c_t*LHXk-tanh*c_t*DHXSk)*dt;      %迭代方程
        F(X(1:3,k))
        HXt=H(X(:,k));
        R=sita_t*(kron(L,Gama)*HXt-LHXk)-tanh*(kron(D,Gama)*(HXt-S)-DHXSk);  
        s(:,k+1)=s(:,k)+f(s(:,k))*dt;   %目标函数的迭代
        S=kron(ones(nodenum,1),s(:,k+1));
        E(:,k+1)=X(:,k+1)-S;       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
        Rt=reshape(R,nodenum,nodedim);
        Et=reshape(E(:,k+1),nodenum,nodedim);
        triggeif=sqrt(sum(Rt'.^2))>a1*sqrt(sum(Et'.^2));
        trinum=trinum+triggeif';
        if sum(triggeif)>0
            HXk(:,triggeif)=H(kron(ones(1,sum(triggeif)),X(:,k)));   
            HSk(:,triggeif)=H(kron(ones(1,sum(triggeif)),S)); 
            LHXk=diagbrock(kron(L,Gama)*HXk,nodedim,1);    
            DHXSk=diagbrock(kron(D,Gama)*(HXk-HSk),nodedim,1);
        end
        u(k+1)=Markov(u(k),P);     %u为马氏链的轨道
        L=LL(:,:,u(k+1));      %马氏链调制耦合结构 
        D(pin,pin)=0;
        pin=round(100*rand(1));
        if pin==0
            pin=1;
        end
        D(pin,pin)=1;
    end
    figure;
    for j=1:nodedim
        for i=1:nodenum
            plot(t,X(nodedim*(j-1)+i,:));
            hold on;
        end
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
function B=diagbrock(A,a,b)
    len=length(A);
    B=zeros(len,b);
    for k=1:len/a
        B((k-1)*a+1:k*a,:)=A((k-1)*a+1:k*a,(k-1)*b+1:k*b);
    end
end %取矩阵的对角块的元素，对角块的形式为a行b列
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
            
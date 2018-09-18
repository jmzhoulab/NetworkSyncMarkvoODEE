function Copy_of_MarkvoODEEolder()
tic;
nodenum=5;      %节点数
nodedim=3;      %节点维数
dt=0.001;        %迭代步长
t=0:dt:0.5;      %迭代时间序列
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
inititspan=[3 1 -1;1 2 -5;1 -4 -1;-1 6 2;8 -3 1];        %各个节点状态的初始值
tanh=3.5;     %控制强度
L1 =[2 -1 0 -1 0;
    -1 2  -1 0 0;
     0 -1  2 -1 0;
    -1 0 -1 3 -1;
     0  0 0 -1 1];
L2=[ 2 -1  0  0 -1;...      %外耦合结构矩阵
    -1  3 -1 -1  0;...
     0 -1  1  0  0;...
     0 -1  0  2 -1;...
    -1  0  0 -1  2];
 L3=[ 1 0 0 -1 0;...      %外耦合结构矩阵
      0 1 0 0 -1;...7
      0 0 1 -1 0;...
     -1 0 -1 3 -1;...
      0 -1 0 -1 2];
L_dim=length(L1);        %外耦合结构矩阵的维数
D1=diag([0 1 0 1 1]);D2=diag([1 0 1 0 1]);D3=diag([1 1 1 0 0]);     %三种不同的控制方式，如D1表示只控制节点1、2
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %马氏链概率转移矩阵
LL=cat(3,L1,L2,L3);     %%外耦合结构矩阵集
DD=cat(3,D1,D2,D3);     %%控制矩阵集
Gama=eye(nodedim);      %内耦合结构矩阵
u=[1];     %马氏链初始状态
L=L1;      %马氏链处于初始状态下的外耦合矩阵
D=D1;      %马氏链处于初始状态下的控制节点集矩阵
c=7;        %耦合强度均值
a1=0.3649;
a=0.2;
b=0.2;
p=0.8;      %随机发生耦合的概率
s=[1,-1,2];         %目标函数的初值
tk=ones(nodenum,1);         %节点的激发时刻
trinum=zeros(nodenum,1);        %用于记录每个节点激发次数
timedim=length(t);      %迭代次数
x=zeros(nodenum,nodedim,timedim);       %生成一个三维的矩阵，用于存放数据，第一维度表示节点，第二维度表示节点状态，第三维度表示时间
w=zeros(nodenum,nodedim,timedim);       %用于存放wi
for i=1:timedim     %把时间t放在数组x的第三维中
    x(:,:,i)=t(i);
end
x(:,:,1)=inititspan;      %赋予初值
e=zeros(nodenum,nodedim,timedim);   %存放误差
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %初始误差
for k=1:timedim-1       %时间方向上的循环
    sita_t=binornd(1,p);         %伯努利随机变量描述随机耦合
    c_t=2*c*rand(1);     %随机耦合强度为[0,2c]上的均匀分布
    for i=1:nodenum
        cp=zeros(1,nodedim);
        for j=1:L_dim        %计算系统中的求和项，即耦合项的和
            q=L(i,j)*g(x(j,:,tk(i)))*Gama; 
            cp=cp+q;
        end
        x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))-sita_t*c_t*cp-tanh*c_t*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama')*dt;   %迭代方程
        %{
        %2、连续情形下的更新规则
        wsum=zeros(1,nodedim);
        for j=1:nodenum    %计算wi(t)中的求和项
            wsum=wsum+L(i,j)*(g(x(j,:,k))-g(x(j,:,tk(i))))*Gama';
        end
        w(i,:,k)=sita_t*wsum-tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:))+g(s(k,:))-g(x(i,:,k)))*Gama';     %wi(t)的表达式
        %if norm(w(i,:,k))>a1*norm(e(i,:,k))    %判断是否满足促发条件1，其中a为常数
        if norm(w(i,:,k))>a*exp(-b*t(k))    %判断是否满足促发条件2，其中a,b为常数
            tk(i)=k;    %更新时间
            trinum(i)=trinum(i)+1;    %记录节点i激发次数
        end
        %}
        %
        %3、离散情形下的更新规则
            %%更新时间的预测
        zsum=zeros(1,nodedim);
        for j=1:nodenum    %计算zi(t)中的求和项
            zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
        end
        z=zeros(nodenum,nodedim,timedim);   %系统中除了f后剩下的所有项
        z(i,:,tk(i))=-sita_t*c_t*zsum-tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama'; 
        ysum=0;y=dt;tkj=tk;yi=0;
        while 1
            for j=1:nodenum
                if j==i
                    continue
                end
                for r=1:tk(j)
                    if t(r)>y+t(tk(i))
                        tkj(j)=r-1;
                    end
                end
                zj=zeros(nodenum,nodedim,timedim);   %系统中除了f后剩下的所有项
                zj(j,:,tkj(j))=-sita_t*c_t*zsum-tanh*D(j,j)*(g(x(j,:,tkj(j)))-g(s(tkj(j),:)))*Gama';
                ysum=ysum-L(i,j)*phi(y,z(i,:,tk(i)),zj(j,:,tkj(j)),x(i,:,tk(i)),x(j,:,tk(i)));
             end
             upwi=ysum+tanh*D(i,i)*phi(y,z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi更新规则中不等号左边的项
             %if upwi>a1*pci(y,z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %误差项控制跟新规则       
             if upwi>a*exp(-b*(y+t(tk(i))))        %指数控制跟新规则
                yi=y-dt;
                break;      %跳出while循环体
             end
             y=y+dt;
             if y>10*dt
                 break;
             end
        end
        for j=1:nodenum       %判断邻居节点是否有更新
            if j==i
                continue;
            end
            if t(k)==t(tk(j))
                judge=1;
            else
                judge=0;
            end
        end
        if yi==dt       %预测区间刚好等于下一个时间点时就发生更新
           tk(i)=k;
           trinum(i)=trinum(i)+1;    %记录节点i激发次数
        elseif judge==1     %邻居节点发生更新时也要更新
               tk(i)=k;
               trinum(i)=trinum(i)+1;    %记录节点i激发次数
        elseif k>1      
            if u(k)~=u(k-1)     %如果马氏链发生切换是也发生更新
               tk(i)=k;
               trinum(i)=trinum(i)+1;    %记录节点i激发次数
            end 
        else
        end 
        %离散情形下更新规则结束
        %}
    end
    s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %目标函数的迭代
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
    u(k+1)=Markov(u(k),P);     %u为马氏链的轨道
    L=LL(:,:,u(k+1));      %马氏链调制耦合结构 
    D=DD(:,:,u(k+1));      %马氏链调制控制节点集 
      %1、无更新规则连续控制
        %tk=(k+1).*ones(nodenum,1);
end

E=zeros(nodenum*nodedim,timedim);       %总误差向量，按列排
for k=1:timedim
    E(:,k)=[e(1,:,k),e(2,:,k),e(3,:,k),e(4,:,k),e(5,:,k)]';
end
%李雅普诺夫函数
V=zeros(timedim,1);
for k=1:timedim
    V(k)=1/2*E(:,k)'*E(:,k);
end

x1=zeros(nodenum,1);x2=zeros(nodenum,1);x3=zeros(nodenum,1);
for k=1:timedim
   x1(:,k)=x(:,1,k);       %节点状态1的数据，其中行表示节点，列表示时间
   x2(:,k)=x(:,2,k);       %节点状态2的数据，其中行表示节点，列表示时间
   x3(:,k)=x(:,3,k);       %节点状态3的数据，其中行表示节点，列表示时间
end
e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
for k=1:timedim
   e1(:,k)=e(:,1,k);       %节点状态1的误差，其中行表示节点，列表示时间
   e2(:,k)=e(:,2,k);       %节点状态2的误差，其中行表示节点，列表示时间
   e3(:,k)=e(:,3,k);       %节点状态3的误差，其中行表示节点，列表示时间
end
ee1=zeros(1,timedim);
ee2=zeros(1,timedim);
ee3=zeros(1,timedim);
for i=1:nodenum
    ee1=e1(i,:).^2+ee1;
    ee2=e2(i,:).^2+ee2;
    ee3=e3(i,:).^2+ee3;
end
ee=((ee1+ee2+ee3)).^(1/2);

DEtime=trinum;
DEV=V;

figure;
   %状态1的轨迹图
    plot(t,x1(1,:),'b');
    hold on;
    plot(t,x2(1,:),'g');
    hold on;  
    plot(t,x3(1,:),'r');
    hold on;
for i=2:nodenum
    plot(t,x1(i,:),'b');
    hold on;
    
end
    %状态2的轨迹图
for i=2:nodenum
    plot(t,x2(i,:),'g');
    hold on;
end
   %状态3的轨迹图
for i=2:nodenum
    plot(t,x3(i,:),'r');
    hold on;
end
legend('x_{i}^{1},the state 1 trajectory of nodes','x_{i}^{2},the state 2 trajectory of nodes','x_{i}^{3},the state 3 trajectory of nodes')
xlabel('t');
ylabel('x_{i}^{k}, i=1,……,5');


figure;    
               %状态1的误差图
    plot(t,e1(1,:),'b');
    hold on;
    plot(t,e2(1,:),'g');
    hold on;  
    plot(t,e3(1,:),'r');
    hold on;          
for i=2:nodenum
    plot(t,e1(i,:));
    hold on;
end

                 %状态2的误差图
for i=2:nodenum
    plot(t,e2(i,:),'g');
    hold on;
end
                %状态3的误差图
for i=2:nodenum
    plot(t,e3(i,:),'r');
    hold on;
end
legend('e_{i}^{1},the state 1 error of nodes','e_{i}^{2},the state 2 error of nodes','e_{i}^{3}，the state 3 error of nodes')
xlabel('t');
ylabel('error');
%{
figure;     %李雅普诺夫函数图形
plot(t,V);
xlabel('t');
ylabel('V(t)');
%}

figure;        %总体误差图
plot(t,ee);
xlabel('t');
ylabel('E(t)');
axis([0 0.5 -1 15]);
%{
T =

        1986         141        4071        4013
        2286         179        4138        4018
        2009         167        4088        4013
        2056         207        4099        4018
        2407         209        4331        4023
W =

         340         176        4131        4022
         370         227        4103        4018
         380         180        4115        4014
         396         266        4112        4022
         328         238        4145        4036
%}
%{
tt=t(1:200:timedim);
CNVV=CNV(1:200:timedim);
CEVV=CEV(1:200:timedim);
DNVV=DNV(1:200:timedim);
DEVV=DEV(1:200:timedim);
noncontrolVV=noncontrolV(1:200:timedim);
figure
plot(tt,CNVV,'r-')     %Continuous monitoring under rule(9)
hold on
plot(tt,CEVV,'g.')     %Continuous monitoring under rule(10)
hold on
plot(tt,DNVV,'co')    %Discrete monitoring under rule(40)
hold on
plot(tt,DEVV,'b+')     %Discrete monitoring under rule(41)
hold on
plot(tt,noncontrolVV,'k-.')    %Continuous updating
hold on
plot(tt,exp(-2*tt),'c*')   %exp(-2t)
%}
    runtime=toc;
    display(strcat('计算时间：',num2str(runtime),'秒'));
end
function fx=f(x)
   %{
    p=1.8;
    q=2;
    m0=-1.5;
    m1=-0.5;
    %}
    p=9.78;
    q=14.97;
    m0=-1.31;
    m1=-0.75;
    
    fx = x;
    fx(1)= p*(-x(1)+x(2)-m1*x(1)+1/2*(m0-m1)*(abs(x(1)+1)-abs(x(1)-1)));   
    fx(2)= x(1)-x(2)+x(3);
    fx(3)= -q*x(2);
    %f(x)的雅可比矩阵
    J1=[-p*(1+m0) p 0;
        1        -1 1;
        0        -q 0];
    J2=[-p*(1+m1) p 0;
        1        -1 1;
        0        -q 0];
end
function  gx=g(x)
 gx=zeros(1,3);
 gx(1)=2*x(1)+0.2*sin(x(1));
 gx(2)=2*x(2)+1;
 gx(3)=3*x(3)+0.5*cos(2*x(3));
end
function a=pci(t,x,y,u0,v0)
    v=1;
    l1=3.0608;
    a1=norm(x-y)/(-v*(2*l1^0.5+v));
    a=(norm(u0-v0)*exp(-2*l1^0.5-v)*t-a1*(exp(-(2*l1^0.5-v)*t)-1))^0.5;
end
function a=phi(t,x,y,u0,v0)
    l1=3.0608;
    l2=3.4438;
    beta=3.25;
    eta1=(norm(x-y)+l1*norm(u0-v0))/l1;
    eta2=(norm(x-y)+l2*norm(u0-v0))/l2;
    a=beta*sqrt(eta1^2*(exp(l1*t)-1)^2+eta2^2*(exp(l2*t)-1)^2);
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
    %v=[v,random,min,max]
end

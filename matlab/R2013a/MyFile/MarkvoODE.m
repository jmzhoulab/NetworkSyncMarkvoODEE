nodenum=8;      %节点数
nodedim=3;      %节点维数
dt=0.01;        %迭代步长
t=0:dt:3;      %迭代时间序列
inititspan=[3 1 -1;1 2 1;1 -1 -1;-1 1 2;2 -3 1;2 -1 3;1 1 -2;2 -1 -1];        %各个节点状态的初始值
L1=[ -2 0 0 1 0 0 1 0;...      %外耦合结构矩阵
     1 -3 1 0 0 1 0 0;...
     1 1 -4 0 1 0 0 1;...
     0 1 0 -3 1 1 0 0;...
     1 0 0 0 -1 0 0 0;...
     0 1 0 0 1 -4 1 1;...
     1 0 0 1 0 1 -2 0;...
     1 0 1 0 0 1 0 -3];
 L2=[-3 1 0 1 0 0 0 1;...      %外耦合结构矩阵
     0 -3 1 1 0 1 0 0;...
     1 0 -2 0 1 0 0 0;...
     0 1 0 -3 0 1 0 1;...
     1 0 1 0 -4 0 1 1;...
     0 1 0 0 1 -2 0 0;...
     1 0 0 0 0 1 -2 0;...
     0 0 1 1 0 1 0 -3];
 L3=[-4 1 0 1 0 1 0 1;...      %外耦合结构矩阵
     0 -2 1 0 0 0 1 0;...
     1 0 -2 0 0 1 0 0;...
     1 1 0 -4 0 1 0 1;...
     0 0 0 0 -2 0 1 1;...
     0 1 0 0 0 -1 0 0;...
     1 0 1 0 1 0 -3 0;...
     1 0 0 0 0 1 0 -2];
 L4=[-3 0 0 1 0 1 0 1;...      %外耦合结构矩阵
     0 -2 0 0 1 0 1 0;...
     1 0 -3 1 0 0 0 1;...
     0 1 1 -4 1 1 0 0;...
     0 1 0 1 -4 0 1 1;...
     1 0 1 0 0 -3 1 0;...
     0 0 0 0 1 0 -1 0;...
     0 1 1 0 0 0 0 -2];
L_dim=length(L1);        %外耦合结构矩阵的维数
D1=diag([0 1 0 0 0 1 0 0]);D2=diag([1 0 1 0 0 1 0 0]);D3=diag([0 1 0 1 0 0 0 1]);D4=diag([0 0 0 1 1 0 1 0]);     %三种不同的控制方式，如D1表示只控制节点2、6
P=[0.2 0.4 0.3 0.1;0.3 0.2 0.2 0.3;0.1 0.3 0.4 0.2;0.4 0.2 0.2 0.2];       %马氏链概率转移矩阵
LL=cat(4,L1,L2,L3,L4);     %%外耦合结构矩阵集
DD=cat(4,D1,D2,D3,D4);     %%控制矩阵集
Gama=eye(nodedim);      %内耦合结构矩阵
u=[1];     %马氏链初始状态
L=L1;      %马氏链处于初始状态下的外耦合矩阵
D=D1;      %马氏链处于初始状态下的控制节点集矩阵
c=1;        %耦合强度均值
tanh=1;      %控制强度
a=0.5;
b=0.3;
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
e=zeros(nodenum,nodedim,timedim);       %定义误差变量
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
        wsum=zeros(1,nodedim);
        for j=1:nodenum    %计算wi(t)中的求和项
            wsum=wsum+L(i,j)*(g(x(j,:,k))-g(x(j,:,tk(i))))*Gama';
        end
        w(i,:,k)=sita_t*wsum+tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:))+g(s(k,:))-g(x(i,:,k)));     %wi(t)的表达式
        %if norm(w(i,:,k))>a*norm(e(i,:,k))    %判断是否满足促发条件1，其中a为常数
        if norm(w(i,:,k))>a*exp(-b*x(1,1,k))    %判断是否满足促发条件2，其中a,b为常数
            tk(i)=k;    %更新时间
            trinum(i)=trinum(i)+1;    %记录节点i激发次数
        end
    end
    s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %目标函数的迭代
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
    
    u(k+1)=Markov(u(k),P);     %u为马氏链的轨道
    L=LL(:,:,u(k+1));      %马氏链调制耦合结构 
    D=DD(:,:,u(k+1));      %马氏链调制控制节点集 
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
figure;
subplot(3,1,1);     %状态1的轨迹图
for i=1:nodenum
    plot(t,x1(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{1}, i=1,……,5');
 
subplot(3,1,2);     %状态2的轨迹图
for i=1:nodenum
    plot(t,x2(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{2}, i=1,……,5');

subplot(3,1,3);     %状态3的轨迹图
for i=1:nodenum
    plot(t,x3(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{3}, i=1,……,5');


figure;    
subplot(3,1,1);     %状态1的误差图
for i=1:nodenum
    plot(t,e1(i,:));
    hold on;
end
ylabel('{e_{i}^{1}}, i=1,……,5');

subplot(3,1,2);      %状态2的误差图
for i=1:nodenum
    plot(t,e2(i,:));
    hold on;
end
ylabel('{e_{i}^{2}}, i=1,……,5');

subplot(3,1,3);     %状态3的误差图
for i=1:nodenum
    plot(t,e3(i,:));
    hold on;
end
xlabel('t');
ylabel('e_{i}^{3}, i=1,……,5');

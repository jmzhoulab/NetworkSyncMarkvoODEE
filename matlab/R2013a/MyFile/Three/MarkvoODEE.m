%function sol =MarkvoODEE
nodenum=5;      %节点数
nodedim=3;      %节点维数
dt=0.001;        %迭代步长
t=0:dt:5;      %迭代时间序列
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
inititspan=[3 1 -1;1 2 -5;-5 -4 -1;-1 6 2;3 -3 1];        %各个节点状态的初始值
tanh=5;     %控制强度
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
      0 1 0 0 -1;...
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
c=2;        %耦合强度均值
eplong=0.02;
w=1.3;
p=1;      %随机发生耦合的概率
s=[1,-1,2];         %目标函数的初值
tk=ones(nodenum,1);         %节点的激发时刻
trinum=zeros(nodenum,1);        %用于记录每个节点激发次数
timedim=length(t);      %迭代次数
x=zeros(nodenum,nodedim,timedim);       %生成一个三维的矩阵，用于存放数据，第一维度表示节点，第二维度表示节点状态，第三维度表示时间
ewan=zeros(nodenum,nodedim,timedim);   %存放测量误差
dW=zeros(nodenum,nodedim,timedim);
gc=zeros(1,timedim);
gi=zeros(nodenum,timedim);
for i=1:timedim     %把时间t放在数组x的第三维中
    x(:,:,i)=t(i);
end
x(:,:,1)=inititspan;      %赋予初值
e=zeros(nodenum,nodedim,timedim);   %存放误差
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %初始误差
for k=1:timedim-1       %时间方向上的循环
    sita_t=binornd(1,p);         %伯努利随机变量描述随机耦合
    c_t=2*c*rand(1);     %随机耦合强度为[0,2c]上的均匀分布
    dw=0.05*randn(nodenum,nodedim);
    dW(:,:,k) =dw;
    for i=1:nodenum
        cp=zeros(1,nodedim);
        ewan(:,:,k)=x(:,:,tk(i))-x(:,:,k);
        for j=1:L_dim        %计算系统中的求和项，即耦合项的和
            q=L(i,j)*g(x(j,:,tk(i)))*Gama; 
            cp=cp+q;
        end
        x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))-sita_t*c_t*cp-tanh*c_t*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama')*dt+exp(-t(k))*dw(i,:);   %迭代方程

        %
        %1、分散式更新规则
        swan=s(k,:)-s(tk(i));
        eik=e(i,:,k);
        ewanik=ewan(i,:,k);
        gi(i,k)=(1+eplong)*(norm(ewanik(:)))^2+eplong*(norm(swan(:)))^2-w*(norm(eik(:)))^2;
        if gi(i,k)>0    %判断是否满足促发条件2，其中a,b为常数
           tk(i)=k;    %更新时间
           trinum(i)=trinum(i)+1;    %记录节点i激发次数
           gi(i,k)=0;
        end
        %}
    end
        %{
        %2、集中式更新规则
        swan=s(k,:)-s(tk(1));
        ek=e(:,:,k);
        ewank=ewan(:,:,k);
        gc(k)=(1+eplong)*(norm(ewank(:)))^2+eplong*(norm(kron(ones(nodenum,1),swan(:))))^2-w*(norm(ek(:)))^2;
        if gc(k)>0    %判断是否满足促发条件2，其中a,b为常数 
            tk=k*ones(nodenum,1);   %更新时间
            trinum=trinum+1;    %记录节点i激发次数
            gc(k)=0;
        end
        %}
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
w1=zeros(nodenum,1);w2=zeros(nodenum,1);w3=zeros(nodenum,1);
e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
s1=zeros(timedim,1);s2=zeros(timedim,1);s3=zeros(timedim,1);
for k=1:timedim
   x1(:,k)=x(:,1,k);       %节点状态1的数据，其中行表示节点，列表示时间
   x2(:,k)=x(:,2,k);       %节点状态2的数据，其中行表示节点，列表示时间
   x3(:,k)=x(:,3,k);       %节点状态3的数据，其中行表示节点，列表示时间
   w1(:,k)=dW(:,1,k);      %zaosheng
   w2(:,k)=dW(:,2,k);      
   w3(:,k)=dW(:,3,k);    
   s1(k)=s(k,1);      %zaosheng
   s2(k)=s(k,2);      
   s3(k)=s(k,3); 
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
legend('x_{i}^{1}(t)','x_{i}^{2}(t)','x_{i}^{3}(t)')
xlabel('t');
ylabel('x_{i}^{k}, i=1,……,5');
figure;
   %状态1的轨迹图
    plot(t,w1(1,:),'b');
    hold on;
    plot(t,w2(1,:),'g');
    hold on;  
    plot(t,w3(1,:),'r');
    hold on;
for i=2:nodenum
    plot(t,w1(i,:),'b');
    hold on;
    
end
    %状态2的轨迹图
% for i=2:nodenum
%     plot(t,w2(i,:),'g');
%     hold on;
% end
%    %状态3的轨迹图
% for i=2:nodenum
%     plot(t,w3(i,:),'r');
%     hold on;
% end
legend('w_{i}^{1}(t)','w_{i}^{2}(t)','w_{i}^{3}(t)')
xlabel('t');
ylabel('w_{i}^{k}, i=1,……,5');
trinum
       
figure;        %总体误差图
plot(t,ee);
xlabel('t');
ylabel('E(t)');
axis([0 5 -1 15]);
figure;        %gc
plot(t,gc);
xlabel('t');
ylabel('g^{c}(t)');
figure;
    plot(t,gi(1,:),'b');
    hold on;
    plot(t,gi(2,:),'g');
    hold on;
    plot(t,gi(3,:),'r');
    hold on;
    plot(t,gi(4,:),'c');
    hold on;
    plot(t,gi(5,:),'k');
    xlabel('t');
    ylabel('g^{i}(t)');

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
%end


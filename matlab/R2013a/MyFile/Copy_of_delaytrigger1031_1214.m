nodenum=5;      %节点数
nodedim=3;      %节点维数
dt=0.01;        %迭代步长
t=-0.02:dt:1;      %迭代时间序列
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
%inititspan=[3 1 -1;8 2 -5;1 -4 -3;-1 6 -3;8 -3 4];        %各个节点状态的初始值
inititspan=10*(rand(5,3)-0.5);
tau=20;    %5.2;     %控制强度
beta=0.8;
c=8.2;   %2.2;        %耦合强度
alpha1=14.0228;
alpha2=0.9530; 
L1 =[2 0 -1  0 -1;
        0 3  -1 -1 0;
       -1 -1  3 -1 0;
        0  -1 -1 3 -1;
       -1  0 0 -1 2];
L2=[ 2  -1   0   0 -1;...      %外耦合结构矩阵
       -1   2   0  -1  0;...
        0  -1   3  -1  -1;...
        0   0  -1   2  -1;...
       -1   0  -1  -1  3];
  L3=[2 -1  0 -1  0;...      %外耦合结构矩阵
       -1  3 -1  0 -1;...
        0 -1  2 -1  0;...
       -1  0 -1  3 -1;...
        0 -1  0 -1  2];
L_dim=length(L1);        %外耦合结构矩阵的维数
D=diag(binornd(1,beta,5,1));    %三种不同的控制方式
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %马氏链概率转移矩阵
Q=[-8 4 4;5 -8 3;1 7 -8]; 
S1=[0 0 3;1 2 3;0 2 0];
S2=[1 2 0;0 0 0;1 0 3];
a=[3.618 4.4812 5.133];
LL=cat(3,L1,L2,L3);     %%外耦合结构矩阵集
Gama=eye(nodedim);      %内耦合结构矩阵
gamma=1;
u=[1,2,1];     %马氏链初始状态
L=L1;      %马氏链处于初始状态下的外耦合矩阵
tk=1;       %控制激发时刻
trinum=0;        %用于记录激发次数
timedim=length(t);      %迭代次数
x=zeros(nodenum,nodedim,timedim);       %生成一个三维的矩阵，用于存放数据，第一维度表示节点，第二维度表示节点状态，第三维度表示时间
w=zeros(nodenum,nodedim,timedim);       %用于存放wi
x(:,:,1)=inititspan;      %赋予初值
x(:,:,2)=10*(rand(5,3)-0.5);      %赋予初值
x(:,:,3)=10*(rand(5,3)-0.5);      %赋予初值
e=zeros(nodenum,nodedim,timedim);   %存放同步误差
ewan=zeros(nodenum,nodedim,timedim);   %存放测量误差
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),1/nodenum*sum(x(:,:,1)));       %初始误差
g=zeros(length(t),1);
for k=1:3
    switch k
        case 2
            L=L2;
        otherwise
            L=L1;
    end
xhat=x(:,:,k)';
dehat=ewan(:,:,k)';
g(k)=norm(kron(L,eye(nodedim))*dehat(:))-0.35*norm(kron(L,eye(nodedim))*xhat(:));
end
for k=3:timedim-1       %时间方向上的循环
    D=diag(binornd(1,beta,5,1));    %三种不同的控制方式
        xhat=x(:,:,k)';
        dehat=ewan(:,:,k)';
        g(k)=norm(kron(L,eye(nodedim))*dehat(:))-0.35*norm(kron(L,eye(nodedim))*xhat(:));
        if g(k)>0
            tk=k;
            trinum=trinum+1;
        end
        ewan(:,:,k)=x(:,:,tk)-x(:,:,k);
    %end更新规则
    for i=1:nodenum
        cop=zeros(1,nodedim);
        for j=1:L_dim        %计算系统中的求和项，即耦合项的和
            q=L(i,j)*x(j,:,tk)*Gama; 
            cop=cop+q;
        end
        x(i,:,k+1)=x(i,:,k)+(delayf(x(i,:,k),x(i,:,k-2))-c*cop-tau*c*D(i,i)*(x(i,:,tk)-1/nodenum*sum(x(:,:,tk)))*Gama')*dt;   %迭代方程
    end
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),1/nodenum*sum(x(:,:,k+1)));      %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
    u(k+1)=Markov(u(k),P);     %u为马氏链的轨道
    L=LL(:,:,u(k+1));      %马氏链调制耦合结构 
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
ee=(ee1+ee2+ee3).^(1/2);

%测量误差
%节点在时刻k时候的向量范数
nodeNurm=zeros(nodenum,timedim);
for k=1:length(t)
    for i=1:nodenum
        nodeNurm(i,k)=norm(ewan(i,:,k));    %第i个节点在时刻k的范数
    end
end
DEtime=trinum;
DEV=V;
color='rgbkc';
figure;
subplot(3,1,1);
   %状态1的轨迹图
for i=1:nodenum
    plot(t(3:end),x1(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{1}');
subplot(3,1,2);
    %状态2的轨迹图
for i=1:nodenum
    plot(t(3:end),x2(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{2}');
subplot(3,1,3);
   %状态3的轨迹图
for i=1:nodenum
    plot(t(3:end),x3(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{3}');
xlabel('t');

e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
for k=1:timedim
   e1(:,k)=e(:,1,k);       %节点状态1的误差，其中行表示节点，列表示时间
   e2(:,k)=e(:,2,k);       %节点状态2的误差，其中行表示节点，列表示时间
   e3(:,k)=e(:,3,k);       %节点状态3的误差，其中行表示节点，列表示时间

end
figure;     %李雅普诺夫函数图形
semilogy(t(3:end),V(3:end),'b');
hold on;
semilogy(t(3:end),exp(-2*t(3:end)),'r');
xlabel('t');
legend('V(t)','e^{-2t}');

figure;        %总体误差图
plot(t(4:end),ee(4:end));
xlabel('t');
ylabel('E(t)');
axis([0 t(end) -1 10]);
  figure;
  for i=1:nodenum
    plot(t(3:end),nodeNurm(i,3:end))
    hold all;
  end
% set(0,'DefaultAxesLineStyleOrder','remove')
% set(0,'DefaultAxesColorOrder','remove')
% figure;
% plot(t,g,'b');
% xlabel('t');
% ylabel('g(t)');




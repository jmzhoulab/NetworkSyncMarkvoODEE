nodenum=5;      %节点数
nodedim=3;      %节点维数
dt=0.01;        %迭代步长
t=0:dt:1.5;      %迭代时间序列
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
inititspan=[3 1 -1;8 2 -5;1 -4 -3;-1 6 -3;8 -3 4];        %各个节点状态的初始值
tau=5.998;    %5.2;     %控制强度
c=2.3;   %2.2;        %耦合强度
alpha1=14.0228;
alpha2=0.9530; 
a0=0.518;   %a0>0.2127
L1 =[2 -1 0 -1 0;
    -1 2  -1 0 0;
     0 -1  2 -1 0;
    -1 0 -1 3 -1;
     0  0 0 -1 1];
L2=[ 2 -1  0  0 -1;...      %外耦合结构矩阵
    -1  2  0 -1  0;...
     0 -1  2  -1 0;...
     0  0 -1  2 -1;...
    -1  0  0 -1  2];
 L3=[ 1 0 0 -1 0;...      %外耦合结构矩阵
      0 1 0 0 -1;...
      0 0 1 -1 0;...
     -1 0 -1 3 -1;...
      0 -1 0 -1 2];
L_dim=length(L1);        %外耦合结构矩阵的维数
D1=diag([0 1 0 1 1]);D2=diag([1 0 1 0 1]);D3=diag([1 1 1 0 0]);     %三种不同的控制方式，如D1表示只控制节点1、2
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %马氏链概率转移矩阵
Q=[-8 4 4;5 -8 3;1 7 -8]; 
S1=[0 0 3;1 2 3;0 2 0];
S2=[1 2 0;0 0 0;1 0 3];
a=[3.618 4.4812 5.133];
LL=cat(3,L1,L2,L3);     %%外耦合结构矩阵集
DD=cat(3,D1,D2,D3);     %%控制矩阵集
Gama=eye(nodedim);      %内耦合结构矩阵
u=[1,2,1];     %马氏链初始状态
L=L1;      %马氏链处于初始状态下的外耦合矩阵
D=D1;      %马氏链处于初始状态下的控制节点集矩阵
tk=1;       %控制激发时刻
trinum=0;        %用于记录激发次数
timedim=length(t);      %迭代次数
x=zeros(nodenum,nodedim,timedim);       %生成一个三维的矩阵，用于存放数据，第一维度表示节点，第二维度表示节点状态，第三维度表示时间
w=zeros(nodenum,nodedim,timedim);       %用于存放wi
x(:,:,1)=inititspan;      %赋予初值
x(:,:,2)=inititspan*rand(1);      %赋予初值
x(:,:,3)=inititspan*rand(1);      %赋予初值
e=zeros(nodenum,nodedim,timedim);   %存放同步误差
ewan=zeros(nodenum,nodedim,timedim);   %存放测量误差
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),1/nodenum*sum(x(:,:,1)));       %初始误差
g=zeros(length(t),1);
k0=zeros(1,3);
k0(1)=c*gamma+tau*c*gamma*0.5858/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+4.1024);
k0(2)=c*gamma+tau*c*gamma*1/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+14.3145);
k0(3)=c*gamma+tau*c*gamma*1/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+14.3145);
for k=1:3
    switch k
        case 2
            L=L2;D=D2;
        otherwise
            L=L1;D=D1;
    end
xhat=x(:,:,k)';
ehat=ewan(:,:,k)';
E=kron((L+tau*D-tau*D*ones(nodenum,nodenum)/nodenum),Gama);
g(k)=norm(E*ehat(:))-k0(u(k))/c*norm(xhat(:));
end
for k=3:timedim-1       %时间方向上的循环
    %begin更新规则
        %{
        egDL=sort(eig(L*D));
        egL=[eig(L1),eig(L2),eig(L3)];
        lambda2=min(egL(2,:));
        for i=1:length(egDL)
            if egDL(i)>0
                lambda=egDL(i);
              break;
            end
        end
        lambdav=[max(eig(L1)),max(eig(L2)),max(eig(L3))];
        sm=0;
        for v=1:length(P)
            if ismember(v,S1(u(k),:))
            sm=sm+P(u(k),v)*(lambdav(v)-a(u(k)));
            end
        end
        piu=1/(2*lambda2^2)*sm;
        gamma=min(eig(Gama));
        k0=c*gamma+tau*c*gamma*lambda/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+piu);
        %}
        
        E=kron((L+tau*D-tau*D*ones(nodenum,nodenum)/nodenum),Gama);
        xhat=x(:,:,k)';
        ehat=ewan(:,:,k)';
        g(k)=norm(E*ehat(:))-k0(u(k))/c*norm(xhat(:));
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
    D=DD(:,:,u(k+1));      %马氏链调制控制节点集 
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
ee=(ee1.^(1/2)+ee2.^(1/2)+ee3.^(1/2));

%测量误差
%节点在时刻k时候的向量范数
nodeNurm=zeros(nodenum,timedim);
for k=1:length(t)
    for i=1:nodenum
        nodeNurm(i,k)=norm(ewan(i,:,k));    %第i个节点在时刻k的范数
    end
end

% ewan1=zeros(nodenum,1);ewan2=zeros(nodenum,1);ewan3=zeros(nodenum,1);
% for k=1:timedim
%    ewan1(:,k)=ewan(:,1,k);       %节点状态1的误差，其中行表示节点，列表示时间
%    ewan2(:,k)=ewan(:,2,k);       %节点状态2的误差，其中行表示节点，列表示时间
%    ewan3(:,k)=ewan(:,3,k);       %节点状态3的误差，其中行表示节点，列表示时间
% end
% eewan1=zeros(1,timedim);
% eewan2=zeros(1,timedim);
% eewan3=zeros(1,timedim);
% for i=1:nodenum
%     eewan1=ewan1(i,:).^2+eewan1;
%     eewan2=ewan2(i,:).^2+eewan2;
%     eewan3=ewan3(i,:).^2+eewan3;
% end

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

figure;     %李雅普诺夫函数图形
plot(t,V,'b');
hold on;
plot(t,exp(-2*t),'r');
xlabel('t');
ylabel('V(t)');

figure;        %总体误差图
plot(t,ee);
xlabel('t');
ylabel('E(t)');

%axis([0 t(end) -1 10]);
% %测量误差图
% eewan1v=eewan1(1:10:timedim);
% eewan2v=eewan2(1:10:timedim);
% eewan3v=eewan3(1:10:timedim);
% figure; plot(tv,eewan2v,'r');
% figure; plot(tv,eewan1v,'r');
% figure; plot(tv,eewan3v,'r');
set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
      'DefaultAxesLineStyleOrder','-|--|:')
 %各个节点的测量误差范数图
  figure;
  for i=1:nodenum
    plot(t,nodeNurm(i,:))
    hold all;
  end
set(0,'DefaultAxesLineStyleOrder','remove')
set(0,'DefaultAxesColorOrder','remove')
figure;
plot(t,g,'b');
xlabel('t');
ylabel('g(t)');




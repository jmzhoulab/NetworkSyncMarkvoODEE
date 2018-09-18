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
      0 1 0 0 -1;...
      0 0 1 -1 0;...
     -1 0 -1 3 -1;...
      0 -1 0 -1 2];
L_dim=length(L1);        %外耦合结构矩阵的维数
D1=diag([0 1 0 1 1]);D2=diag([1 0 1 0 1]);D3=diag([1 1 1 0 0]);     %三种不同的控制方式，如D1表示只控制节点1、2
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %马氏链概率转移矩阵
LL=cat(3,L1,L2,L3);     %%外耦合结构矩阵集
DD=cat(3,D1,D2,D3);     %%控制矩阵集
Gama=3*eye(nodedim);      %内耦合结构矩阵
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
        if norm(w(i,:,k))>a1*norm(e(i,:,k))    %判断是否满足促发条件1，其中a为常数
        %if norm(w(i,:,k))>a*exp(-b*t(k))    %判断是否满足促发条件2，其中a,b为常数
            tk(i)=k;    %更新时间
            trinum(i)=trinum(i)+1;    %记录节点i激发次数
        end
        %}
    end
    u(k+1)=Markov(u(k),P);     %u为马氏链的轨道  
    %
    %3、离散情形下的更新规则
         %%更新时间的预测
    if k==1|u(k+1)~=u(k)    %如果处于初始状态时以及马氏链切换时全部更新预测
    for i=1:nodenum
        zsum=zeros(1,nodedim);
        for j=1:nodenum    %计算zi(t)中的求和项
            zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
        end
        z=zeros(nodenum,nodedim,timedim);   %系统中除了f后剩下的所有项
        z(i,:,tk(i))=-sita_t*c_t*zsum-tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama'; 
        ysum=0;y=dt;tkj=tk;yy=zeros(nodenum,1);
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
             upwi=ysum+tanh*D(i,i)*phi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi更新规则中不等号左边的项
             %if upwi>a1*psi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %误差项控制跟新规则       
             if upwi>a*exp(-b*(y+t(tk(i))))        %指数控制跟新规则
                yi=y-dt;
                break;      %跳出while循环体
             end
             yy(i)=yy(i)+dt;
             if yy(i)>0.5
                 break;
             end
        end
    end
    end
    tri=zeros(nodenum,1);
    for i=1:nodenum
        if yy(i)==dt       %预测区间刚好等于下一个时间点时就发生更新
           tk(i)=k;
           trinum(i)=trinum(i)+1;    %记录节点i激发次数
           tri(i)=1;    
        end
    end
    for i=1:nodenum       %判断邻居节点是否有更新
        if tri(i)==1
           zsum=zeros(1,nodedim);
            for j=1:nodenum    %计算zi(t)中的求和项
                zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
            end
            z=zeros(nodenum,nodedim,timedim);   %系统中除了f后剩下的所有项
            z(i,:,tk(i))=-sita_t*c_t*zsum-tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama'; 
            ysum=0;y=dt;tkj=tk;yy=zeros(nodenum,1);
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
                upwi=ysum+tanh*D(i,i)*phi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi更新规则中不等号左边的项
                %if upwi>a1*psi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %误差项控制跟新规则       
                if upwi>a*exp(-b*(y+t(tk(i))))        %指数控制跟新规则
                    yi=y-dt;
                    break;      %跳出while循环体
                end
                yy(i)=yy(i)+dt;
                if yy(i)>0.5
                    break;
                end
            end
        end
    end   
        %离散情形下更新规则结束
        %}
    s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %目标函数的迭代
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %节点与目标的误差，这里所有的节点以及对应状态误差都一起计算
    L=LL(:,:,u(k+1));      %马氏链调制耦合结构 
    D=DD(:,:,u(k+1));      %马氏链调制控制节点集 
      %1、无更新规则连续控制
      %  tk=(k+1).*ones(nodenum,1);
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
%{
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
%}
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

figure;     %李雅普诺夫函数图形
plot(t,V);
xlabel('t');
ylabel('V(t)');


figure;        %总体误差图
plot(t,ee);
xlabel('t');
ylabel('E(t)');

%{
figure;
subplot(2,1,1);     %状态1的轨迹图
for i=1:nodenum
    plot(t,x1(i,:));
    hold on;
end
ylabel('x_{i}^{1}, i=1,……,5');
subplot(2,1,2);     %状态1的误差图
for i=1:nodenum
    plot(t,e1(i,:));
    hold on;
end
xlabel('t');
ylabel('{e_{i}^{1}}, i=1,……,5');

figure;
subplot(2,1,1);     %状态2的轨迹图
for i=1:nodenum
    plot(t,x2(i,:));
    hold on;
end
ylabel('x_{i}^{2}, i=1,……,5');
subplot(2,1,2);      %状态2的误差图
for i=1:nodenum
    plot(t,e2(i,:));
    hold on;
end
xlabel('t');
ylabel('{e_{i}^{2}}, i=1,……,5');

figure;
subplot(2,1,1);     %状态3的轨迹图
for i=1:nodenum
    plot(t,x3(i,:));
    hold on;
end
ylabel('x_{i}^{3}, i=1,……,5');
subplot(2,1,2);     %状态3的误差图
for i=1:nodenum
    plot(t,e3(i,:));
    hold on;
end
xlabel('t');
ylabel('e_{i}^{3}, i=1,……,5');
%}
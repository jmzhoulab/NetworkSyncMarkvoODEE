nodenum=5;      %�ڵ���
nodedim=3;      %�ڵ�ά��
dt=0.01;        %��������
t=-0.02:dt:1;      %����ʱ������
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
%inititspan=[3 1 -1;8 2 -5;1 -4 -3;-1 6 -3;8 -3 4];        %�����ڵ�״̬�ĳ�ʼֵ
inititspan=10*(rand(5,3)-0.5);
tau=20;    %5.2;     %����ǿ��
beta=0.8;
c=8.2;   %2.2;        %���ǿ��
alpha1=14.0228;
alpha2=0.9530; 
L1 =[2 0 -1  0 -1;
        0 3  -1 -1 0;
       -1 -1  3 -1 0;
        0  -1 -1 3 -1;
       -1  0 0 -1 2];
L2=[ 2  -1   0   0 -1;...      %����Ͻṹ����
       -1   2   0  -1  0;...
        0  -1   3  -1  -1;...
        0   0  -1   2  -1;...
       -1   0  -1  -1  3];
  L3=[2 -1  0 -1  0;...      %����Ͻṹ����
       -1  3 -1  0 -1;...
        0 -1  2 -1  0;...
       -1  0 -1  3 -1;...
        0 -1  0 -1  2];
L_dim=length(L1);        %����Ͻṹ�����ά��
D=diag(binornd(1,beta,5,1));    %���ֲ�ͬ�Ŀ��Ʒ�ʽ
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %����������ת�ƾ���
Q=[-8 4 4;5 -8 3;1 7 -8]; 
S1=[0 0 3;1 2 3;0 2 0];
S2=[1 2 0;0 0 0;1 0 3];
a=[3.618 4.4812 5.133];
LL=cat(3,L1,L2,L3);     %%����Ͻṹ����
Gama=eye(nodedim);      %����Ͻṹ����
gamma=1;
u=[1,2,1];     %��������ʼ״̬
L=L1;      %���������ڳ�ʼ״̬�µ�����Ͼ���
tk=1;       %���Ƽ���ʱ��
trinum=0;        %���ڼ�¼��������
timedim=length(t);      %��������
x=zeros(nodenum,nodedim,timedim);       %����һ����ά�ľ������ڴ�����ݣ���һά�ȱ�ʾ�ڵ㣬�ڶ�ά�ȱ�ʾ�ڵ�״̬������ά�ȱ�ʾʱ��
w=zeros(nodenum,nodedim,timedim);       %���ڴ��wi
x(:,:,1)=inititspan;      %�����ֵ
x(:,:,2)=10*(rand(5,3)-0.5);      %�����ֵ
x(:,:,3)=10*(rand(5,3)-0.5);      %�����ֵ
e=zeros(nodenum,nodedim,timedim);   %���ͬ�����
ewan=zeros(nodenum,nodedim,timedim);   %��Ų������
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),1/nodenum*sum(x(:,:,1)));       %��ʼ���
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
for k=3:timedim-1       %ʱ�䷽���ϵ�ѭ��
    D=diag(binornd(1,beta,5,1));    %���ֲ�ͬ�Ŀ��Ʒ�ʽ
        xhat=x(:,:,k)';
        dehat=ewan(:,:,k)';
        g(k)=norm(kron(L,eye(nodedim))*dehat(:))-0.35*norm(kron(L,eye(nodedim))*xhat(:));
        if g(k)>0
            tk=k;
            trinum=trinum+1;
        end
        ewan(:,:,k)=x(:,:,tk)-x(:,:,k);
    %end���¹���
    for i=1:nodenum
        cop=zeros(1,nodedim);
        for j=1:L_dim        %����ϵͳ�е������������ĺ�
            q=L(i,j)*x(j,:,tk)*Gama; 
            cop=cop+q;
        end
        x(i,:,k+1)=x(i,:,k)+(delayf(x(i,:,k),x(i,:,k-2))-c*cop-tau*c*D(i,i)*(x(i,:,tk)-1/nodenum*sum(x(:,:,tk)))*Gama')*dt;   %��������
    end
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),1/nodenum*sum(x(:,:,k+1)));      %�ڵ���Ŀ������������еĽڵ��Լ���Ӧ״̬��һ�����
    u(k+1)=Markov(u(k),P);     %uΪ�������Ĺ��
    L=LL(:,:,u(k+1));      %������������Ͻṹ 
end
E=zeros(nodenum*nodedim,timedim);       %�����������������
for k=1:timedim
    E(:,k)=[e(1,:,k),e(2,:,k),e(3,:,k),e(4,:,k),e(5,:,k)]';
end
%������ŵ����
V=zeros(timedim,1);
for k=1:timedim
    V(k)=1/2*E(:,k)'*E(:,k);
end

x1=zeros(nodenum,1);x2=zeros(nodenum,1);x3=zeros(nodenum,1);
for k=1:timedim
   x1(:,k)=x(:,1,k);       %�ڵ�״̬1�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x2(:,k)=x(:,2,k);       %�ڵ�״̬2�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x3(:,k)=x(:,3,k);       %�ڵ�״̬3�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
end
e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
for k=1:timedim
   e1(:,k)=e(:,1,k);       %�ڵ�״̬1���������б�ʾ�ڵ㣬�б�ʾʱ��
   e2(:,k)=e(:,2,k);       %�ڵ�״̬2���������б�ʾ�ڵ㣬�б�ʾʱ��
   e3(:,k)=e(:,3,k);       %�ڵ�״̬3���������б�ʾ�ڵ㣬�б�ʾʱ��
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

%�������
%�ڵ���ʱ��kʱ�����������
nodeNurm=zeros(nodenum,timedim);
for k=1:length(t)
    for i=1:nodenum
        nodeNurm(i,k)=norm(ewan(i,:,k));    %��i���ڵ���ʱ��k�ķ���
    end
end
DEtime=trinum;
DEV=V;
color='rgbkc';
figure;
subplot(3,1,1);
   %״̬1�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x1(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{1}');
subplot(3,1,2);
    %״̬2�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x2(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{2}');
subplot(3,1,3);
   %״̬3�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x3(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{3}');
xlabel('t');

e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
for k=1:timedim
   e1(:,k)=e(:,1,k);       %�ڵ�״̬1���������б�ʾ�ڵ㣬�б�ʾʱ��
   e2(:,k)=e(:,2,k);       %�ڵ�״̬2���������б�ʾ�ڵ㣬�б�ʾʱ��
   e3(:,k)=e(:,3,k);       %�ڵ�״̬3���������б�ʾ�ڵ㣬�б�ʾʱ��

end
figure;     %������ŵ����ͼ��
semilogy(t(3:end),V(3:end),'b');
hold on;
semilogy(t(3:end),exp(-2*t(3:end)),'r');
xlabel('t');
legend('V(t)','e^{-2t}');

figure;        %�������ͼ
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




nodenum=5;      %�ڵ���
nodedim=3;      %�ڵ�ά��
dt=0.01;        %��������
t=-0.02:dt:1;      %����ʱ������
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
%inititspan=[3 1 -1;8 2 -5;1 -4 -3;-1 6 -3;8 -3 4];        %�����ڵ�״̬�ĳ�ʼֵ
inititspan=10*(rand(5,3)-0.5);
tau=7;    %5.2;     %����ǿ��
c=3.85;   %2.2;        %���ǿ��
alpha1=14.0228;
alpha2=0.9530; 
a0=0.518;   %a0>0.2127
L1 =[2 -1 0 -1 0;
    -1 2  -1 0 0;
     0 -1  2 -1 0;
    -1 0 -1 3 -1;
     0  0 0 -1 1];
L2=[ 2 -1  0  0 -1;...      %����Ͻṹ����
    -1  2  0 -1  0;...
     0 -1  2  -1 0;...
     0  0 -1  2 -1;...
    -1  0  0 -1  2];
 L3=[ 1 0 0 -1 0;...      %����Ͻṹ����
      0 1 0 0 -1;...
      0 0 1 -1 0;...
     -1 0 -1 3 -1;...
      0 -1 0 -1 2];
L_dim=length(L1);        %����Ͻṹ�����ά��
D1=diag([0.2 0.8 0.1 0.7 0.8]);D2=diag([0.8 0.1 0.6 0.3 0.7]);D3=diag([0.8 0.6 0.6 0.1 0.8]);     %���ֲ�ͬ�Ŀ��Ʒ�ʽ����D1��ʾֻ���ƽڵ�1��2
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %����������ת�ƾ���
Q=[-8 4 4;5 -8 3;1 7 -8]; 
S1=[0 0 3;1 2 3;0 2 0];
S2=[1 2 0;0 0 0;1 0 3];
a=[3.618 4.4812 5.133];
LL=cat(3,L1,L2,L3);     %%����Ͻṹ����
DD=cat(3,D1,D2,D3);     %%���ƾ���
Gama=eye(nodedim);      %����Ͻṹ����
gamma=1;
u=[1,2,1];     %��������ʼ״̬
L=L1;      %���������ڳ�ʼ״̬�µ�����Ͼ���
D=D1;      %���������ڳ�ʼ״̬�µĿ��ƽڵ㼯����
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
k0=zeros(1,3);
egL=[eig(L1),eig(L2),eig(L3)];
lambda2=min(egL(2,:));
egDL=sort([eig(L1*D1),eig(L2*D2),eig(L3*D3)]);
egDL21=egDL(2,1);
egDL22=egDL(2,2);
egDL23=egDL(2,3);
k0(1)=c*gamma+tau*c*gamma*egDL21/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+4.1024);
k0(2)=c*gamma+tau*c*gamma*egDL22/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+ 11.0937);
k0(3)=c*gamma+tau*c*gamma*egDL23/lambda2^2-(a0/(2*lambda2)+alpha1/lambda2+alpha2/2+4.6234);
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
for k=3:timedim-1       %ʱ�䷽���ϵ�ѭ��
    %begin���¹���
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
        
        E=kron((L+tau*D+tau*D*ones(nodenum,nodenum)/nodenum),Gama);
        xhat=x(:,:,k)';
        ehat=ewan(:,:,k)';
        g(k)=norm(E*ehat(:))-k0(u(k))/c*norm(xhat(:));
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
    D=DD(:,:,u(k+1));      %���������ƿ��ƽڵ㼯 
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
e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
xx=zeros(10,timedim);
for k=1:timedim
   x1(:,k)=x(:,1,k);       %�ڵ�״̬1�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x2(:,k)=x(:,2,k);       %�ڵ�״̬2�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x3(:,k)=x(:,3,k);       %�ڵ�״̬3�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   e1(:,k)=e(:,1,k);       %�ڵ�״̬1���������б�ʾ�ڵ㣬�б�ʾʱ��
   e2(:,k)=e(:,2,k);       %�ڵ�״̬2���������б�ʾ�ڵ㣬�б�ʾʱ��
   e3(:,k)=e(:,3,k);       %�ڵ�״̬3���������б�ʾ�ڵ㣬�б�ʾʱ��
   lag=1;
    for i=1:5-1
        for j=i+1:5
            xx(lag,k)=norm(x(i,:,k)-x(j,:,k))^2;
            lag=lag+1;
        end
    end
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

% ewan1=zeros(nodenum,1);ewan2=zeros(nodenum,1);ewan3=zeros(nodenum,1);
% for k=1:timedim
%    ewan1(:,k)=ewan(:,1,k);       %�ڵ�״̬1���������б�ʾ�ڵ㣬�б�ʾʱ��
%    ewan2(:,k)=ewan(:,2,k);       %�ڵ�״̬2���������б�ʾ�ڵ㣬�б�ʾʱ��
%    ewan3(:,k)=ewan(:,3,k);       %�ڵ�״̬3���������б�ʾ�ڵ㣬�б�ʾʱ��
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
color='rgbkc';
figure;
subplot(3,1,1);
   %״̬1�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x1(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{1}');
%legend('x_{1}^{1}','x_{2}^{1}','x_{3}^{1}','x_{4}^{1}','x_{5}^{1}');
%xlabel('t');
subplot(3,1,2);
    %״̬2�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x2(i,3:end),strcat(color(i)));
    hold on;
end
ylabel('x_{i}^{2}');
%legend('x_{1}^{2}','x_{2}^{2}','x_{3}^{2}','x_{4}^{2}','x_{5}^{2}');
%xlabel('t');
subplot(3,1,3);
   %״̬3�Ĺ켣ͼ
for i=1:nodenum
    plot(t(3:end),x3(i,3:end),strcat(color(i)));
    hold on;
end
%legend('x_{1}^{3}','x_{2}^{3}','x_{3}^{3}','x_{4}^{3}','x_{5}^{3}');
ylabel('x_{i}^{3}');
xlabel('t');

e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
for k=1:timedim
   e1(:,k)=e(:,1,k);       %�ڵ�״̬1���������б�ʾ�ڵ㣬�б�ʾʱ��
   e2(:,k)=e(:,2,k);       %�ڵ�״̬2���������б�ʾ�ڵ㣬�б�ʾʱ��
   e3(:,k)=e(:,3,k);       %�ڵ�״̬3���������б�ʾ�ڵ㣬�б�ʾʱ��

end
% e1sum=zeros(timedim,1);
% e2sum=zeros(timedim,1);
% e3sum=zeros(timedim,1);
% for k=1:timedim
%     e1sum(k)=sqrt(sum(e1(:,k).^2));
%     e2sum(k)=sqrt(sum(e2(:,k).^2));
%     e3sum(k)=sqrt(sum(e3(:,k).^2));
% end
% figure;
% subplot(3,1,1);plot(t(3:end),e1sum(3:end));        %״̬1�����ͼ 
% subplot(3,1,2);plot(t(3:end),e2sum(3:end));     %״̬2�����
% subplot(3,1,3);plot(t(3:end),e3sum(3:end));     %״̬3�����ͼ
% legend('e_{i}^{1},the state 1 error of nodes','e_{i}^{2},the state 2 error of nodes','e_{i}^{3}��the state 3 error of nodes')
% xlabel('t');
% ylabel('');
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
% %�������ͼ
% eewan1v=eewan1(1:10:timedim);
% eewan2v=eewan2(1:10:timedim);
% eewan3v=eewan3(1:10:timedim);
% figure; plot(tv,eewan2v,'r');
% figure; plot(tv,eewan1v,'r');
% figure; plot(tv,eewan3v,'r');
% set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
%       'DefaultAxesLineStyleOrder','-|--|:')
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
% 
% e11=sqrt(sum(e1-1/5*ones(length(e1),1)*sum(e1)));
% e21=sqrt(sum(e2-1/5*ones(length(e2),1)*sum(e2)));
% e31=sqrt(sum(e3-1/5*ones(length(e3),1)*sum(e3)));
% figure;
% plot(t(3:end),e11(3:end),'r');
% hold on;
% plot(t(3:end),e21(3:end),'g');
% plot(t(3:end),e31(3:end),'b');
% xlabel('t');
figure;
plot(t(3:end),xx(:,3:end));
        



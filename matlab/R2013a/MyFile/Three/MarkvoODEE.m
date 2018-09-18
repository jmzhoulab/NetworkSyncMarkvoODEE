%function sol =MarkvoODEE
nodenum=5;      %�ڵ���
nodedim=3;      %�ڵ�ά��
dt=0.001;        %��������
t=0:dt:5;      %����ʱ������
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
inititspan=[3 1 -1;1 2 -5;-5 -4 -1;-1 6 2;3 -3 1];        %�����ڵ�״̬�ĳ�ʼֵ
tanh=5;     %����ǿ��
L1 =[2 -1 0 -1 0;
    -1 2  -1 0 0;
     0 -1  2 -1 0;
    -1 0 -1 3 -1;
     0  0 0 -1 1];
L2=[ 2 -1  0  0 -1;...      %����Ͻṹ����
    -1  3 -1 -1  0;...
     0 -1  1  0  0;...
     0 -1  0  2 -1;...
    -1  0  0 -1  2];
 L3=[ 1 0 0 -1 0;...      %����Ͻṹ����
      0 1 0 0 -1;...
      0 0 1 -1 0;...
     -1 0 -1 3 -1;...
      0 -1 0 -1 2];
L_dim=length(L1);        %����Ͻṹ�����ά��
D1=diag([0 1 0 1 1]);D2=diag([1 0 1 0 1]);D3=diag([1 1 1 0 0]);     %���ֲ�ͬ�Ŀ��Ʒ�ʽ����D1��ʾֻ���ƽڵ�1��2
P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %����������ת�ƾ���
LL=cat(3,L1,L2,L3);     %%����Ͻṹ����
DD=cat(3,D1,D2,D3);     %%���ƾ���
Gama=eye(nodedim);      %����Ͻṹ����
u=[1];     %��������ʼ״̬
L=L1;      %���������ڳ�ʼ״̬�µ�����Ͼ���
D=D1;      %���������ڳ�ʼ״̬�µĿ��ƽڵ㼯����
c=2;        %���ǿ�Ⱦ�ֵ
eplong=0.02;
w=1.3;
p=1;      %���������ϵĸ���
s=[1,-1,2];         %Ŀ�꺯���ĳ�ֵ
tk=ones(nodenum,1);         %�ڵ�ļ���ʱ��
trinum=zeros(nodenum,1);        %���ڼ�¼ÿ���ڵ㼤������
timedim=length(t);      %��������
x=zeros(nodenum,nodedim,timedim);       %����һ����ά�ľ������ڴ�����ݣ���һά�ȱ�ʾ�ڵ㣬�ڶ�ά�ȱ�ʾ�ڵ�״̬������ά�ȱ�ʾʱ��
ewan=zeros(nodenum,nodedim,timedim);   %��Ų������
dW=zeros(nodenum,nodedim,timedim);
gc=zeros(1,timedim);
gi=zeros(nodenum,timedim);
for i=1:timedim     %��ʱ��t��������x�ĵ���ά��
    x(:,:,i)=t(i);
end
x(:,:,1)=inititspan;      %�����ֵ
e=zeros(nodenum,nodedim,timedim);   %������
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %��ʼ���
for k=1:timedim-1       %ʱ�䷽���ϵ�ѭ��
    sita_t=binornd(1,p);         %��Ŭ�������������������
    c_t=2*c*rand(1);     %������ǿ��Ϊ[0,2c]�ϵľ��ȷֲ�
    dw=0.05*randn(nodenum,nodedim);
    dW(:,:,k) =dw;
    for i=1:nodenum
        cp=zeros(1,nodedim);
        ewan(:,:,k)=x(:,:,tk(i))-x(:,:,k);
        for j=1:L_dim        %����ϵͳ�е������������ĺ�
            q=L(i,j)*g(x(j,:,tk(i)))*Gama; 
            cp=cp+q;
        end
        x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))-sita_t*c_t*cp-tanh*c_t*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama')*dt+exp(-t(k))*dw(i,:);   %��������

        %
        %1����ɢʽ���¹���
        swan=s(k,:)-s(tk(i));
        eik=e(i,:,k);
        ewanik=ewan(i,:,k);
        gi(i,k)=(1+eplong)*(norm(ewanik(:)))^2+eplong*(norm(swan(:)))^2-w*(norm(eik(:)))^2;
        if gi(i,k)>0    %�ж��Ƿ�����ٷ�����2������a,bΪ����
           tk(i)=k;    %����ʱ��
           trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
           gi(i,k)=0;
        end
        %}
    end
        %{
        %2������ʽ���¹���
        swan=s(k,:)-s(tk(1));
        ek=e(:,:,k);
        ewank=ewan(:,:,k);
        gc(k)=(1+eplong)*(norm(ewank(:)))^2+eplong*(norm(kron(ones(nodenum,1),swan(:))))^2-w*(norm(ek(:)))^2;
        if gc(k)>0    %�ж��Ƿ�����ٷ�����2������a,bΪ���� 
            tk=k*ones(nodenum,1);   %����ʱ��
            trinum=trinum+1;    %��¼�ڵ�i��������
            gc(k)=0;
        end
        %}
    s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %Ŀ�꺯���ĵ���
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %�ڵ���Ŀ������������еĽڵ��Լ���Ӧ״̬��һ�����
    u(k+1)=Markov(u(k),P);     %uΪ�������Ĺ��
    L=LL(:,:,u(k+1));      %������������Ͻṹ 
    D=DD(:,:,u(k+1));      %���������ƿ��ƽڵ㼯 
      %1���޸��¹�����������
        %tk=(k+1).*ones(nodenum,1);
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
w1=zeros(nodenum,1);w2=zeros(nodenum,1);w3=zeros(nodenum,1);
e1=zeros(nodenum,1);e2=zeros(nodenum,1);e3=zeros(nodenum,1);
s1=zeros(timedim,1);s2=zeros(timedim,1);s3=zeros(timedim,1);
for k=1:timedim
   x1(:,k)=x(:,1,k);       %�ڵ�״̬1�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x2(:,k)=x(:,2,k);       %�ڵ�״̬2�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   x3(:,k)=x(:,3,k);       %�ڵ�״̬3�����ݣ������б�ʾ�ڵ㣬�б�ʾʱ��
   w1(:,k)=dW(:,1,k);      %zaosheng
   w2(:,k)=dW(:,2,k);      
   w3(:,k)=dW(:,3,k);    
   s1(k)=s(k,1);      %zaosheng
   s2(k)=s(k,2);      
   s3(k)=s(k,3); 
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
ee=((ee1+ee2+ee3)).^(1/2);

DEtime=trinum;
DEV=V;

figure;
   %״̬1�Ĺ켣ͼ
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
    %״̬2�Ĺ켣ͼ
for i=2:nodenum
    plot(t,x2(i,:),'g');
    hold on;
end
   %״̬3�Ĺ켣ͼ
for i=2:nodenum
    plot(t,x3(i,:),'r');
    hold on;
end
legend('x_{i}^{1}(t)','x_{i}^{2}(t)','x_{i}^{3}(t)')
xlabel('t');
ylabel('x_{i}^{k}, i=1,����,5');
figure;
   %״̬1�Ĺ켣ͼ
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
    %״̬2�Ĺ켣ͼ
% for i=2:nodenum
%     plot(t,w2(i,:),'g');
%     hold on;
% end
%    %״̬3�Ĺ켣ͼ
% for i=2:nodenum
%     plot(t,w3(i,:),'r');
%     hold on;
% end
legend('w_{i}^{1}(t)','w_{i}^{2}(t)','w_{i}^{3}(t)')
xlabel('t');
ylabel('w_{i}^{k}, i=1,����,5');
trinum
       
figure;        %�������ͼ
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


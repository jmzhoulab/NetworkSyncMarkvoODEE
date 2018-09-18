nodenum=5;      %�ڵ���
nodedim=3;      %�ڵ�ά��
dt=0.001;        %��������
t=0:dt:0.5;      %����ʱ������
%inititspan=[1 1 -1;1 1 -1;1 1 -1;1 1 -1;1 1 -1];  
inititspan=[3 1 -1;1 2 -5;1 -4 -1;-1 6 2;8 -3 1];        %�����ڵ�״̬�ĳ�ʼֵ
tanh=3.5;     %����ǿ��
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
Gama=3*eye(nodedim);      %����Ͻṹ����
u=[1];     %��������ʼ״̬
L=L1;      %���������ڳ�ʼ״̬�µ�����Ͼ���
D=D1;      %���������ڳ�ʼ״̬�µĿ��ƽڵ㼯����
c=7;        %���ǿ�Ⱦ�ֵ
a1=0.3649;
a=0.2;
b=0.2;
p=0.8;      %���������ϵĸ���
s=[1,-1,2];         %Ŀ�꺯���ĳ�ֵ
tk=ones(nodenum,1);         %�ڵ�ļ���ʱ��
trinum=zeros(nodenum,1);        %���ڼ�¼ÿ���ڵ㼤������
timedim=length(t);      %��������
x=zeros(nodenum,nodedim,timedim);       %����һ����ά�ľ������ڴ�����ݣ���һά�ȱ�ʾ�ڵ㣬�ڶ�ά�ȱ�ʾ�ڵ�״̬������ά�ȱ�ʾʱ��
w=zeros(nodenum,nodedim,timedim);       %���ڴ��wi
for i=1:timedim     %��ʱ��t��������x�ĵ���ά��
    x(:,:,i)=t(i);
end
x(:,:,1)=inititspan;      %�����ֵ
e=zeros(nodenum,nodedim,timedim);   %������
e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %��ʼ���
for k=1:timedim-1       %ʱ�䷽���ϵ�ѭ��
    sita_t=binornd(1,p);         %��Ŭ�������������������
    c_t=2*c*rand(1);     %������ǿ��Ϊ[0,2c]�ϵľ��ȷֲ�
    for i=1:nodenum
        cp=zeros(1,nodedim);
        for j=1:L_dim        %����ϵͳ�е������������ĺ�
            q=L(i,j)*g(x(j,:,tk(i)))*Gama; 
            cp=cp+q;
        end
        x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))-sita_t*c_t*cp-tanh*c_t*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:)))*Gama')*dt;   %��������
        %{
        %2�����������µĸ��¹���
        wsum=zeros(1,nodedim);
        for j=1:nodenum    %����wi(t)�е������
            wsum=wsum+L(i,j)*(g(x(j,:,k))-g(x(j,:,tk(i))))*Gama';
        end
        w(i,:,k)=sita_t*wsum-tanh*D(i,i)*(g(x(i,:,tk(i)))-g(s(tk(i),:))+g(s(k,:))-g(x(i,:,k)))*Gama';     %wi(t)�ı��ʽ
        if norm(w(i,:,k))>a1*norm(e(i,:,k))    %�ж��Ƿ�����ٷ�����1������aΪ����
        %if norm(w(i,:,k))>a*exp(-b*t(k))    %�ж��Ƿ�����ٷ�����2������a,bΪ����
            tk(i)=k;    %����ʱ��
            trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
        end
        %}
    end
    u(k+1)=Markov(u(k),P);     %uΪ�������Ĺ��  
    %
    %3����ɢ�����µĸ��¹���
         %%����ʱ���Ԥ��
    if k==1|u(k+1)~=u(k)    %������ڳ�ʼ״̬ʱ�Լ��������л�ʱȫ������Ԥ��
    for i=1:nodenum
        zsum=zeros(1,nodedim);
        for j=1:nodenum    %����zi(t)�е������
            zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
        end
        z=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
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
                zj=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
                zj(j,:,tkj(j))=-sita_t*c_t*zsum-tanh*D(j,j)*(g(x(j,:,tkj(j)))-g(s(tkj(j),:)))*Gama';
                ysum=ysum-L(i,j)*phi(y,z(i,:,tk(i)),zj(j,:,tkj(j)),x(i,:,tk(i)),x(j,:,tk(i)));
             end
             upwi=ysum+tanh*D(i,i)*phi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi���¹����в��Ⱥ���ߵ���
             %if upwi>a1*psi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %�������Ƹ��¹���       
             if upwi>a*exp(-b*(y+t(tk(i))))        %ָ�����Ƹ��¹���
                yi=y-dt;
                break;      %����whileѭ����
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
        if yy(i)==dt       %Ԥ������պõ�����һ��ʱ���ʱ�ͷ�������
           tk(i)=k;
           trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
           tri(i)=1;    
        end
    end
    for i=1:nodenum       %�ж��ھӽڵ��Ƿ��и���
        if tri(i)==1
           zsum=zeros(1,nodedim);
            for j=1:nodenum    %����zi(t)�е������
                zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
            end
            z=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
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
                    zj=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
                    zj(j,:,tkj(j))=-sita_t*c_t*zsum-tanh*D(j,j)*(g(x(j,:,tkj(j)))-g(s(tkj(j),:)))*Gama';
                    ysum=ysum-L(i,j)*phi(y,z(i,:,tk(i)),zj(j,:,tkj(j)),x(i,:,tk(i)),x(j,:,tk(i)));
                end
                upwi=ysum+tanh*D(i,i)*phi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi���¹����в��Ⱥ���ߵ���
                %if upwi>a1*psi(yy(i),z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %�������Ƹ��¹���       
                if upwi>a*exp(-b*(y+t(tk(i))))        %ָ�����Ƹ��¹���
                    yi=y-dt;
                    break;      %����whileѭ����
                end
                yy(i)=yy(i)+dt;
                if yy(i)>0.5
                    break;
                end
            end
        end
    end   
        %��ɢ�����¸��¹������
        %}
    s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %Ŀ�꺯���ĵ���
    e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %�ڵ���Ŀ������������еĽڵ��Լ���Ӧ״̬��һ�����
    L=LL(:,:,u(k+1));      %������������Ͻṹ 
    D=DD(:,:,u(k+1));      %���������ƿ��ƽڵ㼯 
      %1���޸��¹�����������
      %  tk=(k+1).*ones(nodenum,1);
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
ee=((ee1+ee2+ee3)).^(1/2);
%{
figure;
subplot(3,1,1);     %״̬1�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x1(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{1}, i=1,����,5');
 
subplot(3,1,2);     %״̬2�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x2(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{2}, i=1,����,5');

subplot(3,1,3);     %״̬3�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x3(i,:));
    hold on;
end
xlabel('t');
ylabel('x_{i}^{3}, i=1,����,5');
%}
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
legend('x_{i}^{1},the state 1 trajectory of nodes','x_{i}^{2},the state 2 trajectory of nodes','x_{i}^{3},the state 3 trajectory of nodes')
xlabel('t');
ylabel('x_{i}^{k}, i=1,����,5');


figure;    
               %״̬1�����ͼ
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

                 %״̬2�����ͼ
for i=2:nodenum
    plot(t,e2(i,:),'g');
    hold on;
end
                %״̬3�����ͼ
for i=2:nodenum
    plot(t,e3(i,:),'r');
    hold on;
end
legend('e_{i}^{1},the state 1 error of nodes','e_{i}^{2},the state 2 error of nodes','e_{i}^{3}��the state 3 error of nodes')
xlabel('t');
ylabel('error');

figure;     %������ŵ����ͼ��
plot(t,V);
xlabel('t');
ylabel('V(t)');


figure;        %�������ͼ
plot(t,ee);
xlabel('t');
ylabel('E(t)');

%{
figure;
subplot(2,1,1);     %״̬1�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x1(i,:));
    hold on;
end
ylabel('x_{i}^{1}, i=1,����,5');
subplot(2,1,2);     %״̬1�����ͼ
for i=1:nodenum
    plot(t,e1(i,:));
    hold on;
end
xlabel('t');
ylabel('{e_{i}^{1}}, i=1,����,5');

figure;
subplot(2,1,1);     %״̬2�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x2(i,:));
    hold on;
end
ylabel('x_{i}^{2}, i=1,����,5');
subplot(2,1,2);      %״̬2�����ͼ
for i=1:nodenum
    plot(t,e2(i,:));
    hold on;
end
xlabel('t');
ylabel('{e_{i}^{2}}, i=1,����,5');

figure;
subplot(2,1,1);     %״̬3�Ĺ켣ͼ
for i=1:nodenum
    plot(t,x3(i,:));
    hold on;
end
ylabel('x_{i}^{3}, i=1,����,5');
subplot(2,1,2);     %״̬3�����ͼ
for i=1:nodenum
    plot(t,e3(i,:));
    hold on;
end
xlabel('t');
ylabel('e_{i}^{3}, i=1,����,5');
%}
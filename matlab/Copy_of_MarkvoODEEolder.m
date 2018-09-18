function Copy_of_MarkvoODEEolder()
tic;
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
      0 1 0 0 -1;...7
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
        %if norm(w(i,:,k))>a1*norm(e(i,:,k))    %�ж��Ƿ�����ٷ�����1������aΪ����
        if norm(w(i,:,k))>a*exp(-b*t(k))    %�ж��Ƿ�����ٷ�����2������a,bΪ����
            tk(i)=k;    %����ʱ��
            trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
        end
        %}
        %
        %3����ɢ�����µĸ��¹���
            %%����ʱ���Ԥ��
        zsum=zeros(1,nodedim);
        for j=1:nodenum    %����zi(t)�е������
            zsum=zsum+L(i,j)*g(x(j,:,tk(i)))*Gama';
        end
        z=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
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
                zj=zeros(nodenum,nodedim,timedim);   %ϵͳ�г���f��ʣ�µ�������
                zj(j,:,tkj(j))=-sita_t*c_t*zsum-tanh*D(j,j)*(g(x(j,:,tkj(j)))-g(s(tkj(j),:)))*Gama';
                ysum=ysum-L(i,j)*phi(y,z(i,:,tk(i)),zj(j,:,tkj(j)),x(i,:,tk(i)),x(j,:,tk(i)));
             end
             upwi=ysum+tanh*D(i,i)*phi(y,z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:));        %pi���¹����в��Ⱥ���ߵ���
             %if upwi>a1*pci(y,z(i,:,tk(i)),0,x(i,:,tk(i)),s(tk(i),:))       %�������Ƹ��¹���       
             if upwi>a*exp(-b*(y+t(tk(i))))        %ָ�����Ƹ��¹���
                yi=y-dt;
                break;      %����whileѭ����
             end
             y=y+dt;
             if y>10*dt
                 break;
             end
        end
        for j=1:nodenum       %�ж��ھӽڵ��Ƿ��и���
            if j==i
                continue;
            end
            if t(k)==t(tk(j))
                judge=1;
            else
                judge=0;
            end
        end
        if yi==dt       %Ԥ������պõ�����һ��ʱ���ʱ�ͷ�������
           tk(i)=k;
           trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
        elseif judge==1     %�ھӽڵ㷢������ʱҲҪ����
               tk(i)=k;
               trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
        elseif k>1      
            if u(k)~=u(k-1)     %��������������л���Ҳ��������
               tk(i)=k;
               trinum(i)=trinum(i)+1;    %��¼�ڵ�i��������
            end 
        else
        end 
        %��ɢ�����¸��¹������
        %}
    end
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
%{
figure;     %������ŵ����ͼ��
plot(t,V);
xlabel('t');
ylabel('V(t)');
%}

figure;        %�������ͼ
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
    display(strcat('����ʱ�䣺',num2str(runtime),'��'));
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
    %f(x)���ſɱȾ���
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
    random=rand(1);  %����[0��1]�Ͼ��ȷֲ�������������ڿ̻�ת�Ƶ������
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

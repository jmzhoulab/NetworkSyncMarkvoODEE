function [trinum,X,E,tt,ss]=Fistpaper_parfor(L1,L2,L3,varargin)     % varargin �ǿɱ����������Ϊ���ǿ�ȡ�����ǿ�ȡ�ʱ�䳤�ȣ��Լ��������� 
    global nodenum nodedim t s;
    c=7;    tanh=3.5;  tn=2;  dt=0.001;     %Ĭ�����ǿ�ȡ�����ǿ�ȡ�ʱ�䳤�ȡ��Լ��������� 
    if numel(varargin)==0
    elseif numel(varargin)==1
            c=varargin{1}; 
    elseif numel(varargin)==2
            c=varargin{1}; tanh=varargin{2};  
    elseif numel(varargin)==3
            c=varargin{1}; tanh=varargin{2};    tn=varargin{3};
    elseif numel(varargin)==4
            c=varargin{1}; tanh=varargin{2};  tn=varargin{3};  dt=varargin{4};
    else
    disp('Parameter is wrong, please input again!');
    return;
    end;
    t=0:dt:tn;       %����ʱ������
    nodenum=length(L1);    nodedim=3;    timedim=length(t);
    P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %����������ת�ƾ���
    Gama=eye(nodedim);      %����Ͻṹ����   
    pinnum=round(nodenum/5);
    D1=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D1);  diagD(pin)=ones(length(pin),1); D1=diag(diagD)+D1-diag(diag(D1));
    D2=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D2);  diagD(pin)=ones(length(pin),1); D2=diag(diagD)+D2-diag(diag(D2));
    D3=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D3);  diagD(pin)=ones(length(pin),1); D3=diag(diagD)+D3-diag(diag(D3));
    u=[1];     L=L1;    D=D1;       %��������ʼ״̬     %���������ڳ�ʼ״̬�µ�����Ͼ���
    a1=0.3649;    a=0.02;    b=0.2;    p=0.8;      %���������ϵĸ���
    trinum=zeros(nodenum,1);    %��¼�����ڵ�ļ�������
    X=zeros(nodenum*nodedim,timedim);     X(:,1)=4*(rand(nodenum*nodedim,1)-0.5);     %����״̬���󲢽��г�ʼ��
    E=X;  s=zeros(nodedim,timedim);  s(:,1)=4*(rand(nodedim,1)-0.5);    S=kron(ones(nodenum,1),s(:,1));       %������Ķ����Ŀ�꺯���ĳ�ʹ��
    Xk=kron(ones(1,nodenum),X(:,1));     Sk=kron(ones(1,nodenum),S(:,1));       %��ʼ������ʱ�̽ڵ�״̬����
    Xksize = size(Xk);
    Sksize = size(Sk);
    HXk=H(Xk);    HSk=H(Sk);        %�����Ժ���h()�����µļ���ʱ�̽ڵ�״̬����
    LHXk=diagbrock(kron(L,Gama),HXk);    DHXSk=diagbrock(kron(D,Gama),HXk-HSk);    
    for k=1:timedim-1
        tic;
        sita_t=binornd(1,p);         %��Ŭ�������������������
        c_t=2*c*rand(1);     %������ǿ��Ϊ[0,2c]�ϵľ��ȷֲ�
        [FXk, HXt]=FH(X(:,k));
        X(:,k+1)=X(:,k)+(FXk-sita_t*c_t*LHXk-tanh*c_t*DHXSk)*dt;      %��������
        %X(:,k+1)=X(:,k)+(F(X(:,k))-sita_t*c_t*kron(L,Gama)*H(X(:,k))-tanh*c_t*kron(D,Gama)*(X(:,k)-S))*dt;
        clear FXk;
        R=sita_t*(kron(L,Gama)*HXt-LHXk)-tanh*(kron(D,Gama)*(HXt-S)-DHXSk);     %д��kron���������
        s(:,k+1)=s(:,k)+f(s(:,k))*dt;   %Ŀ�꺯���ĵ���
        S=kron(ones(nodenum,1),s(:,k+1));
        E(:,k+1)=X(:,k+1)-S;       %�ڵ���Ŀ������������еĽڵ��Լ���Ӧ״̬��һ�����
        Rt=reshape(R,nodenum,nodedim);
        Et=reshape(E(:,k+1),nodenum,nodedim);
        triggeif=sqrt(sum(Rt'.^2))>a1*sqrt(sum(Et'.^2));            %���������µļ�������1
        %triggeif=sqrt(sum(Rt'.^2))>a*exp(-b*t(k));              %���������µļ�������2
        triggeif=sqrt(sum(Rt'.^2))>-1;
        trinum=trinum+triggeif';
        if sum(triggeif)>0
            [XH,SH]=XSH(kron(ones(1,sum(triggeif)),X(:,k+1)),kron(ones(1,sum(triggeif)),S));
           	HXk(:,triggeif)=XH;       %���¶�Ӧ�����ڵ����Ϣ
            HSk(:,triggeif)=SH;        %����Ŀ��״̬����Ϣ
            LHXk =diagbrock(kron(L,Gama),XH);
            DHXSk=diagbrock(kron(D,Gama),XH-HSk);
        end
        u(k+1)=Markov(u(k),P);     %uΪ�������Ĺ��
        switch u(k+1)
            case 1,
                L=L1;   D=D1;
            case 2,
                L=L2;   D=D2;
            otherwise,
                L=L3;   D=D3;
        end
        toc
    end
    figplot(X,E);
    runtime=toc;
    display(strcat('����ʱ�䣺',num2str(runtime),'��'));
    ss=s;
    tt=t;
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
    %v=[v,random,min,max];
end
function D=diagbrock(A,B)   %��nodenum��nodedim��ȡ����A��B�ĳ˻���ĶԽǿ飬����ɾ���
    global nodenum nodedim;
    D=zeros(nodenum*nodedim,1);
    for k=1:nodenum
        D((k-1)*nodedim+1:k*nodedim,:)=A((k-1)*nodedim+1:k*nodedim,:)*B(:,k);
    end
end 
function [F,H]=FH(X)
    nodenum=100; nodedim=3;
    X=reshape(X, nodedim, nodenum);
    F=zeros(nodedim, nodenum);
    H=zeros(nodedim, nodenum);
    for j=1:nodenum
        F(:,j)=f(X(:,j));
        H(:,j)=g(X(:,j));
    end
    F=F(:);
    H=H(:);
end
function H=H(X)
    nodenum=100; nodedim=3;
    sizexk=size(X);
    H=zeros(sizexk(1),sizexk(2));
    for j=1:sizexk(2)
        Xj=reshape(X(:,j),nodedim, nodenum);
        Hj=zeros(nodedim,nodenum);
        for k=1:nodenum
            Hj(:,k)=g(Xj(:,k));
        end
        H(:,j)=Hj(:);
    end
end
function [XH,SH]=XSH(X,S)
    nodenum=100; nodedim=3;
    sizexk=size(X);
    XH=zeros(sizexk(1),sizexk(2));
    SH=zeros(sizexk(1),sizexk(2));
    for j=1:sizexk(2)
        Xj=reshape(X(:,j),nodedim, nodenum);
        Sj=reshape(S(:,j),nodedim, nodenum);
        XHj=zeros(nodedim,nodenum);
        SHj=zeros(nodedim,nodenum);
        for k=1:nodenum
            XHj(:,k)=g(Xj(:,k));
            SHj(:,k)=g(Xj(:,k));
        end
        XH(:,j)=XHj(:);
        SH(:,j)=SHj(:);
    end
end
function fx=f(x)
    p=9.78;
    q=14.97;
    m0=-1.31;
    m1=-0.75;
    fx = x;
    fx(1)= p*(-x(1)+x(2)-m1*x(1)+1/2*(m0-m1)*(abs(x(1)+1)-abs(x(1)-1)));   
    fx(2)= x(1)-x(2)+x(3);
    fx(3)= -q*x(2);
end
function gx=g(x)
     gx=zeros(1,3);
     gx(1)=2*x(1)+0.2*sin(x(1));
     gx(2)=2*x(2)+1;
     gx(3)=3*x(3)+0.5*cos(2*x(3));
end
function figplot(X,E)
    global nodenum nodedim t s;
    color='rgb';
    figure;     %�ڵ��״̬�켣ͼ
    for j=1:nodedim
       plot(t,X(j:nodedim:nodenum*nodedim,:),color(j));
       hold on;
       plot(t,s(j,:),color(j));
       hold on;
    end
    figure;     %�ڵ����켣ͼ
    for j=1:nodedim
       plot(t,sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)./nodenum),color(j));
       hold on;
    end
    ee=zeros(nodenum,length(t));
    for j=1:nodedim
       ee(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)./nodenum);
    end
    E=sum(ee)/nodedim;
    figure;
    plot(t,E);
end
            
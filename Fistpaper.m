function [trinum,X,E,tt,ss]=Fistpaper(L1,L2,L3,varargin)     % varargin ÊÇ¿É±ä²ÎÊý£¬ÒÀ´ÎÎªñîºÏÇ¿¶È¡¢¿ØÖÆÇ¿¶È¡¢Ê±¼ä³¤¶È£¬ÒÔ¼°µü´ú²½³¤ 
    tic;
    global nodenum nodedim t s;
    c=7;    tanh=3.5;  tn=2;  dt=0.001;     %Ä¬ÈÏñîºÏÇ¿¶È¡¢¿ØÖÆÇ¿¶È¡¢Ê±¼ä³¤¶È¡¢ÒÔ¼°µü´ú²½³¤ 
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
    t=0:dt:tn;       %µü´úÊ±¼äÐòÁÐ
    nodenum=length(L1);    nodedim=3;    timedim=length(t);
    P=[0.2 0.4 0.4;0.5 0.2 0.3;0.1 0.7 0.2];       %ÂíÊÏÁ´¸ÅÂÊ×ªÒÆ¾ØÕó
    Gama=eye(nodedim);      %ÄÚñîºÏ½á¹¹¾ØÕó   
    pinnum=round(nodenum/5);
    D1=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D1);  diagD(pin)=ones(length(pin),1); D1=diag(diagD)+D1-diag(diag(D1));
    D2=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D2);  diagD(pin)=ones(length(pin),1); D2=diag(diagD)+D2-diag(diag(D2));
    D3=zeros(nodenum,nodenum);    pin=randperm(nodenum,pinnum);    diagD=diag(D3);  diagD(pin)=ones(length(pin),1); D3=diag(diagD)+D3-diag(diag(D3));
    u=[1];     L=L1;    D=D1;       %ÂíÊÏÁ´³õÊ¼×´Ì¬     %ÂíÊÏÁ´´¦ÓÚ³õÊ¼×´Ì¬ÏÂµÄÍâñîºÏ¾ØÕó
    a1=0.3649;    a=0.02;    b=0.2;    p=0.8;      %Ëæ»ú·¢ÉúñîºÏµÄ¸ÅÂÊ
    trinum=zeros(nodenum,1);    %¼ÇÂ¼¸÷¸ö½ÚµãµÄ¼¤·¢´ÎÊý
    X=zeros(nodenum*nodedim,timedim);     X(:,1)=4*(rand(nodenum*nodedim,1)-0.5);     %¶¨Òå×´Ì¬¾ØÕó²¢½øÐÐ³õÊ¼»¯
    E=X;  s=zeros(nodedim,timedim);  s(:,1)=4*(rand(nodedim,1)-0.5);    S=kron(ones(nodenum,1),s(:,1));       %Îó²î¾ØÕóµÄ¶¨ÒåºÍÄ¿±êº¯ÊýµÄ³õÊ¹»¯
    Xk=kron(ones(1,nodenum),X(:,1));     Sk=kron(ones(1,nodenum),S(:,1));       %³õÊ¼»¯¼¤·¢Ê±¿Ì½Úµã×´Ì¬¾ØÕó
    HXk=H(Xk);    HSk=H(Sk);        %·ÇÏßÐÔº¯Êýh()×÷ÓÃÏÂµÄ¼¤·¢Ê±¿Ì½Úµã×´Ì¬¾ØÕó
    LHXk=diagbrock(kron(L,Gama),HXk);    DHXSk=diagbrock(kron(D,Gama),HXk-HSk);    
    for k=1:timedim-1
        sita_t=binornd(1,p);         %²®Å¬ÀûËæ»ú±äÁ¿ÃèÊöËæ»úñîºÏ
        c_t=2*c*rand(1);     %Ëæ»úñîºÏÇ¿¶ÈÎª[0,2c]ÉÏµÄ¾ùÔÈ·Ö²¼
        X(:,k+1)=X(:,k)+(F(X(:,k))-sita_t*c_t*LHXk-tanh*c_t*DHXSk)*dt;      %µü´ú·½³Ì
        %X(:,k+1)=X(:,k)+(F(X(:,k))-sita_t*c_t*kron(L,Gama)*H(X(:,k))-tanh*c_t*kron(D,Gama)*(X(:,k)-S))*dt;
        HXt=H(X(:,k));
        R=sita_t*(kron(L,Gama)*HXt-LHXk)-tanh*(kron(D,Gama)*(HXt-S)-DHXSk);     %Ð´³Ékron»ýºóµÄÓàÏî
        s(:,k+1)=s(:,k)+f(s(:,k))*dt;   %Ä¿±êº¯ÊýµÄµü´ú
        S=kron(ones(nodenum,1),s(:,k+1));
        E(:,k+1)=X(:,k+1)-S;       %½ÚµãÓëÄ¿±êµÄÎó²î£¬ÕâÀïËùÓÐµÄ½ÚµãÒÔ¼°¶ÔÓ¦×´Ì¬Îó²î¶¼Ò»Æð¼ÆËã
        Rt=reshape(R,nodenum,nodedim);
        Et=reshape(E(:,k+1),nodenum,nodedim);
        triggeif=sqrt(sum(Rt'.^2))>a1*sqrt(sum(Et'.^2));            %Á¬ÐøÇéÐÎÏÂµÄ¼¤·¢¹æÔò1
        %triggeif=sqrt(sum(Rt'.^2))>a*exp(-b*t(k));              %Á¬ÐøÇéÐÎÏÂµÄ¼¤·¢¹æÔò2
        triggeif=sqrt(sum(Rt'.^2))>-1;
        trinum=trinum+triggeif';
        if sum(triggeif)>0
            HXk(:,triggeif)=H(kron(ones(1,sum(triggeif)),X(:,k+1)));        %¸üÐÂ¶ÔÓ¦¼¤·¢½ÚµãµÄÐÅÏ¢
            HSk(:,triggeif)=H(kron(ones(1,sum(triggeif)),S));          %¸üÐÂÄ¿±ê×´Ì¬µÄÐÅÏ¢
            LHXk=diagbrock(kron(L,Gama),HXk);    
            DHXSk=diagbrock(kron(D,Gama),HXk-HSk);
        end
        u(k+1)=Markov(u(k),P);     %uÎªÂíÊÏÁ´µÄ¹ìµÀ
        switch u(k+1)
            case 1,
                L=L1;   D=D1;
            case 2,
                L=L2;   D=D2;
            otherwise,
                L=L3;   D=D3;
        end
    end
    figplot(X,E);
    runtime=toc;
    display(strcat('¼ÆËãÊ±¼ä£º',num2str(runtime),'Ãë'));
    ss=s;
    tt=t;
end

function v=Markov(u,p)
    %v=[0];
    n=length(p);
    s=1:n;
    pp=[p(u,:);s];
    random=rand(1);  %²úÉú[0£¬1]ÉÏ¾ùÔÈ·Ö²¼µÄËæ»úÊý£¬ÓÃÓÚ¿Ì»­×ªÒÆµÄËæ»úÐÔ
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
function D=diagbrock(A,B)   %ÒÔnodenumÐÐnodedimÁÐÈ¡¾ØÕóA£¬BµÄ³Ë»ýºóµÄ¶Ô½Ç¿é£¬²¢×é³É¾ØÕó
    global nodenum nodedim;
    D=zeros(nodenum*nodedim,1);
    for k=1:nodenum
        D((k-1)*nodedim+1:k*nodedim,:)=A((k-1)*nodedim+1:k*nodedim,:)*B(:,k);
    end
end 
function F=F(X)
global nodenum nodedim;
F=zeros(length(X),1);
for i=1:nodenum
    F((i-1)*nodedim+1:i*nodedim)=f(X((i-1)*nodedim+1:i*nodedim));
end
end
function H=H(X)
global nodenum nodedim;
sizexk=size(X);
H=zeros(sizexk(1),sizexk(2));
for j=1:sizexk(2)
    for i=1:nodenum
        H((i-1)*nodedim+1:i*nodedim,j)=g(X((i-1)*nodedim+1:i*nodedim,j));
    end
end
function gx=g(x)
 gx=zeros(1,3);
 gx(1)=2*x(1)+0.2*sin(x(1));
 gx(2)=2*x(2)+1;
 gx(3)=3*x(3)+0.5*cos(2*x(3));
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
function figplot(X,E)
    global nodenum nodedim t s;
    color='rgb';
    figure;     %½Úµã¸÷×´Ì¬¹ì¼£Í¼
    for j=1:nodedim
       plot(t,X(j:nodedim:nodenum*nodedim,:),color(j));
       hold on;
       plot(t,s(j,:),color(j));
       hold on;
    end
    figure;     %½ÚµãÎó²î¹ì¼£Í¼
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
            
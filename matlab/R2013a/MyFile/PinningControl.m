function [x,e,t,s] = PinningControl(L1)
    tic;
    global nodenum nodedim timedim t s x e;
    kapa=0.3
    delta=0.6
    tn=2;
    dt=0.001;     %Ĭ�����ǿ�ȡ�����ǿ�ȡ�ʱ�䳤�ȡ��Լ��������� 
    t=0:dt:tn;       %����ʱ������
    nodenum=length(L1);    nodedim=3;

    inititspan=2*(2*rand(nodenum,nodedim)-0.5);
    s=zeros(nodedim,timedim);  s(:,1)=4*(rand(nodedim,1)-0.5);
    timedim=length(t);      %��������
    x=zeros(nodenum,nodedim,timedim);       %����һ����ά�ľ������ڴ�����ݣ���һά�ȱ�ʾ�ڵ㣬�ڶ�ά�ȱ�ʾ�ڵ�״̬������ά�ȱ�ʾʱ��
    for i=1:timedim     %��ʱ��t��������x�ĵ���ά��
        x(:,:,i)=t(i);
    end
    x(:,:,1)=inititspan;      %�����ֵ
    Sigmma = zeros(nodenum,nodedim,timedim);
    for i=1:nodenum
        for j=1:nodenum
            Sigmma(i,j,1)=norm(x(i,:,1)-x(j,:,1));
        end
    end
    e=zeros(nodenum,nodedim,timedim);   %������
    q=zeros(nodenum,timedim);
    e(:,:,1)=x(:,:,1)-kron(ones(nodenum,1),s(1,:));       %��ʼ���
    for i=1:nodenum
        q(i,1)=kapa*norm(e(i,:,1));
    end
    for k=1:timedim-1       %ʱ�䷽���ϵ�ѭ��
        for i=1:nodenum
            cp=zeros(1,nodedim);
            for j=1:nodenum        %����ϵͳ�е������������ĺ� 
                q=Sigmma(i,j,k)*g(x(j,:,tk(i))); 
                cp=cp+q;
                Sigmma(i,j,k)=Sigmma(i,j,k)+alpha*norm(x(i,:,1)-x(j,:,1));
            end
            x(i,:,k+1)=x(i,:,k)+(f(x(i,:,k))-cp-normrnd(delta,1)*q(i,k)*(x(i,:,k)-s(k,:)))*dt;   %��������
            q(i,k)=q(i,k)+kapa*norm(e(i,:,k));
        end
        s(k+1,:)=s(k,:)+f(s(k,:))*dt;   %Ŀ�꺯���ĵ���
        e(:,:,k+1)=x(:,:,k+1)-kron(ones(nodenum,1),s(k+1,:));       %�ڵ���Ŀ������������еĽڵ��Լ���Ӧ״̬��һ�����

        if mod(k, 10) == 0
            runtime=toc;
            display(strcat('����ʱ�䣺',num2str(runtime),'��'));
        end
    end
    figplot()
    runtime=toc;
    display(strcat('�ܼ���ʱ�䣺',num2str(runtime),'��'));

end

function gx=g(x)
     gx=zeros(1,3);
     gx(1)=2*x(1)+0.2*sin(x(1));
     gx(2)=2*x(2)+1;
     gx(3)=3*x(3)+0.5*cos(2*x(3));
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

function figplot()
    global nodenum nodedim t s X E;
    color='rgb';
    figure;     %�ڵ��״̬�켣ͼ
    for j=1:nodedim
        plot(t,s(j,:),color(j));
        hold on;
    end
    for j=1:nodedim
       plot(t,X(j:nodedim:nodenum*nodedim,:),color(j));
       hold on;
    end
    se=zeros(nodedim,length(t));
    figure;     %�ڵ����켣
    for j=1:nodedim
       plot(t,sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)),color(j));
       hold on;
    end
    for j=1:nodedim
        se(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2));
    end
    talle=sum(se);
    figure;
    plot(t,talle);
end


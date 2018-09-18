function adaptfigplot(nodenum, nodedim, t, s, X, E, q, Sigmma)
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
%     figure;     %�ڵ����켣
%     for j=1:nodedim
%        plot(t,sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2)),color(j));
%        hold on;
%     end
    for j=1:nodedim
        se(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2));
    end
    talle=sum(se);
    figure;
    for i=1:nodenum
        plot(t,q(i,:),'r');
        hold on;
    end
    figure;
    for i=1:length(Sigmma(:,1))
        plot(t,Sigmma(i,:),'r');
        hold on;
    end
    figure;
    plot(t,talle);
end

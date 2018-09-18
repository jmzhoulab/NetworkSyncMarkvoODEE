function figplot(X,E,nodenum,nodedim,t,s)
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

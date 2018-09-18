tau=5.2;
c=2.1;
k1=(1+2.1764*tau)*c-25.404
k2=(1+3.71534*tau)*c-25.9818
k3=(1+3.7153*tau)*c-23.6702
figure; plot(tv,eewan2v,'g');axis([0 t(end) 0 0.08]);
figure; plot(tv,eewan1v,'b');axis([0 t(end) 0 0.08]);
figure; plot(tv,eewan3v,'r');axis([0 t(end) 0 0.08]);
nodeNurm=zeros(nodenum,timedim);
for k=1:length(t)
    for i=1:nodenum
        nodeNurm(i,k)=norm(ewan(i,:,k));    %第i个节点在时刻k的范数
    end
end
set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
      'DefaultAxesLineStyleOrder','-|--|:')
  figure;
  for i=1:nodenum
    plot(t,nodeNurm(i,:))
    hold all;
end
set(0,'DefaultAxesLineStyleOrder','remove')
set(0,'DefaultAxesColorOrder','remove')
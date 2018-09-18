ee=zeros(nodenum,length(t));
for j=1:nodedim
  ee(j,:)=sqrt(sum(E(j:nodedim:nodenum*nodedim,:).^2));
end
E=sum(ee);
figure;
plot(t,E);
hold on
E1=[e4;e2;e7;e8];
E2=[e3;e1;e5;e6];
color='rgbk';
figure;
for i=1:4
    plot(t,E1(i,:),color(i));
    hold on;
    plot(t,E2(i,:),strcat('--',color(i)));
    hold on;
end

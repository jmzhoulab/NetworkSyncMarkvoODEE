%¿Ó—≈∆’≈µ∑Ú∫Ø ˝
V0=zeros(1,timedim);
V9=zeros(1,timedim);
V10=zeros(1,timedim);
V24=zeros(1,timedim);
V25=zeros(1,timedim);
for k=1:timedim
    V0(k)=1/2*egen(:,k)'*egen(:,k);
    V9(k)=1/2*e2223(:,k)'*e2223(:,k);
    V10(k)=1/2*e1612(:,k)'*e1612(:,k);
    V24(k)=1/2*e2212(:,k)'*e2212(:,k);
    V25(k)=1/2*e3840(:,k)'*e3840(:,k);
end
% V=[V9;V10;V24;V25;V0;exp(-2*t)];
% color='rgbkcm';
% figure;
% for i=1:6
%     plot(t,V(i),color(i));
%     hold on;
% end
d=1;
tt=t(1:d:timedim);
CNVV=V9(1:d:timedim);
CEVV=V10(1:d:timedim);
DNVV=V24(1:d:timedim);
DEVV=V25(1:d:timedim);
noncontrolVV=V0(1:d:timedim);
figure
plot(tt,CNVV,'r-')     %Continuous monitoring under rule(9)
hold on
plot(tt,CEVV,'g-')     %Continuous monitoring under rule(10)
hold on
plot(tt,DNVV,'b-')    %Discrete monitoring under rule(40)
hold on
plot(tt,DEVV,'k-')     %Discrete monitoring under rule(41)
hold on
plot(tt,noncontrolVV,'b--')    %Continuous updating
hold on
plot(tt,exp(-2*tt),'k--')   %exp(-2t)
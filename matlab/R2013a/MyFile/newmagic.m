function G=newmagic(M)
a=size(M,1);
n=floor(a/2);
b=randperm(n);
c=a+1-b;
d=rand(1,n)>0.5;
e=~d;
e=[d;e];
f=[b;c];
F=[f(e);ceil(a/2);flipud(f(~e))];
M1=M(F,:);
G=M1(:,F);
end
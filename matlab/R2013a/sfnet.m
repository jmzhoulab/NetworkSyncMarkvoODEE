%%%%%%%%%%%%The BA scalefree model %%%%%%%
%%%%%%%%%%%%Write by Wu Ye %%%%%%%%%
n=500; % The total number of the net
a=zeros(n,n);
m=4; % The mean degree
% The initial random network
n0=5;
p0=0.8;
for i=1:n0
    for j=i+1:n0
        if rand(1,1)<p0
            a(i,j)=1;
            a(j,i)=1;
        end
    end
end
for i=1:n0
    deg(i,1)=sum(a(i,:));
end
%%%%%%%%%%%%%%%%%%%%
for i=n0:n-1
     b=zeros(m,1);
    [b]=scalefree(i,m,deg);
    for j=1:m
        a(b(j,1),i+1)=1;
        a(i+1,b(j,1))=1;
        deg(b(j,1),1)=deg(b(j,1),1)+1;
    end
    deg(i+1,1)=m;
end

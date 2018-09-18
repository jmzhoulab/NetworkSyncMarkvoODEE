function fv=V(S,D)
D1=zeros(1,length(D));D2=zeros(1,length(D));
u=zeros(1,length(S));v=zeros(1,length(S));
for i=1:length(S)
    u(i)=S{1,i}(1,1);
    v(i)=S{1,i}(1,2);
end
for i=1:length(D)
    D1(i)=D{1,i}(1,1);
    D2(i)=D{1,i}(1,2);
end
M=sparse(D1,ones(1,length(D)),D2);
index=1:length(S);re=0;
for k=1:length(index)-1
    p=combntns(index,k);
    for i=1:length(p)
        pi=p(i,:);
        p_i=setdiff(index,pi);
        foot=0;
        for j=1:length(pi)
            foot=foot+pi(j)*10^(length(pi)-j);
        end
        re=re+prod(u(pi))*prod(v(p_i))*M(foot,1);
    end
    if k==length(index)-1
       re=re+prod(u(index))*M(end,1);
    end
end
fv=re;
end
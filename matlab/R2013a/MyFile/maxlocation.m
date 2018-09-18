function w=maxlocation(A)
m=length(A);
B=zeros(1,m);
for j=1:m
    B(j)=sum(A(:,j))-A(j,j);
end
b=B(1);n=1;
for i=1:m
    for j=i+1:m
        if B(j)>=B(i)&&B(j)>=b
            n=j;
            b=B(j);
        else
            q=1;
        end
    end
end
w=n;
end
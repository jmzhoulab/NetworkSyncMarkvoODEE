function onetri=onetri(a,b)
A=[];
for i=1:a
    for j=1:b
        if i==j
            A(i,j)=1;
        else if abs(i-j)==1
            A(i,j)=1;
        end
    end
end
end
onetri=A;
function B=diagbrock(A,a,b)
    len=length(A);
    B=zeros(len,b);
    for k=1:len/a
        B((k-1)*a+1:k*a,:)=A((k-1)*a+1:k*a,(k-1)*b+1:k*b);
    end
end %取矩阵的对角块的元素，对角块的形式为a行b列
function B=diagbrock(A,a,b)
    len=length(A);
    B=zeros(len,b);
    for k=1:len/a
        B((k-1)*a+1:k*a,:)=A((k-1)*a+1:k*a,(k-1)*b+1:k*b);
    end
end %ȡ����ĶԽǿ��Ԫ�أ��Խǿ����ʽΪa��b��
A=[5 2 1;-1 4 2;2 -3 10];
b=[-12 20 3]';
 D=diag(diag(A));
 L=-tril(A,-1);
 U=-triu(A,1);
 B=inv(D-L)*U;  %迭代矩阵
 f=inv(D-L)*b;
 x=[0 0 0]';
 k=0;
 r=1;
 w=x;
 while r>10^(-4) %迭代停止条件，即误差小于10^（-4）
xx=x;
x=B*x+f;%迭代公式
k=k+1;%记录迭代次数
w=[w,x];%将每次迭代的结果以列向量存储在矩阵w中
r=norm(x-xx,inf);%求相邻两次迭代的无穷范数
 end
 k,w=vpa(w,8)
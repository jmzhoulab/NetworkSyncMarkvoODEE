A=[-6 1 1;1 -6 1;1 1 -6];
b=[-1130 -3200 -4200];
 D=diag(diag(A));
 B=eye(3)-inv(D)*A;%迭代矩阵
 f=inv(D)*b';
 q=norm(B,inf);
 x=[0 0 0]';
 k=0;
 r=1;
 X=x;
 while r>0.01 %迭代停止条件，即误差小于10^（-4）
xx=x;
x=B*x+f;%迭代公式
k=k+1;%记录迭代次数
X=[X,x];%将每次迭代的结果以列向量存储在矩阵X中
r=q/(1-q)*norm(x-xx,inf);%求相邻两次迭代的无穷范数
 end
 k,X=vpa(X(:,k+1),6)
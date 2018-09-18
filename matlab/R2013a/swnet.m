%%%%%%%%%The parameter of the small-world model%%%%%%%
%%%%%%%%% Write by Wu Ye %%%%%%%%%%%%%%% 
n=100;   % The total number of the network
m=2;     % The mean degree is 2*m
p=0.5;   % The remove probability
a=zeros(n,n); % The connection matrix
%%%%%%% The initial regular network
for i=1:n
    for j=1:m
        i0=i+j;
        i1=i-j;
        if i0>n
            i0=i0-n;
        end
        a(i,i0)=1;
        a(i0,i)=1;
    end
end
%%%%%%%%%%%%%%%%%%
a_record=zeros(n,n); 
for i=1:n
    for j=i+1:n
        if a(i,j)==1&&a_record(i,j)==0
            if rand(1,1)<p
                a(i,j)=0;       % if you want to use NW model,please remove the sentence
                a(j,i)=0;       % also this one           
                pd=1;
                while pd>0
                    i0=fix(rand(1,1)*n)+1;
                    if (i0-i)*(i0-j)==0
                        pd=1;
                    else
                        pd=0;
                    end
             end
                a(i,i0)=1;
                a(i0,i)=1; 
                a_record(i,i0)=1;
                a_record(i0,i)=1;   % mark a(i,i0) is a new edge which did not need to remove again
            end
        end
    end
end
%%%%%%%%%%% The view of the network
re=zeros(n*m,2);
k=0;
for i=1:n
    for j=1:n
        if a(i,j)==1
            k=k+1;
            re(k,1)=i;
            re(k,2)=j;
        end
    end
end
plot(re(:,1),re(:,2),'.')
%%%%%%%%%%%%%%%%%%%

    
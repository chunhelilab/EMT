function dx=EMTode_plus(t,x,flag,s,a,aa,b,k,n,mat)
%This is an implementation of solving the stable states of the ODEs

%%x:current location
%%a, aa, b, sa, sb, k, n: the parameters
%%mat:the matrix representing the relationship between genes


d=0.1;
num=16;
dx=zeros(2*num,1);


for i=1:num
    for j=1:num 
        if mat(j,i)==1
            dx(i)=dx(i)+a*x(j)^n/(s(j,i)^n+x(j)^n);
            if j==i
                dx(i+num)=dx(i+num)+2*(aa*n*x(j)^(n-1)*s(j,i)^n/(s(j,i)^n+x(j)^n)^2)*x(i+num);
                dx(i)=dx(i)+aa*x(j)^n/(s(j,i)^n+x(j)^n);
            end
        elseif mat(j,i)==-1
            if j>10 && i<11
                dx(i)=dx(i)+b*s(j,i)^n/(s(j,i)^n+x(j)^n);
            else    
                dx(i)=dx(i)+b*s(j,i)^n/(s(j,i)^n+x(j)^n);
                if j==i
                    dx(i+num)=dx(i+num)-2*(b*n*x(j)^(n-1)*s(j,i)^n/(s(j,i)^n+x(j)^n)^2)*x(i+num);
                end
            end
        end
    end
    dx(i)=dx(i)-k*x(i);
    dx(i+num)=dx(i+num)+2*d-2*k*x(i+num); 
end











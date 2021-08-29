function s = EMT_action(x,mat,dd,para) 
%This is the implementation of calculating the path action

%%x:transition path
%%mat:the matrix representing the relationship between genes
%%dd:diffusion coefficient
%%para: parameter n, a, aa, b, sa, sb, k

s=0;
ss=0;
lave=0;
num=16;
M=20;
diag=zeros(num,num);

for j=1:num
    for n=1:num
        if n==j
            diag(n,j)=1;
        else
            diag(n,j)=0;
        end
    end
end

for i=1:M-1
    lave=lave+dl(diag,x(i,:),x(i+1,:));
end

lave=lave/(M-1);
[vi,f]=veff(x(1,:),mat,dd,para);
[vf,f]=veff(x(M,:),mat,dd,para);
e_eff=0.5;

for i=1:M-1
    a=dshj(x(i,:),x(i+1,:),lave,e_eff,mat,dd,para);
    s=s+a;    %Total transition action
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = dl(inv, xi, xf)
%\delta l_{n,n+1} term in Equation (12)

l=0;
num=16;
for i=1:num
    for j=1:num
        l=l+inv(i,j)*(xf(i)-xi(i))*(xf(j)-xi(j));    
    end
end
l=sqrt(l);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds= dshj(xi, xf, lave,e_eff,mat,dd,para)
%Transition action of a line segment

lambda=10e-8;
num=16;
invd=zeros(num,num);
xc=zeros(num,1);
ds=0; fdl=0;
dx=zeros(num,1); fd=zeros(num,1);
mu=100;
diag=zeros(num,num);
for j=1:num
    for n=1:num
        if n==j
            diag(n,j)=1;
            invd(n,j)=1/dd;
        end
    end
end
for i=1:num
    xc(i)=(xf(i)+xi(i))/2;
end
[vi,f]=veff(xc,mat,dd,para);
for i=1:num
    dx(i)=xf(i)-xi(i);
end
for i=1:num
    for j=1:num
        fdl=fdl+invd(i,j)*f(j)*dx(i);
    end
end
ds=dl(invd,xi,xf)*sqrt(max(0,vi+e_eff))+lambda*(dl(diag,xi,xf)-lave)^2-0.5*fdl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vf,f]=veff(x,mat,dd,para)
%V(n) term in Equation (12)

num=16;
f=zeros(num,1);
df=zeros(num,num);
b1=0; b2=0;

f=zeros(num,1);df=f;
n=para(1);a=para(2);aa=para(3);b=para(4);sa=para(5);sb=para(6);k=para(7);

for i=1:num
    for j=1:num
        if j>10 && i<11
            s(j,i)=sb;
        else
            s(j,i)=sa;
        end
    end
end

for i=1:num    
    for j=1:num
        if mat(j,i)==1           
            f(i)=f(i)+a*x(j)^n/(s(j,i)^n+x(j)^n);
            if j==i
                df(i)=df(i)+(aa*n*x(j)^(n-1)*s(j,i)^n/(s(j,i)^n+x(j)^n)^2);
                f(i)=f(i)+aa*x(j)^n/(s(j,i)^n+x(j)^n);
            end
        elseif mat(j,i)==-1
            f(i)=f(i)+b*s(j,i)^n/(s(j,i)^n+x(j)^n);
            if j==i
                df(i)=df(i)-(b*n*x(j)^(n-1)*s(j,i)^n/(s(j,i)^n+x(j)^n)^2);
            end
        end
    end
    

    f(i)=f(i)-k*x(i);
    df(i)=df(i)-k;
end

invd=zeros(num,num);
diff=zeros(num,num);
for j=1:num
    for n=1:num
        if n==j
            diff(n,j)=dd;
            invd(n,j)=1/dd;
        end
    end
end

for i=1:num
    for j=1:num
        b1=b1+invd(i,j)*f(i)*f(j);
    end
end
for i=1:num
    b2=b2+df(i);
end
vf=0.25*b1+0.5*b2;

end
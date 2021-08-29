%This is an implementation of of simulating EMT, 
%including establishing energy landscape and calculating minimum action path

clear;

%% parameters
times1=1000;                     %Number of iterations 1
times=1000;                      %Number of iterations 2
mat=importdata('matrix.txt');    %Relationship between genes
num=16;                          %Number of the genes
k=1;                             %Degradation rate
n=4;                             %Hill coefficient
mnos=10;                         %Maximum number of solution
nos=1;                           %Number of solution(at least one)
w=zeros(1,mnos);                 %Weight of every solution
xx=zeros(2*num,1); 
sol=zeros(2*num,mnos);           %The solution of the ODE
iniodsol=zeros(2*num,mnos);      %Initial condition of different solutions
ini=zeros(2*num,1);              %Initial value
dd=1;                            %Diffusion coefficient
N=300;                           %Initial point selection times
s=zeros(num,num);                %Threshold matrix
a=2;aa=7.9;b=4;sa=5.8;sb=0.21;   %Tunable parameters

para=zeros(1,7);
para(1)=n;para(2)=a;para(3)=aa;para(4)=b;para(5)=sa;para(6)=sb;para(7)=k;
for i=1:num
    for j=1:num
        if j>10 && i<11
            s(j,i)=sb;
        else
            s(j,i)=sa;
        end
    end
end

%% Calculate the stable states

for i=1:num
    ini(i)=25*rand;    % Take a random initial condition
end
nn=1;


while 1     % Calculate the first solution             
    [T,Y]=ode15s('EMTode_plus',[0 0.05:nn*300],ini,[],s,a,aa,b,k,n,mat);
    xx=Y(length(Y),:);
    xx1=Y(length(Y)-1,:);
    xxx=0;
    dis=abs(xx-xx1);
    for jj=1:num                                   
        if dis(jj)>xx(jj)/1000 && dis(jj)>0.001    % Make sure the solution is stable
            xxx=1;
            break;
        end
    end
    if xxx==0
        break;
    else
        nn=nn+1;
    end
end
sol(:,1)=xx';    % The first solution

for i=1:N    % And then we put other solution into it
    for ii=1:num
        ini(ii)=25*rand;      % Take a random initial condition
    end

    nn=1;
    while 1                    
        [T,Y]=ode15s('EMTode_plus',[0 0.05:nn*300],ini,[],s,a,aa,b,k,n,mat);
        xx=Y(length(Y),:);
        xx1=Y(length(Y)-1,:);
        xxx=0;
        dis=abs(xx-xx1);
        for jj=1:num
            if dis(jj)>xx(jj)/1000 && dis(jj)>0.001    % Make sure the solution is stable
                xxx=1;
                break;
            end
        end
        
        if xxx==0
            break;
        else
            nn=nn+1;
        end
    end

    j=0;
    for kk=1:nos    % Put the new solution into it
        xxx=0;
        dist=abs(xx'-sol(:,kk));
        for jj=1:num
            if dist(jj)>xx(jj)/100 && dist(jj)>0.001
                xxx=1;
                break;
            end
        end

        if xxx==0               % If the solution is in the sol
            w(kk)=w(kk)+1;      % Count the number of every solution
            j=1;
        end
    end
    if j==0
        iniodsol(:,nos)=xx';    % Save the initial condition of different solutions
        nos=nos+1;
        sol(:,nos)=xx';
        w(nos)=w(nos)+1;
    end
end
ww=w/N;
w=ww(1:nos);    % The weight of solutions
xt=sol(1:num,1:nos);    % The solutions of the ODE
sigma=sol((num+1):2*num,1:nos);    % The standard deviation of solutions

%% Sort these solutions

para=zeros(2*num+1,nos);  
para(1:num,:)=xt;
para(num+1:2*num,:)=sigma;
para(end,:)=w;
para=para';
para=sortrows(para,num);
para=para';
xt=para(1:num,:);
sigma=para(num+1:2*num,:);
w=para(2*num+1,:);

%% Calculate the minimum action path (MAP)

M=20; %Number of segments of MAP
SS=xt;
actions1=zeros(5,times1);    %EMT path
actions2=zeros(5,times1);    %MET path
for j=1:5
    fprintf('Optimizing the %d segment of the path (total 5 segments)\n',j)
    fprintf('Optimizing...\n')
    for jj=1:2    
        path=zeros(num,M);  
        if jj==1
            iiia=j;iiib=j+1;    %EMT path
        else
            iiia=j+1;iiib=j;    %MET path
        end          
        for i=1:num
            path(i,:)=linspace(SS(i,iiia),SS(i,iiib),M);
        end

        a=SS(:,iiia);
        b=SS(:,iiib);
        lb=0.5*min(a,b);
        ub=1.5*max(a,b);

        % Simulated annealing
        t=10^(-5);
        itgl=0;itglr=0;
        rho=0.05;
        point=1;
        psr=path;
        itgl=EMT_action(path',mat,dd,para);
        amplitude=0.05;
        for i=1:times1
            a_num=0;
            if rho>0.4 
                amplitude=amplitude*1.2;
            end
            if rho<0.02
                amplitude=amplitude*0.8;
            end
            for kk=1:times
                psr=path;
                pj=2+fix((num-2)*rand);
                for mk=1:10000
                    for nk=2:(M-1)
                        psr(pj,nk)=path(pj,nk)+(max(path(pj,:))-min(path(pj,:)))*amplitude*(1-2*rand);
                    end
                    xxx=0;
                    for mi=2:(M-1)
                        if psr(pj,mi)>ub(pj) || psr(pj,mi)<lb(pj)
                            xxx=1;
                        end
                    end
                    if xxx==0
                        break;
                    end
                end
                itglr=EMT_action(psr',mat,dd,para);
                if itglr<=itgl
                    path=psr;
                    a_num=a_num+1;
                    itgl=itglr;
                else
                    l=exp((itgl-itglr)/t);

                    if rand<=l
                        path=psr;
                        itgl=itglr;
                        a_num=a_num+1;
                    end
                end    
            end
            rho=a_num/times;
            t=t*0.95;
            if jj==1
                actions1(j,i)=itgl;
            else
                actions2(j,i)=itgl;
            end
        end
        % End of simulated annealing
        
        if j==1    % Record paths
            if jj==1
                p12=path;
            else
                p21=path;
            end
        elseif j==2
            if jj==1
                p23=path;
            else
                p32=path;
            end
        elseif j==3            
            if jj==1
                p34=path;
            else
                p43=path;
            end
        elseif j==4
            if jj==1
                p45=path;
            else
                p54=path;
            end
        elseif j==5
            if jj==1
                p56=path;
            else
                p65=path;
            end
        end       
    end
end

%% Print the landscape

% the coordinates of the projected landscape
s1=2;    % ZEB2
s2=16;    % CDH1
% The sequence number is consistent with the gene sequence in Figure S8
% 1 for TGF-beta
% 2 for ZEB1
% 3 for ZEB2
% 4 for SNAI1
% 5 for SNAI2
% 6 for TWIST1
% 7 for FOXC2
% 8 for GSC
% 9 for TCF3
% 10 for VIM
% 11 for miR-145
% 12 for miR-141
% 13 for miR-200
% 14 for miR-34a
% 15 for Ovol2
% 16 for CDH1

gx=50;gy=50;
tt=500;
tx=gx/tt;
ty=gy/tt;
P=zeros(gx/tx+1,gy/ty+1);
fprintf('Drawing the landscape...\n')

for i=1:nos
    x=0:tx:gx;
    y=0:ty:gy;
    [x,y]=meshgrid(x,y);
    pi_sigma=1;
    for j=1:num
        pi_sigma=pi_sigma*sigma(j,i);
    end
    
    z=exp(-(x-xt(s1,i)).^2/(2*sigma(s1,i))-(y-xt(s2,i)).^2/(2*sigma(s2,i)))/((2*pi)^(num/2)*sqrt(pi_sigma));
    P=P+w(i)*z;

end

P=P./sum(sum(P));    
U=-log(P);           % Potential energy
height=200;          % The height of landscape
l=length(U);
for i=1:1:l      
    for j=1:1:l
        if U(i,j)>height
            U(i,j)=height;
        end
    end
end     

%% Print the surface

figure;
set(gcf,'outerposition', [100 80 750 600]);
kx=0:tx:gx;
ky=0:ty:gy;
surf(kx,ky,U);hold on
shading interp;
set(gca,'FontSize',22);
xlabel('ZEB1'),ylabel(['CDH1']),zlabel('U'),title('');

%% Generate mesh
x1=0;x2=30;y1=0;y2=30;
gridwidth=1;
[X,Y] = meshgrid(linspace(x1, x2, fix((x2-x1)/gridwidth)),linspace(y1, y2, fix((y2-y1)/gridwidth)));
for i=1:fix((y2-y1)/gridwidth)
    ppt=griddata(kx,ky,U,EMT_refine(X(i,:)),EMT_refine(Y(i,:)));
    plot3(EMT_refine(X(i,:)),EMT_refine(Y(i,:)),ppt,'Color',[.4 .4 .4],'LineWidth',.01);hold on;
end
for i=1:fix((x2-x1)/gridwidth)
    ppt=griddata(kx,ky,U,EMT_refine(X(:,i)'),EMT_refine(Y(:,i)'));
    plot3(EMT_refine(X(:,i)'),EMT_refine(Y(:,i)'),ppt,'Color',[.4 .4 .4],'LineWidth',.01);hold on;
end

%% Print the paths

hold on;
pp21=griddata(kx,ky,U,EMT_refine(p21(s1,:)),EMT_refine(p21(s2,:)));
plot3(EMT_refine(p21(s1,:)),EMT_refine(p21(s2,:)),pp21,'g','LineWidth',2);hold on;
pp32=griddata(kx,ky,U,EMT_refine(p32(s1,:)),EMT_refine(p32(s2,:)));
plot3(EMT_refine(p32(s1,:)),EMT_refine(p32(s2,:)),pp32,'g','LineWidth',2);hold on;
pp43=griddata(kx,ky,U,EMT_refine(p43(s1,:)),EMT_refine(p43(s2,:)));
plot3(EMT_refine(p43(s1,:)),EMT_refine(p43(s2,:)),pp43,'g','LineWidth',2);hold on;
pp54=griddata(kx,ky,U,EMT_refine(p54(s1,:)),EMT_refine(p54(s2,:)));
plot3(EMT_refine(p54(s1,:)),EMT_refine(p54(s2,:)),pp54,'g','LineWidth',2);hold on;
pp65=griddata(kx,ky,U,EMT_refine(p65(s1,:)),EMT_refine(p65(s2,:)));
plot3(EMT_refine(p65(s1,:)),EMT_refine(p65(s2,:)),pp65,'g','LineWidth',2);hold on;

pp12=griddata(kx,ky,U,EMT_refine(p12(s1,:)),EMT_refine(p12(s2,:)));
plot3(EMT_refine(p12(s1,:)),EMT_refine(p12(s2,:)),pp12,'r','LineWidth',2);hold on;
pp23=griddata(kx,ky,U,EMT_refine(p23(s1,:)),EMT_refine(p23(s2,:)));
plot3(EMT_refine(p23(s1,:)),EMT_refine(p23(s2,:)),pp23,'r','LineWidth',2);hold on;
pp34=griddata(kx,ky,U,EMT_refine(p34(s1,:)),EMT_refine(p34(s2,:)));
plot3(EMT_refine(p34(s1,:)),EMT_refine(p34(s2,:)),pp34,'r','LineWidth',2);hold on;
pp45=griddata(kx,ky,U,EMT_refine(p45(s1,:)),EMT_refine(p45(s2,:)));
plot3(EMT_refine(p45(s1,:)),EMT_refine(p45(s2,:)),pp45,'r','LineWidth',2);hold on;
pp56=griddata(kx,ky,U,EMT_refine(p56(s1,:)),EMT_refine(p56(s2,:)));
plot3(EMT_refine(p56(s1,:)),EMT_refine(p56(s2,:)),pp56,'r','LineWidth',2);hold on;

xlim([0 30])
ylim([10,30])
zlim([0,height])
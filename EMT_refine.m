function path = EMT_refine(pp)
%This is the implementation of refining the transition path
%It does not affect the transition path, but only beautifies it

%%pp:transition path

[x,y]=size(pp);
path=zeros(x,10*y-9);
for i=1:y
    path(:,10*i-9)=pp(:,i);
end
for i=1:(y-1)
    for j=1:x
        path(j,(10*i-9):(10*i+1))=linspace(pp(j,i),pp(j,i+1),11);
    end
end

end


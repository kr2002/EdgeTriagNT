% This function ordered piecewise constant graphons 
% according to the value of p_{ii}.
function [CC GG]=SortGraphon(C,G)

Nc=length(C);
for k=1:Nc
    if C(k)<1e-7
        for j=1:Nc
            G(k,j)=0.0; 
            G(j,k)=0.0;
        end
    end
end
GG=zeros(Nc,Nc);
CC=zeros(1,Nc);
[DD Ind]=sort(diag(G),'descend');
for j1=1:Nc
    for j2=1:Nc
        GG(j1,j2)=G(Ind(j1),Ind(j2));
    end
    CC(j1)=C(Ind(j1));
end

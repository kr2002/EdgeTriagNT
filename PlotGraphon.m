% This function visualizes a given graphon
function h=PlotGraphon(C,G);
Nc=length(C);
GG=zeros(Nc,Nc);
CC=zeros(1,Nc);
[DD Ind]=sort(diag(G),'descend');
for j1=1:Nc
    for j2=1:Nc
        GG(j1,j2)=G(Ind(j1),Ind(j2));
    end
    CC(j1)=C(Ind(j1));
end

X=zeros(1,Nc+1);
for jj=1:Nc
    X(jj+1)=X(jj)+CC(jj);
end
Y=X;
Z=zeros(Nc+1,Nc+1);
Z(1:Nc,1:Nc)=GG(:,:);
figure;
h=pcolor(X,Y,Z); colorbar;
caxis([0 1]);
axis square; axis tight;
axis off;
drawnow;

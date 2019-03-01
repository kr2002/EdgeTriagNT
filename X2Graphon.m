% Return graphon from the optimization variable X 
function [C G s]=X2Graphon(X,Nc)
Ng=Nc*(Nc+1)/2;
Ncg=Nc+Ng; % number of unknowns
C=X(1:Nc);
G=zeros(Nc,Nc);
for k1=1:Nc
	for k2=1:k1
		k12=Nc+(k1-1)*k1/2+k2;
		G(k1,k2)=X(k12);
	end
end
G=G+G'-diag(diag(G));
if nargout > 2
	s=gs(G,C);
end

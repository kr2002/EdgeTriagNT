%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function computes the rate function, 
% i.e. the negative entropy density for a given x 
% for the edge-triangle model
function [f g]=EdgeTriagNT_Obj(x,Nc)
C=x(1:Nc); % the first Nc elements of x contains c values
G=zeros(Nc,Nc);
for k1=1:Nc
    for k2=1:k1
        kg=Nc+(k1-1)*k1/2+k2;
        G(k1,k2)=x(kg);
    end
end
G=G+G'-diag(diag(G));
f=-gs(G,C);

if nargout > 1 % calculate gradient
	Ncg=length(x);
	g=zeros(Ncg,1);
	for k=1:Nc
		s=0.0;
		for j=1:Nc
			s=s+gS0(G(k,j))*C(j);
		end
		g(k)=-2*s;
	end
	for k1=1:Nc
		for k2=1:k1
			kg=Nc+(k1-1)*k1/2+k2;
			s=-gS0d(G(k1,k2))*C(k1)*C(k2);
			if k1==k2
				g(kg)=s;
			else
				g(kg)=2*s;
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function computes the e=e0 constraints
function [c ceq gc gceq] = ETNT_Con(x,e0,Nc)

[C G]=X2Graphon(x,Nc);
c=[];
ceq=[gE(G,C)-e0];

if nargout > 2 % calculate gradient
	gc=[];
	Ncg=length(x);
	gceq=zeros(Ncg,1);
	for p=1:Nc
		s1=0.0;
		for j=1:Nc
			s1=s1+C(j)*G(j,p);
		end
		gceq(p)=2*s1;
	end
	for k1=1:Nc
		for k2=1:k1
			kg=Nc+(k1-1)*k1/2+k2;
			s1=C(k1)*C(k2);
			if k1==k2
				gceq(kg)=s1;
			else
				gceq(kg)=2*s1;
			end
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
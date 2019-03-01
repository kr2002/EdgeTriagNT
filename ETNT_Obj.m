%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function computes the triangle density

function [f g]=ETNT_Obj(x,t0,Nc)

[C G]=X2Graphon(x,Nc);
t=gTTriag(G,C,Nc);
f=(t-t0)^2;

if nargout > 1 % calculate gradient
	Ncg=length(x);
	g=zeros(Ncg,1);
	for p=1:Nc
		s=0.0;
		for j=1:Nc
			for k=1:Nc
				s=s+C(j)*C(k)*G(j,k)*G(p,j)*G(k,p);
			end
		end
		g(p)=3*s;
	end
	for k1=1:Nc
		for k2=1:k1
			kg=Nc+(k1-1)*k1/2+k2;
			s=0.0;
			for k=1:Nc
				s=s+G(k2,k)*G(k,k1)*C(k1)*C(k2)*C(k);
			end
			if k1==k2
				g(kg)=3*s;
			else
				g(kg)=6*s;
			end
		end
    end
    g=2*(t-t0)*g;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
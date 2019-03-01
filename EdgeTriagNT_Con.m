%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function computes the (e, t) constraints for 
% the minimization problem of edge-triangle model
function [c ceq gc gceq] = EdgeTriagNT_Con(x,e0,t0,Nc)
C=x(1:Nc);
G=zeros(Nc,Nc);
for k1=1:Nc
    for k2=1:k1
        kg=Nc+(k1-1)*k1/2+k2;
        G(k1,k2)=x(kg);
    end
end
G=G+G'-diag(diag(G));
c=[];
ceq=[gE(G,C)-e0 gTTriag(G,C,Nc)-t0];

if nargout > 2 % calculate gradient
	gc=[];
	Ncg=length(x);
	gceq=zeros(Ncg,2);
	for p=1:Nc
		s1=0.0; s2=0.0;
		for j=1:Nc
			s1=s1+C(j)*G(j,p);
			for k=1:Nc
				s2=s2+C(j)*C(k)*G(j,k)*G(p,j)*G(k,p);
			end
		end
		gceq(p,1)=2*s1;
		gceq(p,2)=3*s2;
	end
	for k1=1:Nc
		for k2=1:k1
			kg=Nc+(k1-1)*k1/2+k2;
			s1=C(k1)*C(k2);	s2=0.0;
			for k=1:Nc
				s2=s2+G(k2,k)*G(k,k1)*C(k);
            end
            s2=s2*s1;
			if k1==k2
				gceq(kg,1)=s1;
				gceq(kg,2)=3*s2;
			else
				gceq(kg,1)=2*s1;
				gceq(kg,2)=6*s2;
			end
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

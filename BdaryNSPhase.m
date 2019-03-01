%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is used to find boundary of the symmetric N-podal phase
% in the edge-triangle model.
% The analytical formula for the minimizer is used.
%
% Last updated: 09/01/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X Y Id]=BdaryNSPhase(Nc)
Ec=[0.01 1/2 2/3 3/4 4/5 5/6];
E0=Ec(Nc-1):0.001:Ec(Nc+1);
NE=length(E0);
NT=800;

Id=zeros(NE,NT+1);
X=Id; Y=Id;
for ke=1:NE   
    E00=E0(ke);
    tmin=LBdaryNSEdgeTriag(Nc,E00);
    tmax=E00^3;
    dt=(tmax-tmin)/NT;
    T0=tmin:dt:tmax;    
    NTT=length(T0);
    if NTT~=NT+1
	disp('Inconsistency in NY');
	break;
    end
    for kt=1:NTT          

        T00=T0(kt);

        X(ke,kt)=E00;
        Y(ke,kt)=T00;
        ff=((-T00+E00^3)/(Nc-1))^(1/3);
        a=E00-(Nc-1)*ff; %
        b=E00+ff;
	a=max(eps,a); a=min(a,1-eps);
	b=max(eps,b); b=min(b,1-eps);

        %a11=4*(I0(a)-I0(b)-b*dI0(b)+a*dI0(a));
        %a12=-2*(dI0(b)-dI0(a))*(a+b);
        %a13=-4*(Nc-2)*(dI0(b)-dI0(a))*a;
        %a22=2*(ddI0(b)*(b-a)^2-2*(dI0(b)-dI0(a))*b);
        %a23=-4*(Nc-2)*(dI0(b)-dI0(a))*a;
        %a33=4*(Nc-2)*(ddI0(a)*(b-a)^2-(dI0(b)-dI0(a))*(2*b+(Nc-4)*a));

        a11=-2*(gS0(a)-gS0(b)-a*gS0d(a)+b*gS0d(b));
        a12=-(a+b)*(gS0d(b)-gS0d(a));
        a13=(2*b)*(gS0d(b)-gS0d(a));
        a22=-(gS0dd(a)*(a-b)^2+2*a*(gS0d(b)-gS0d(a)));
        a23=(2*b)*(gS0d(b)-gS0d(a));
        a33=-(2*gS0dd(b)*(a-b)^2+(4*a-2*b)*(gS0d(b)-gS0d(a)));

        A=[a11 a12 a13;a12 a22 a23; a13 a23 a33];
        %A=diag([ddI0(a) ddI0(a) ddI0(b) 2*ddI0(b) 2*ddI0(b) 2*ddI0(b)]);
        B=eig(A);
        C=min(B);
        %C=det(A);
        %if C<0 && E00>0.66 && E00<0.72 && T00<E00^3 && T00>tmin+4*abs(E00-0.66)^2
        if C>1e-4
            Id(ke,kt)=1;
        end
        
    end

end

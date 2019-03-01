%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LBdaryNAEdgeTriag: Find lower natural boundary of the 
%                    Nc-podal phase in the edge-triangle model.

% Author: Kui Ren
% Last Updated: 06/16/2016

% Nc: Number of intervals in c. The constraint is: c(1)+c(2)+...+c(Nc)=1.
%
% E0: given value of E

% C : A row vector of length Nc that contains values of the c's
% G : A symmetric matrix of size NcxNc that contains values of the 
%     multipodal graphons

function t=LBdaryNAEdgeTriag(Nc,e0)

ec=zeros(1,30); % starting point of scallop # Nc-1
for k=1:30
	ec(k)=k/(k+1);
end
if Nc==2
    if e0<=0.50
        t=0;
    elseif e0==1.0
        t=1.0;
    else
        P=[-2 8 -3*e0-6 0 2*e0.^2+4*e0 -4*e0.^2 e0.^3];
        X=roots(P);
        C=zeros(1,6);
        T=zeros(1,6);
        K=0;
        for k=1:6
            if imag(X(k))==0 & real(X(k))<1 & real(X(k))>0
                K=K+1;
                C(K)=1-X(k);
                T(K)=(e0-2*C(K)*(1-C(K)))^3/(1-C(K))^3+3*C(K)*(e0-2*C(K)*(1-C(K)));
            end
        end
        [t I]=min(T(1:K));
    end
else
	if e0<=ec(Nc-1)
		t=LBdaryScallopEdgeTriag(e0);
    else
        tS=LBdaryNSEdgeTriag(Nc,e0); % lower boundary of NS phase        
        NMax=10;
        Ng=Nc*(Nc+1)/2;
        Ncg=Nc+Ng; % number of unknowns
        nm=1;
        t=e0^3;
        while t > tS | nm<NMax
		[fval X]=LBdary_Search(e0,Nc,'off','on');
            if fval<t
                t=fval;
            end
            nm=nm+1;            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
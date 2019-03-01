%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EdgeTriagNT:  Finding maximizing multipodal graphons using Newton's 
%               method for the edge-triangle model with edge and triangle 
%               constraints

% Author: Kui Ren
% Last Updated: 06/16/2016

% Nc: Number of intervals in c. The constraint is: c(1)+c(2)+...+c(Nc)=1.
% C : A row vector of length Nc that contains values of the c's

% G : A symmetric matrix of size NcxNc that contains values of the 
%     multipodal graphons

% E : (Normalized) # of edges
% E0: given value of E

% T : (Normalized) # of k-star
% T0: given value of T

% S: The function S(g)=-[g log g-(1-g)log(1-g)]/2 (negative of the true S)
% s: The integral of S(g) over the unit square.
% I: The rate function, i.e. -s, to be minimized.

clear all; close all;

%addpath 'UserDefined/'

tic; tb=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup with parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExitMATLAB='NO'; % Exit MATLAB after calculation? YES or NO

% This line is used to determine whether it is a calculation for fixed t or
% a calculation for fixed e
FixedE = 'YES'; % Fixed E and vary t? 'YES' or 'NO'

Directory='TestCases/'; % Directory to store data

E0=[0.60000000:0.0001:0.60000000];

Nc=8;

NMax=5; % number of local/global maximizers to obtain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ng=Nc*(Nc+1)/2;
Ncg=Nc+Ng; % number of unknowns
NE=length(E0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop on e and t values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KE=1;
while KE<=NE
    
    E00=E0(KE);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup t values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tmax=EREdgeTriag(E00); % ER curve	
    %tmin=LBdaryScallopEdgeTriag(E00); % bottom boundary in scallop part
    %tmin=tmin+39*(tmax-tmin)/40;
    tmax=0.155;
    tmin=0.147;
    dt=(tmax-tmin)/800;
    T0B=tmax:-dt:tmin;
    T0=T0B(1:10);
    %T0=0.5^3-0.1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NT=length(T0);        
    if dt<=1e-7
        disp('The dt is set too small. There are problems with file names!');
        break;            
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine if t values are above natural boundary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t0=LBdaryNAEdgeTriag(Nc,E00);
    if t0>T0(1)
        disp('The t values are all below the natural boundary!');
        break;
    else
        T0BB=T0(1):-dt:t0;
        NTT=length(T0BB);
        if NTT<NT
            NT=NTT;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    KT=1;

    while KT<=NT

        T00=T0(KT);
        E00int=round(1e8*E00);  Echar=charE(E00int);
        T00int=round(1e8*T00);  Tchar=charT(T00int);
        if strcmp(FixedE,'YES')
            CaseName=strcat(Directory,Echar,'/',Echar,'-',Tchar);
        else
            CaseName=strcat(Directory,Tchar,'/',Tchar,'-',Echar);
        end        
        FilesCGET=strcat(CaseName,'-sCGET');        
        FileSummary=strcat(CaseName,'-Summary.txt');
        
        
        E=zeros(1,NMax); T=E; smax=E; C=zeros(NMax,Nc); G=zeros(Nc,Nc,NMax);
        X=zeros(NMax,Ncg); Lambda=zeros(2,NMax);
        
        for nm=1:NMax

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Solve the minimization problem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            disp(' ');
            disp(['Case: ' 'e=' num2str(E00) ';  ' 't=' num2str(T00)... 
                '; ' 'sample #  ' num2str(nm)]);
            disp(' ');
            
            [X(nm,:),fval,lambdaa]=EdgeTriagNT_Search(E00,T00,Nc,'Optim','off','on');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Post-process minimization results
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            smax(nm)=-fval;
            [C(nm,:) G(:,:,nm)]=X2Graphon(X(nm,:),Nc);
            E(nm)=gE(G(:,:,nm),C(nm,:));
            T(nm)=gTTriag(G(:,:,nm),C(nm,:),Nc);
            Lambda(:,nm)=lambdaa.eqnonlin;
            
            clear fval lambdaa;
                       
        end
        
        save(FilesCGET,'E','T','smax','Lambda','C','G');
        
        FF=zeros(NMax,Ncg+5);
        for nm=1:NMax
            FF(nm,:)=[E(nm) T(nm) smax(nm) Lambda(:,nm)' X(nm,:)];
        end
        save(FileSummary,'FF','-ascii');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        clear E T smax Lambda C G X FF;
        
        KT=KT+1;
        
    end
    
    KE=KE+1;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop on e and t values (end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

te=toc;
disp(' ');
disp(' ');
disp(['The code run for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

rmpath 'UserDefined/'

if strcmp(ExitMATLAB,'YES')
    exit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

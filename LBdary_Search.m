%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fval X]=LBdary_Search(e0,Nc,DISP,GRAD)

Ng=Nc*(Nc+1)/2;
Ncg=Nc+Ng; % number of unknowns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set minimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aeq=[ones(1,Nc) zeros(1,Ng)];
beq=1;
LB=zeros(1,Ncg)+eps; UB=ones(1,Ncg)-eps; % upper and lower bounds

%The difference between options = optimoptions('fmincon','algorithm','sqp')
%    and options = optimset('algorithm','sqp') is that the former requires
%    solver name but the later does not.
options = optimset('algorithm','sqp','maxfunevals',5000,'maxiter',100);
options = optimset(options,'tolx',1e-9,'tolcon',1e-9,'tolfun',1e-6);
if strcmp(DISP,'on')
	options = optimset(options,'display','iter-detailed');
else
	options = optimset(options,'display','off');
end
if strcmp(GRAD,'on')
	options = optimset(options,'GradObj','on','GradConstr','on');
elseif strcmp(GRAD,'central')
	options = optimset(options,'FinDiffType','central');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the minimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IND=100; % indicator for rank-deficiency
while IND>=1
    exitflag=100; % indicator for optimization status
    while exitflag~=1
        %X0=rand(1,Ncg); X0=[X0(1:Nc)/sum(X0(1:Nc)) X0(Nc+1:Ncg)];
        X0=rand(1,Ncg); 
        [Ct Gt]=X2Graphon(X0,Nc); Ct=Ct/(sum(Ct)); Et=gE(Gt,Ct);
        X0=[Ct X0(Nc+1:Ncg)/Et*e0];
        clear Ct Gt Et;
        [X,fval,exitflag,output]=fmincon(@(x) LBdary_Obj(x,Nc),X0,...
            [],[],Aeq,beq,LB,UB,@(x) LBdary_Con(x,e0,Nc),options);
        if exitflag==2
            X0=X;
            [X,fval,exitflag,output]=fmincon(@(x) LBdary_Obj(x,Nc),X0,...
                [],[],Aeq,beq,LB,UB,@(x) LBdary_Con(x,e0,Nc),options);
        end
    end
    IND=length(strfind(output.message,'deficient'));
    if IND>=1
        disp('Rank deficiency detected, re-run simulation');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
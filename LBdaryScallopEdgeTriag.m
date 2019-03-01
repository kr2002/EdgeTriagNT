% For a given edge density e, this function finds the t value on the 
% lower boundary of the phase space, i.e. the scallops, following 
% the g_3 function given in the Razborove-CPC08 paper
function t=LBdaryScallopEdgeTriag(e)
if e<=0.5
    t=0.0;
elseif e==1.0
    t=0.0;
else
    p=floor(1/(1-e));
    t=(p-1)*(p-2*sqrt(p*(p-e*(p+1)))).*(p+sqrt(p*(p-e*(p+1)))).^2;
    t=t/(p^2*(1+p)^2);
end

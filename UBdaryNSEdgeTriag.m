% For given edge density e, this function finds the 
% upper natural boundary of the symmetric N-podes 
% phase in the edge-triangle model
function t=UBdaryNSEdgeTrig(N,e)

a=N*e.^3;
b=e.^3+(1-e).^3/(N-1)^2;
t=min(a,b);

end

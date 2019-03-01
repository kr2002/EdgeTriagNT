% For a given edge density e, this function finds
% t on the lower boundary of symmetric N-pod phase 
% in the edge-triangle model

function t=LBdaryNSEdgeTriag(N,e)

a=e.^3-(N-1)*(1-e).^3;
b=e.^3*(1-1/(N-1)^2);
t=max(a,b);

end

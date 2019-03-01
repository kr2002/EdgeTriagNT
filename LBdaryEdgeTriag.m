% LBdaryAllPhase: Find lower boundary of the edge-triagle phase space

% Author: Kui Ren
% Last Updated: 06/16/2016

function t=LBdaryEdgeTriag(e0)

if e0<=0.5
	t=0.0;
else
	t=e0*(2*e0-1);
end

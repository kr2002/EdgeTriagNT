% This function computes Boltzmann entropy
function f=gS0(g)
f=-(g.*log(g)+(1-g).*log(1-g))/2;

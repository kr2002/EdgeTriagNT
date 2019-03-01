% This function compute dS0/dg
function f=gS0d(g)
f=-(log(g)-log(1-g))/2;

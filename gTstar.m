% This function computes k(>=2)-star density
function f=gTstar(g,c,Nc,Kstar)
f=0.0;
d=zeros(1,Nc);
for l=1:Nc
	d(l)=0;
	for m=1:Nc
		d(l)=d(l)+g(l,m)*c(m);
	end
end
for m=1:Nc
	f=f+d(m)^Kstar*c(m);
end

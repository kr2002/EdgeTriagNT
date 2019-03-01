% For a given graphon, this function computes the 
% triangle density
function f=gTTriag(g,c,Nc)
f=0.0;
for l=1:Nc
    for m=1:Nc
        for n=1:Nc
            f=f+g(l,m)*g(m,n)*g(n,l)*c(l)*c(m)*c(n);
        end
    end
end

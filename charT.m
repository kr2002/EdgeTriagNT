% This function converts a integer into a filename
function f=charT(T00int)
if T00int<1e1
    f=strcat('T0p0000000',num2str(T00int));
elseif T00int<1e2
    f=strcat('T0p000000',num2str(T00int));
elseif T00int<1e3
    f=strcat('T0p00000',num2str(T00int));
elseif T00int<1e4
    f=strcat('T0p0000',num2str(T00int));
elseif T00int<1e5
    f=strcat('T0p000',num2str(T00int));
elseif T00int<1e6
    f=strcat('T0p00',num2str(T00int));
elseif T00int<1e7
    f=strcat('T0p0',num2str(T00int));
else
    f=strcat('T0p',num2str(T00int));
end
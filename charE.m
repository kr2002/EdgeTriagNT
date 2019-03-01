% This function convert integer into a filename
function f=charE(E00int)
if E00int<1e1
    Echar=strcat('E0p0000000',num2str(E00int));
elseif E00int<1e2
    f=strcat('E0p000000',num2str(E00int));
elseif E00int<1e3
    f=strcat('E0p00000',num2str(E00int));
elseif E00int<1e4
    f=strcat('E0p0000',num2str(E00int));
elseif E00int<1e5
    f=strcat('E0p000',num2str(E00int));
elseif E00int<1e6
    f=strcat('E0p00',num2str(E00int));
elseif E00int<1e7
    f=strcat('E0p0',num2str(E00int));
else
    f=strcat('E0p',num2str(E00int));
end
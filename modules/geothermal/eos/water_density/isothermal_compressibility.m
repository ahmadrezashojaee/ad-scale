%% Fonction that computes the isothermal compressibility
%% equations (7) and (14) in Spivey et al., 2004

function [ I ] = isothermal_compressibility(p,p0,E,F)

a = E.*(p./p0) + F;
tol = 1e-16;
a(a<tol) = tol;
%I       = 1./E .* log(E.*(p./p0) + F);
I       = 1./E .* log(a);

end


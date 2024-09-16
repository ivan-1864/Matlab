%%%% Make RE and RN %%%%
%
%
% V.Vyazmin, NavLab
% Last revision: 12.11.2015

function [Re, Rn] = Make_Re_Rn(Phi,Alt)
% phi - deg, vector 1 x N or N x 1
% alt - m, vector 1 x N or N x 1
%ESTIM_global_var;
common_c;

SinPhi = sin( deg2rad(Phi) );

dC = 1.0 - GRS_eks*SinPhi.*SinPhi;
dCs = sqrt(dC);

Re = GRS_a ./ dCs + Alt;
Rn = GRS_a*(1.0 - GRS_eks) ./ ( dC.*dCs)  + Alt;

return
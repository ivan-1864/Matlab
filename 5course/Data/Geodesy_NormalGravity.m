% NORMAL GRAVITY 
%
% V.Vyazmin, NavLab, MSU
% Last revision: 15.05.2020

function [g0] = Geodesy_NormalGravity(PhiDegree,h)
% function [g0,gN] = Geodesy_NormalGravity(PhiDegree,h)
% phi, h - geographic lat [deg], height [m]
% g0 = [m/s^2]

%common_c;
%9global GRS_Ge GRS_a
PhiRad = deg2rad(PhiDegree);

 GE=9.78049;
 GRS_a=6378137.0;


%----------------- Normal value at h=0. Helmert formula ----------------
g0 = 9.78030*(1 + 0.005302*(sin(PhiRad)).^2 - 0.000007*(sin(2*PhiRad)).^2)  - 0.00014;  % [m/s2]



%------------------------ Height correction------------------------------------
% Standard free-air correction:
%g_h   =  - 0.3086 * 10^(-5) * h; % m/s/s
W02=GE/GRS_a;
g_h=-2.0*W02*h;
%omega_sh=GRS_Ge/GRS_a;
%g_h   =  - 2*( 1.533439937085077e-06)^2  * h; % m/s/s


% % Torge formula (see TORGE Eq.4.63):
% g_h  =  -( 0.30877 - 0.00043 * (sin(PhiRad)).^2 ).*h * 10^(-5); % Linear term
% g_h  =  g_h  + 0.72*10^(-12)*(h.*h);  % Quadratic term 



%------------------ Normal gravity at height h --------------------
g0   = g0 + g_h;   % m/s/s


return
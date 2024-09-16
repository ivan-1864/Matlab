%%%% Orientation matrix  A_xi_eta: eta --> xi   %%%%
% 
%
% V.Vyazmin, NavLab
% Last revision: 18.08.2017

function [ A_xi_eta ] = Make_A_xi_eta_Matrix( t )

    common_c;

    CosUT = cos(GRS_omega*t);
    SinUT = sin(GRS_omega*t);

	A_xi_eta = zeros(3);
	
	A_xi_eta(1,1) = CosUT;
	A_xi_eta(2,1) = SinUT;

	A_xi_eta(1,2) = -SinUT;
	A_xi_eta(2,2) =  CosUT;

	A_xi_eta(3,3) = 1.0;

return;

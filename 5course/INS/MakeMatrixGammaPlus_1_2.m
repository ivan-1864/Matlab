% V.Vyazmin, NavLab, MSU
% Last revision: 18.08.2017

function [Matr] = MakeMatrixGammaPlus_1_2(gamma)
% gamma -  1x3, [rad]

    C1 = 0.5;
    C2 = 0.125; 
	
Matr = zeros(3);

GammaModule = sqrt( gamma*gamma' );

if GammaModule > 1.0e-14 	       
    C1 = sin( 0.5*GammaModule )/GammaModule; 	
    C2 = (1.0 - cos( 0.5*GammaModule ) )/( GammaModule*GammaModule );	
end

%---------Rodriguez rotation formula------------
Matr(1,1) = 1.0 + C2*( -gamma(2)*gamma(2) - gamma(3)*gamma(3) );
Matr(2,2) = 1.0 + C2*( -gamma(1)*gamma(1) - gamma(3)*gamma(3) );
Matr(3,3) = 1.0 + C2*( -gamma(1)*gamma(1) - gamma(2)*gamma(2) );

Matr(1,2) =   C1*gamma(3) + C2*gamma(1)*gamma(2);	
Matr(2,1) =  -C1*gamma(3) + C2*gamma(1)*gamma(2);

Matr(1,3) =  -C1*gamma(2) + C2*gamma(1)*gamma(3);
Matr(3,1) =   C1*gamma(2) + C2*gamma(1)*gamma(3);

Matr(2,3) =   C1*gamma(1) + C2*gamma(2)*gamma(3);
Matr(3,2) =  -C1*gamma(1) + C2*gamma(2)*gamma(3);


return;

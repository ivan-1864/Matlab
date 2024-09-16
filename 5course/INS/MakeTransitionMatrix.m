
% V.Vyazmin, NavLab, MSU
% Last revision: 18.08.2017

function [Matr] = MakeTransitionMatrix(gamma)
% gamma -  1x3, [rad]

Matr = zeros(3);

%OmegaModule = sqrt( gamma*gamma' );
% if OmegaModule > 1.0e-10 	       
%     C1 = sin( OmegaModule*TimeStep )/OmegaModule; 	
%     C2 = (1.0 - cos( OmegaModule*TimeStep ) )/( OmegaModule*OmegaModule );	
% else    
%     C1 = TimeStep;	
%     C2 = 0.5*TimeStep*TimeStep;
% end

GammaModule = sqrt( gamma*gamma' );
if GammaModule > 1.0e-14 	       
    C1 = sin( GammaModule )/GammaModule; 	
    C2 = (1.0 - cos( GammaModule ) )/( GammaModule*GammaModule );	
else    
    C1 = 1.0;	
    C2 = 0.5;
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

return

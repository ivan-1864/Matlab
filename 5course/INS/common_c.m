%     INPUT PARAMETERS FOR THE DESIRED GEODETIC REFERENCE SYSTEM
%     BASICALLY WGS-84
%
% Last update: 27.08.2019

global GRS_C20 GRS_C40 GRS_C60 GRS_C80 GRS_C100 %GRS_J20 
global GRS_a GRS_gm GRS_omega GRS_finv GRS_eks GRS_Ge GRS_Gp


% WGS-84 ellipsoid:     
       GRS_a        =  6378137.0d0;     % m, major semiaxis, meters, WGS-84 
       GRS_gm       =  3.986004418d14;  % m3/s^2       % GM  
       GRS_finv     =  298.257223563;   % flattening of WGS 84 ellipsoid, 1/finv = (a-b)/a 
       GRS_eks      =  6.6943799901413e-3;% WGS84 ellipsoide eks squared */        
       GRS_omega    =  7292115.d-11;    % earth angular velocity, rad/sec       
       
       GRS_Ge       =  9.7803253359; % Normal gravity at the equator, m/s/s
       GRS_Gp       =  9.8321849378; % Normal gravity at the North pole, m/s/s
       
            
% WGS-84 Normal gravity harmonic coefficients
       GRS_C20  = -0.484166774985d-03;
       GRS_C40  =  0.790303733511d-06;
       GRS_C60  = -0.168724961151d-08;
       GRS_C80  =  0.346052468394d-11;
       GRS_C100 = -0.265002225747d-14;                    
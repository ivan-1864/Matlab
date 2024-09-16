%%%% Orientation matrix  A_eta_x: x --> eta   %%%%
% 
%
% V.Vyazmin, NavLab
% Last revision: 18.08.2017

function [ A_eta_x ] = Make_A_eta_x_Matrix( Lat, Lon, Eps ) 

	   SinEps = sin(Eps); 
	   CosEps = cos(Eps); 
	   SinLat = sin(Lat); 
	   CosLat = cos(Lat); 
	   SinLon = sin(Lon); 
	   CosLon = cos(Lon);

A_eta_x = zeros(3);

A_eta_x(1,1) =  -SinLon*CosEps - CosLon*SinLat*SinEps;
A_eta_x(2,1) =  -SinLat*SinLon*SinEps + CosLon*CosEps;
A_eta_x(3,1) =   CosLat*SinEps;

A_eta_x(1,2) =   SinLon*SinEps - CosLon*SinLat*CosEps;
A_eta_x(2,2) =  -SinLat*SinLon*CosEps - CosLon*SinEps;
A_eta_x(3,2) =   CosLat*CosEps;

A_eta_x(1,3) =   CosLat*CosLon;
A_eta_x(2,3) =   CosLat*SinLon;
A_eta_x(3,3) =   SinLat;

return;
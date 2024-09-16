function [Matr] = Make_C_Matrix(Lat, u, h, Re, Rn, v1, v2)

Matr = zeros(3);

Matr(1,1)=0;
Matr(2,2)=0;
Matr(3,3)=0;

Matr(1,2)=(v1*tan(Lat))/(Re+h)+2*u*sin(Lat);
Matr(1,3)=-(v1/(Re+h))+2*u*cos(Lat);
Matr(2,3)= -v2/(Rn+h);

Matr(2,1)=-((v1*tan(Lat))/(Re+h)+2*u*sin(Lat));
Matr(3,1)=-((v1/(Re+h))+2*u*cos(Lat));
Matr(3,2)= v2/(Rn+h);

return;
function [ L_zx ]= Make_L_zx_Matrix(psi, gamma, nu) 
L_zx=zeros(3);

L_zx(1,1)=cos(psi)*cos(gamma)+sin(psi)*sin(nu)*sin(gamma);
L_zx(1,2)=-sin(psi)*cos(gamma)+cos(psi)*sin(nu)*sin(gamma);
L_zx(1,3)=-cos(nu)*sin(gamma);
L_zx(2,1)=sin(psi)*cos(nu);
L_zx(2,2)=cos(psi)*cos(nu);
L_zx(2,3)=sin(nu);
L_zx(3,1)=cos(psi)*sin(gamma)-sin(psi)*sin(nu)*cos(gamma);
L_zx(3,2)=-sin(psi)*sin(gamma)-cos(psi)*sin(nu)*cos(gamma);
L_zx(3,3)=cos(nu)*cos(gamma);

return
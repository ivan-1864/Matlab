function[Matr]=Make_Cross_Matrix (v);

Matr=zeros(3);

Matr(1,1)=0;
Matr(2,2)=0;
Matr(3,3)=0;

Matr(1,2)=v(3);
Matr(1,3)=-v(2);
Matr(2,3)=v(1);

Matr(2,1)=-v(3);
Matr(3,1)=v(2);
Matr(3,2)=-v(1);

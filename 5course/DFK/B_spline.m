function [b_spline]=B_spline(t) % t принимает ненулевые значения от 0 до 4
b_spline=0;
if t>=0 && t<1 
    b_spline=(t^3)/6;
end
if t>=1 && t<2
    b_spline=(4-3*t*(t-2)^2)/6;
end
if t>=2 && t<3
    b_spline=(4+3*(t-4)*(t-2)^2)/6;
end  
if t>=3 && t<4
    b_spline=((4-t)^3)/6;
end
   return

%%%% Kalman Filter %%%%
%    using Square-Root Information form (see T.Kailath)
%
%
% V.Vyazmin, NavLab
% Last revision: 03.12.2016

function [Y_pred,P_pred,Tmp1,Tmp2,Residual] = KF_forward(F_t,J_t,Q_sqrt,Z_t,H_t,R_sqrt,R_sqrt_inv,Y_last,P_last)
% SYSTEM: 
% x(t+1) = F(t)*x(t) + J(t)*q(t)
% z(t) = H(t)*x(t) + r(t)
%   M[q(t)*q'(t)] = Q
%   M[r(t)*r'(t)] = R
%   M[q(t)*r'(t)] = 0
% INPUT:
%   Y_last - square root inform. vector at t-1, Nx1
%   P_last - square root inform. matrix at t-1, NxN  (upper triang) 
%   Q_sqrt - upper triang. square root of Q, qxq 
%   R_sqrt - upper triang. square root inverse of R, rxr 
%   Z_t    - measurements at t, rx1
N        =  size(P_last,1);
nq       =  size(Q_sqrt,1); 
nr       =  size(R_sqrt,1); 


%------ Construct matrix and make QR----------------
CompMatr1      = [ R_sqrt,  H_t*P_last, zeros(nr,nq); 
                   zeros(N,nr), F_t*P_last,  J_t*Q_sqrt];

[Unitar,Upper] = qr(CompMatr1');


CompMatr = [ CompMatr1;             
             zeros(N,nr), eye(N), zeros(N,nq);
             Z_t'*R_sqrt_inv, -Y_last', zeros(1,nq)];
        
Res      =   CompMatr * Unitar; % Upper Traingular complemented


Tmp1  = Res(nr+N+1:nr+2*N,1:nr);      %  N x nr
Tmp2  = Res(nr+N+1:nr+2*N,nr+1:nr+N); %  N x N

Residual =  Res(end,1:nr)';           %  nr x 1

% Predicted results
Y_pred   = -Res(end,nr+1:nr+N)';      %  N x 1
P_pred   =  Res(nr+1:nr+N,nr+1:nr+N); %  N x N

end
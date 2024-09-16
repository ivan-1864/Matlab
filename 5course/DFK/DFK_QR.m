% close all
clear

timetime=clock;
% --------------загрузка данных--------------------------------------------
D               = importdata('../../INS/Output_data/output_INS_data.txt');  
time            = D.data(:, 1);
dVx_array       = D.data(:, 2:4);
Vx_array        = D.data(:, 5:7);
Omega_x_array   = D.data(:, 8:10);
u_E_array       = D.data(:, 23:25);
L_zx_array      = [D.data(:, 11),D.data(:, 12),D.data(:, 13),...
                   D.data(:, 14),D.data(:, 15),D.data(:, 16),...
                   D.data(:, 17),D.data(:, 18),D.data(:, 19)];
g0_array        = D.data(:, 20:22); 
omega_x_array   = Omega_x_array+u_E_array;
clear D

D = importdata('../../Data/Output_data/GPS_data.txt');
lon = D.data(:, 2);
lat = D.data(:, 3);
clear D

cfg;
% ---------------глобальные параметры--------------------------------------
Length          = size(time,1);
time_step       = mean(diff(time));
% -------------------------------------------------------------------------


lat_start=min(lat);
lat_end=max(lat);
% period_spline=Length/(num_spline-3); % кол-во строк, соответвующее одному сплайну
koef=(lat_end-lat_start)/(num_spline-3);% коэффициент пропорциональности м/у сплайнами и широтой


% ----------создание матриц P, Q, R----------------------------------------
d_dv=(Std_dV)*ones(1,3);
d_beta=(deg2rad(Std_beta/60))*ones(1,3);
d_f=[(Std_f*10^(-5))*ones(1,2), 10^(-5)];
d_nu=(deg2rad(Std_nu)/3600)*ones(1,3);
d_c=(Std_c)*ones(1,3*num_spline);
P_last_diag     = [d_dv,d_beta,d_f,d_nu,d_c]; % diagonal of sqrt of cov. matrix at time t=0
P_last          = diag(P_last_diag);
% -------------------------------------------------------------------------
d_ax            = (10^(-5)*Std_AX)*ones(1,3);
d_dus           = (deg2rad(Std_DUS)/3600)*ones(1,3);
Q_sqrt_diag     = [d_ax, d_dus];
Q_sqrt          = diag(Q_sqrt_diag);
% Q_sqrt          = Q_sqrt*(time_step)^(0.5);%приведение на нужную частоту
% -------------------------------------------------------------------------
d_vel           = (Std_GPS)*ones(1,3);
R_sqrt_diag     = d_vel;
R_sqrt          = diag(R_sqrt_diag);
R_sqrt_inv      = diag(1./R_sqrt_diag);

% --------------------------вектор состояния и "память"--------------------
dimX            = 3*num_spline+12;
Y_last          = zeros(dimX,1);
%----------------- Declare arrays for results ---------------------- 
Y_forw          = zeros(TimeEnd, dimX);
P_forw          = cell(1,TimeEnd);
P_forw{1}       = P_last;

% ------------------алгоритм ДФК-------------------------------------------
disp('начало работы алгоритма')
% for i = 1:Length-1
for i = 1:TimeEnd
%    --------------инициализация-------------------------------------------
    g_0         = Make_Cross_Matrix(g0_array(i,:));
    Vx          = Make_Cross_Matrix(Vx_array(i, :));
    Z_t         = dVx_array(i, :)';
    omega_x     = Make_Cross_Matrix(omega_x_array(i, :));
    omega_pl_u  = Make_Cross_Matrix(Omega_x_array(i, :))+2*Make_Cross_Matrix(u_E_array(i, :));
%     ---------------------------------------------------------------------
    L_zx        = [L_zx_array(i,1), L_zx_array(i,2), L_zx_array(i,3);
                   L_zx_array(i,4), L_zx_array(i,5), L_zx_array(i,6);
                   L_zx_array(i,7), L_zx_array(i,8), L_zx_array(i,9)];
    
%     ---------------------------------------------------------------------
    B           = zeros(3, dimX-12);
    time_spline=(lat(i)-lat_start)/koef;  
    B1          = B_spline(time_spline-floor(time_spline)+3)*eye(3);
    B2          = B_spline(time_spline-floor(time_spline)+2)*eye(3);
    B3          = B_spline(time_spline-floor(time_spline)+1)*eye(3);
    B4          = B_spline(time_spline-floor(time_spline))*eye(3);
    if time_spline ~= num_spline-3 % костыль
        B_diag      = [B1, B2, B3, B4];
        currInterval= 3*floor(time_spline)+1:3*floor(time_spline)+12;
    else
        B_diag=[B1, B2, B3];
        currInterval =3*floor(time_spline)+1:3*floor(time_spline)+9;
    end
    B(:,currInterval) = B_diag;
%     ---------------------------------------------------------------------
     A         = [omega_pl_u, g_0, L_zx', Vx*L_zx', -B;
                   zeros(3), omega_x, zeros(3), L_zx', zeros(3, dimX-12);
                   zeros(dimX-6, dimX)];
    H_t         = [eye(3), -Vx, zeros(3, dimX-6)];
   F_t=eye(size(A))+A*time_step;
%     ---------------------------------------------------------------------
    J_t         =[L_zx', Vx*L_zx';
                  zeros(3,3),L_zx';     
                  zeros(dimX-6,6)];
    
    [Y_pred,P_pred,Tmp1,Tmp2,Resid] = KF_forward(F_t,J_t,Q_sqrt,Z_t,H_t,R_sqrt,R_sqrt_inv,Y_last,P_last);
    Y_last = Y_pred; 
    P_last = P_pred;  
    
    % Copy into arrays
    Y_forw(i+1,:) = Y_pred';% 1 x Nx
    P_forw{i+1}   = P_pred;


%    -----------------запись-----------------------------------------------
    if i/100-floor(i/100)==0 
        disp(['i=', num2str(i), ' of ', num2str(Length), ' (' , num2str(i/Length*100), '%)']);
    end
end
%----------Estimates from Forward KF ------------------------------------
X_j=zeros(dimX, TimeEnd);
for j = 1 : TimeEnd
       
    Y_current = Y_forw(j, :)';    
    P_current = P_forw{j};    
    
    P_j  = P_current * P_current';  % Covariance matrix    
    X_j(:, j)  = P_current * Y_current;   % State vector estimate  
    if j/1000-floor(j/1000)==0 
        disp(['i=', num2str(j), ' of ', num2str(Length), ' (' , num2str(j/Length*100), '%)']);
    end
end
% ---------------------конец работы алгоритма ДФК--------------------------

% ---------------------запись данных в файл--------------------------------

disp('запись оценок')
File = fopen('../Output_data/DFK_estimation_QR.txt', 'w');
fprintf(File,'%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s \n','time', 'delta V_1','delta V_2','delta V_1',...
         'betta1', 'betta2', 'betta3', 'Delta f_1', 'Delta f_2', 'Delta f_3', 'nu1', 'nu2', 'nu3');
for i =1:TimeEnd
   fprintf(File,'%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f \n',time (i, 1), X_j(1:12, i));
end

disp('запись коэффициентов')
File = fopen('../Output_data/spline_coef_QR.txt', 'w');
for i =13:length(X_j(:, TimeEnd))
    fprintf(File, '%20.10f\n', X_j(i, TimeEnd));
end

File = fopen('../Output_data/Log.txt', 'w');
    fprintf(File,'%20s %20s \n', 'koef', 'num_spline');
    fprintf(File, '%20.10f %20.10f \n', koef, num_spline);

disp(['время работы программы: ', num2str(etime(clock, timetime)), ' секунд'])
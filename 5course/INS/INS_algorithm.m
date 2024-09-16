% описание
% водные параметры: 
% gps: время1, координаты, скорости
% инс: время2, ориентация, удельные силы, угловые скорости
% данные поступают в цикле по галсам. 
% все данные приводятся ко времени инс
% вычисляется приборная скорость как результат работы алгоритма ИНС (как численное решение модельных диф уравнений)
% полученные скорости прореживаются на частоту gps
% вычисляются ошибка скоростей и ошибка ускорений (численным диф-ем)
% получившиеся данные записываются в разные файлы (по галсам)


close all
clear
timetime=clock;
clc;
common_c;
global GRS_omega %#ok<*GVMIS> 

D = importdata('../../Data/Output_data/IMU_data.txt');    
IMU_data = D.data;
clear D;

D = importdata('../../Data/Output_data/GPS_data.txt');
GPS_data=D.data;
clear D;


[Length_IMU,Width_IMU]=size(IMU_data);
[Length_GPS,Width_GPS]=size(GPS_data);
ImuTimeStep=mean(diff(IMU_data(:,1)),1);     
GpsTimeStep=mean(diff(GPS_data(:,1)),1);
method  = 'spline';
GPS_data_interp = interp1(GPS_data(:,1), GPS_data, IMU_data(:,1) ,method);

plot(GPS_data(:, 1), GPS_data(:, 3));
%--------- Initial coord. & velocity ----------------------------------
MHei = GPS_data_interp(:,4);
MLat = GPS_data_interp(:,3);
MLon = GPS_data_interp(:,2);
MVe  = GPS_data_interp(:,5);
MVn  = GPS_data_interp(:,6);
%--------- Precompute Angular Velocities ---------------------------------------------
u_E       =   GRS_omega * [zeros(size(MLat)),cos(deg2rad(MLat)),sin(deg2rad(MLat))]; % [rad/s], Length x 3
[Re,Rn]   =   Make_Re_Rn( MLat, MHei ); 
Omega_x   =  [-MVn./Rn, MVe./Re, MVe.*sin(deg2rad(MLat))./cos(deg2rad(MLat))./Re];% [rad/s], relative ang. vel
%------------- IMU Data -----------------
Fs            = IMU_data(:,5:7);%arr(:, 11:13);  Nx3, acc 
omega_s       = deg2rad([IMU_data(:,8),  IMU_data(:,9), IMU_data(:,10) ])/3600; % Nx3, gyros [deg/h]->[rad/s]
MHeading_t0   = 0;
MRoll_t0      = 0;
MPitch_t0=0;
%--------- Initialize at t0 ----------------------------------
%Bsx_t0         =  Make_Bsx_Matrix( MHeading_t0, MRoll_t0, MPitch_t0 ); % Ideal Alignment, Mx->Ms
Bsx_t0         =eye(3);
Bsx_model_last = Bsx_t0;
g_array        = zeros(Length_IMU,3);
Bsx_array      =  zeros(Length_IMU, 9);
Bsx_array(1,:) =  [Bsx_t0(1,1:3),Bsx_t0(2,1:3),Bsx_t0(3,1:3)];
Fx_array       = zeros(Length_IMU,3);
Fx_array(1,:)  = IMU_data(1, 5:7)*Bsx_t0;
V_array        =  zeros(Length_IMU,3);
V_array(1,:)   =  GPS_data(1,5:7);
V_last = GPS_data(1,5:7)';
%--------------------- MAIN CYCLE ------------------------------------
for i = 1 : Length_IMU-1

if i/5000 - floor(i/5000) == 0 
  disp(['i = ',num2str(i), ' of ', num2str(Length_IMU), ' (' , num2str(i/Length_IMU*100), '%)']);
end
%------------- Project Fs from Ms to Mx ---------------------------------
% Make gamma_z:
	gamma_z = omega_s(i+1,:)*ImuTimeStep;     % At i+1
%MakeMatrixGammaMinus_1_2:
    C_Minus_1_2 = MakeMatrixGammaMinus_1_2(gamma_z); %computing matrix for solving Puason eq
%MatrixVectorProduct:
    w = C_Minus_1_2 * Fs(i+1,:)';  % At i+1 -   
%Make_A_x_s_Matrix:
	A_s_x_0 = Bsx_model_last;
%     A_s_x_0 = Make_Bsx_Matrix(GyroHeading,Roll,Pitch);
    A_x_s_0 = A_s_x_0'; 
%MatrixVectorProduct:
    p = A_x_s_0 * w; 
% Make gamma_x:
	gamma_x =(Omega_x(i,:) + u_E(i,:))*ImuTimeStep;   
%MakeMatrixGammaPlus_1_2:
    C_Plus_1_2 = MakeMatrixGammaPlus_1_2( gamma_x );   
%MatrixVectorProduct:
    Fx = C_Plus_1_2 * p; %- находим силу 
% Ideal Variant:   
   % Fx = Bsx_model_last' * Fs(i+1,:)';  % at i+1  
% Copy into large array:
   Fx_array(i+1,:) = Fx; % in Mx
   


%------------------- Transition Matrix Bsx -----------------------
A_eta_x_0 = Make_A_eta_x_Matrix(deg2rad(MLat(i)),deg2rad(MLon(i)),0);
A_x_eta_0 = A_eta_x_0';
   Gamma_x = Omega_x(i,:)*ImuTimeStep;
    W = MakeTransitionMatrix( Gamma_x );
   A_x_eta_1 = W * A_x_eta_0;    
% Make angle Chi:
    %Chi = atan(  A_x_eta_1(1,3)/A_x_eta_1(2,3) );
	%Lon = atan(  A_x_eta_1(3,2)/A_x_eta_1(3,1) );
	%Lat = asin(   A_x_eta_1(3,1) );	           
A_s_eta_0 = A_s_x_0 * A_x_eta_0;
A_s_xi_0 = A_s_eta_0;
W = MakeTransitionMatrix(gamma_z);
A_s_xi_1 = W * A_s_xi_0;
A_xi_eta =  Make_A_xi_eta_Matrix(ImuTimeStep);
A_s_eta_1 =  A_s_xi_1 * A_xi_eta;
A_eta_x_1 = Make_A_eta_x_Matrix(deg2rad(MLat(i+1)),deg2rad(MLon(i+1)),0);
A_s_x_1 =  A_s_eta_1 * A_eta_x_1; 

% OUTPUT - Transition matrix Bsx at (i+1)-th time stamp:   
Bsx_model_t = A_s_x_1;
%C=zeros(3);  
g_0=[0, 0, -Geodesy_NormalGravity(MLat(i),MHei(i))];
Ang_vel=Make_Cross_Matrix(Omega_x(i,:)+2*u_E(i,:));
V_array(i+1,:)=(V_last+((Ang_vel*V_last+g_0'+Fx)*ImuTimeStep))';


%---------Reinitialize cycle-----------------------------------------------
V_last=V_array(i+1,:)'; 
Bsx_model_last  =  A_s_x_1;
  % Copy into large array:
  g_array(i, :)=g_0;
  Bsx_array(i+1,:)=  [ A_s_x_1(1,1),A_s_x_1(1,2),A_s_x_1(1,3), ...
      A_s_x_1(2,1),A_s_x_1(2,2),A_s_x_1(2,3),...
      A_s_x_1(3,1),A_s_x_1(3,2),A_s_x_1(3,3), ];
end
%конец алгоритма ИНС

tsin = timeseries(V_array,IMU_data(:,1));
tsout = resample(tsin,GPS_data(:,1));
new_V_array =tsout.Data;


dV=new_V_array-GPS_data(:,5:7);%ошибка скоростей проверить!!!

% omega_x=Omega_x+u_E;

data_100Hz=[Omega_x,Bsx_array, g_array, u_E];

tsin = timeseries(data_100Hz,IMU_data(:,1)); 
tsout = resample(tsin,GPS_data(:,1));
data_10Hz =tsout.Data;

output_data=[GPS_data(:, 1), dV, new_V_array, data_10Hz];


 F = fopen('../Output_data/output_INS_data.txt', 'w');
    fprintf(F,'%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n','time[s]','dV1[m/s]','dV2[m/s]','dV3[m/s]','V1[m/s]','V2[m/s]','V3[m/s]', ...
        'omega_x1[rad/s]','omega_x2[rad/s]','omega_x3[rad/s]', ...
        'L_zx11','L_zx12','L_zx13','L_zx21','L_zx22','L_zx23','L_zx31','L_zx32','L_zx33', ...
        'g0_1','g0_2','g0_3'); 
 for i=1:Length_GPS
    fprintf(F,'%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n', ...
        output_data(i, :));
 end

disp(['время работы программы: ', num2str(etime(clock, timetime)), ' секунд'])
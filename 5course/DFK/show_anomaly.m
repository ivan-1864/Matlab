GRS_a        =  6378.137;     % m, major semiaxis, kilometers, WGS-84 
TimeEnd1stGals=29000;

D=importdata('../../Data/Output_data/anomaly_data.txt');
true_anomaly=D.data;
clear D

D=importdata('../Output_data/DFK_estimation.txt');
X=D.data;
clear D

D=importdata('../Output_data/DFK_estimation_QR.txt');
X_QR=D.data;
clear D


D=importdata('../Output_data/Log.txt');
koef=D.data(1);
num_spline=D.data(2);
clear D

Length=length(X_QR);

D=importdata('../Output_data/spline_coef.txt');
x=D;
clear D

D=importdata('../Output_data/spline_coef_QR.txt');
x_QR=D;
clear D

D = importdata('../../Data/Output_data/GPS_data.txt');
% lon = D.data(:, 2);
lat = D.data(:, 3);
clear D

%beta

h=figure(1);
hold
% plot(X(:, 1),rad2deg(X(:, 5:7))*3600, '--')
plot(X_QR(:, 1),rad2deg(X_QR(:, 5:7))*3600, 'LineWidth', 1.25)
title('Ошибки определения ориентации')
legend('\beta_1','\beta_2','\beta_3');
xlabel('time, [s]')
ylabel('beta, [arg.sec]')
% saveas(h,'beta','fig')

 
% delta f
h=figure(2);
hold
plot(X_QR(:, 1),10^(5)*X_QR(:, 8:10), 'LineWidth', 1.25)
title('Оценки смещения нулей АКС')
legend('\Delta_{f_1}','\Delta_{f_2}','\Delta_{f_3}');
xlabel('time, [s]')
ylabel('\Delta f^0, [mGal]')
% saveas(h,'delta_f','fig')

% nu
h=figure(3);
hold
plot(X_QR(:, 1),rad2deg(X_QR(:, 11:13))*3600, 'LineWidth', 1.25)
title('Оценки смещения дрейфа ДУС')
legend('\nu_1','\nu_2','\nu_3');
xlabel('time, [s]')
ylabel('\nu^0, [deg/h]')
% saveas(h,'nu','fig')

lat_start = min(lat);

approx_anomaly=zeros(Length, 3);
approx_anomaly_QR=zeros(Length, 3);
for i = 1:Length
    B=zeros(3, 3*num_spline);
    time_spline=(lat(i)-lat_start)/koef; 
    B1=B_spline(time_spline-floor(time_spline)+3)*eye(3);
    B2=B_spline(time_spline-floor(time_spline)+2)*eye(3);
    B3=B_spline(time_spline-floor(time_spline)+1)*eye(3);
    B4=B_spline(time_spline-floor(time_spline))*eye(3);
    if abs(time_spline-num_spline+3)>10^(-7) % костыль
        B14=[B1, B2, B3, B4];
        B(:,3*floor(time_spline)+1:3*floor(time_spline)+12)=B14;
    else
        B14=[B1, B2, B3];
        B(:,3*floor(time_spline)+1:3*floor(time_spline)+9)=B14;
    end
    approx_anomaly_QR(i, :)=(B*x_QR)';
end

sGals=deg2rad(lat(1:TimeEnd1stGals))*GRS_a;
% 
% dg1
h=figure(4);
hold
plot(sGals, 10^5*true_anomaly(1:10:10*TimeEnd1stGals,2), 'LineWidth', 1.25);
plot(sGals, 10^5*approx_anomaly_QR(1:TimeEnd1stGals,1), 'LineWidth', 1.25);
title('УОЛ (восточная компонента)');
legend('true \delta g_1', 'approx \delta g_1');
xlabel('s, [km]')
ylabel('dg1, [mGal]')
% saveas(h,'anomaly_1','fig');
% hold off

% dg1
h=figure(5);
hold
title('УОЛ (северная компонента)');
plot(sGals, 10^5*true_anomaly(1:10:10*TimeEnd1stGals,3), 'LineWidth', 1.25);
plot(sGals, 10^5*approx_anomaly_QR(1:TimeEnd1stGals,2), 'LineWidth', 1.25);
legend('true \delta g_2', 'approx \delta g_2');
xlabel('s, [km]')
ylabel('dg2, [mGal]')


% dg1


h=figure(6);
hold
plot(sGals, 10^5*true_anomaly(1:10:10*TimeEnd1stGals,4), 'LineWidth', 1.25);
plot(sGals, 10^5*approx_anomaly_QR(1:TimeEnd1stGals,3),'LineWidth', 1.25);
title('УОЛ (вертикальная компонента)');
legend('true \delta g_3', 'approx \delta g_3');
xlabel('s, [km]')
ylabel('dg3, [mGal]')



% dg1
h=figure(8);
hold
title('Ошибка определения УОЛ');
plot(sGals, abs(10^5*(true_anomaly(1:10:10*TimeEnd1stGals,2)-approx_anomaly_QR(1:TimeEnd1stGals,1))),  'LineWidth', 1.25);
plot(sGals, abs(10^5*(true_anomaly(1:10:10*TimeEnd1stGals,3)-approx_anomaly_QR(1:TimeEnd1stGals,2))), 'LineWidth', 1.25);
xlabel('s, [km]')
ylabel('dg, [mGal]')
legend('dg1', 'dg2')
% saveas(h,'anomaly_2','fig');
% hold off
% 
% dg1


h=figure(9);
hold
plot(sGals, abs(10^5*(true_anomaly(1:10:10*TimeEnd1stGals,4)-approx_anomaly_QR(1:TimeEnd1stGals,3))), 'LineWidth', 1.25 );
title('Ошибка определения аномалии (вертикальная компонента)');
xlabel('s, [km]')
ylabel('dg3, [mGal]')
% saveas(h,'anomaly_3','fig');
% 


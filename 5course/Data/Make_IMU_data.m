
close all
clear
timetime=clock;

cfg

% ---------показания �?НС----------------------------------------------------
% загружаем идеальные данные
D               = importdata('../Input_data/Trajectory.txt');    
arr             = D.data;
clear D;
Length          = size(arr,1);


if Add_Anomal==1
% -------аномалия-----------------------------------------------------------
    AnomalyAll=zeros(Length, 3);
    fid  = fopen('../Input_data/anomaly/XGM2019e_2159_2100.dat'); 
    data = textscan(fid,'%f %f %f %f %f %f %f','HeaderLines', 34);
    data = cell2mat(data);
    fclose(fid);

% ------перевод аномалии на нужную частоту и промежуток--------------------


    AnomalyAll_x         = zeros(Length,4);
    NumGals         = 7; 
    
    Finish1stGals=16130.05;
    Finish2ndGals=18549.85;
    Time2ndGals=find(data(:, 1)>=Finish1stGals & data(:, 1)<Finish2ndGals);
    
    anomaly1gals    = [data(Time2ndGals, 6),data(Time2ndGals, 7),data(Time2ndGals, 5)] ;
    
    anomaly1gals    =  anomaly1gals-anomaly1gals(1, :); %1st record eq to zero
    
    new_time_anomal=linspace(data(Time2ndGals(1), 1), data(Time2ndGals(end), 1), 237500)';
    InterpAnomal = interp1(data(Time2ndGals, 1), anomaly1gals, new_time_anomal, 'spline');
    
    for i = 1:NumGals
    
        StartImu=650+(i-1)*2375;
        FinishImu=3025+(i-1)*2375;
        
        TimeImu=find(arr(:, 1)>StartImu & arr(:, 1) <= FinishImu);
        
        if i/2-floor(i/2)==0
            AnomalyAll(TimeImu, :)    = flipud(InterpAnomal);
        %     inverse array 
        else
            AnomalyAll(TimeImu, :)    = InterpAnomal;
        end
    end


for i = 1:Length
   L_zx         = Make_L_zx_Matrix(deg2rad(arr(i, 10)),deg2rad(arr(i, 9)),deg2rad(arr(i, 8)));
   AnomalyAll(i, 1)   = -deg2rad(AnomalyAll(i, 1)/3600)*Geodesy_NormalGravity(arr(i, 3), arr(i, 4)); %arcsec->m/s/s
   AnomalyAll(i, 2)   = -deg2rad(AnomalyAll(i, 2)/3600)*Geodesy_NormalGravity(arr(i, 3), arr(i, 4)); %arcsec->m/s/s
   AnomalyAll(i, 3)   = 10^(-5)*AnomalyAll(i, 3); %mGal->m/s/s
   AnomalyAll_x(i, :) = [arr(i, 1), AnomalyAll(i, :)];
   AnomalyAll(i, :)   = AnomalyAll(i, :)*L_zx'; %Mx->Mz
end

end


if Add_Bias==1
    Delta_AX        = Add_Bias*[30,-40,0]*10^(-5); 
    Delta_DUS       = Add_Bias*[-3, 3, 1]*10^(-3); % deg/h

end

if Add_Noise==1
    Std_DUS         = 0.3; %[deg/h]  100Hz
    Std_Acc         = 30*10^(-5); % [m/s/s] 100Hz
    tstep           = 0.01; %timestep
% -------------------------------------------------------------------------


    Std_DUS = Std_DUS*(tstep/2)^(0.5);
    
    rand1=randn(1,Length+1)';

    WxNoise=Std_DUS*(rand1(2:end)-rand1(1:end-1))./tstep; %!!!! this is bad!
    rand2=randn(1,Length+1)';
    WyNoise=Std_DUS*(rand2(2:end)-rand2(1:end-1))./tstep; %std=1*10^(-4) rad/s = 2.1 deg/h
    
    rand3=randn(1,Length+1)';
    WzNoise=Std_DUS*(rand3(2:end)-rand3(1:end-1))./tstep; %model x=a\dot{q}, 100 Hz

    AxNoise=Std_Acc*randn(1, Length)'; %std=30 mGal, model x=aq
    AyNoise=Std_Acc*randn(1, Length)'; %std=30 mGal, model x=aq, 100 Hz
    AzNoise=Std_Acc*randn(1, Length)'; %std=30 mGal, model x=aq, 

end

% adding the errors in the data
IMU_data=zeros(Length,10);
IMU_data(:,1)   = arr(:,1); %time
IMU_data(:,2)   = arr(:,8); %pitch 
IMU_data(:,3)   = arr(:,9); %roll
IMU_data(:,4)   = arr(:,10); %heading
IMU_data(:,5)   = arr(:,11); %fs1 
IMU_data(:,6)   = arr(:,12); %fs2
IMU_data(:,7)   = arr(:,13); %fs3
IMU_data(:,8)   = arr(:,14); %omega_s1
IMU_data(:,9)   = arr(:,15); %omega_s2
IMU_data(:,10)  = arr(:,16); %omega_s3

if Add_Bias==1
    IMU_data(:,5)   = IMU_data(:,5)+Delta_AX(1); %fs1 
    IMU_data(:,6)   = IMU_data(:,6)+Delta_AX(2); %fs2
    IMU_data(:,7)   = IMU_data(:,7)+Delta_AX(3); %fs3
    IMU_data(:,8)   = IMU_data(:,8)-Delta_DUS(1); %omega_s1
    IMU_data(:,9)   = IMU_data(:,9)-Delta_DUS(2); %omega_s2
    IMU_data(:,10)  = IMU_data(:,10)-Delta_DUS(3); %omega_s3
end

if Add_Noise==1
    IMU_data(:,5)   = IMU_data(:,5)+AxNoise; %fs1 
    IMU_data(:,6)   = IMU_data(:,6)+AyNoise; %fs2
    IMU_data(:,7)   = IMU_data(:,7)+AzNoise; %fs3
    IMU_data(:,8)   = IMU_data(:,8)-WxNoise; %omega_s1
    IMU_data(:,9)   = IMU_data(:,9)-WyNoise; %omega_s2
    IMU_data(:,10)  = IMU_data(:,10)-WzNoise; %omega_s3
end

if Add_Anomal==1
    IMU_data(:,5)   = IMU_data(:,5)-AnomalyAll(:,1); %fs1 
    IMU_data(:,6)   = IMU_data(:,6)-AnomalyAll(:,2); %fs2
    IMU_data(:,7)   = IMU_data(:,7)-AnomalyAll(:,3); %fs3
end

% -------------recording---------------------------------------------------
StartNumImu        = find(arr(:, 1)==StartTime); 
if AllGals==7
    EndNumImu          = find(arr(:, 1)==EndTime7Gals);
elseif AllGals==5 
    EndNumImu          = find(arr(:, 1)==EndTime5Gals);
elseif AllGals==3 
    EndNumImu          = find(arr(:, 1)==EndTime3Gals);
elseif AllGals==2
    EndNumImu          = find(arr(:, 1)==EndTime2Gals);    
elseif AllGals==1
    EndNumImu          = find(arr(:, 1)==EndTime1Gals);
else
    disp('Error: check "AllGals" parametr')
end

if Add_Anomal==1
    F = fopen('../Output_data/anomaly_data.txt', 'w');
    fprintf(F,'%20s %20s %20s %20s\n ','time[s]', 'dg1[m/s^2]', 'dg2[m/s^2]', 'dg3[m/s^2]');
    for i=StartNumImu:EndNumImu
        fprintf(F,'%20.15f %20.15f %20.15f %20.15f \n', AnomalyAll_x(i, :));
    end
    fclose(F);
end

    
F = fopen('../Output_data/IMU_data.txt', 'w');
fprintf(F,'%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n','Time[s]','Pitch[d]','Roll[d]','TrueHeading[d]','Fs1[m/s^2]','Fs2[m/s^2]','Fs3[m/s^2]','omega_s1[deg/h]','omega_s2[deg/h]','omega_s3[deg/h]');
for i=StartNumImu:EndNumImu
    fprintf(F,'%20.4f %20.11f %20.11f %20.11f %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f\n', IMU_data(i,:));
    if (i/10000-floor(i/10000)==0)
        disp(num2str(i))
    end
end
fclose(F);

disp(['����� ������ ���������: ', num2str(etime(clock, timetime)), ' ������'])

fid  = fopen('../Input_data/anomaly/XGM2019e_2159_2100.dat'); 
data = textscan(fid,'%f %f %f %f %f %f %f','HeaderLines', 34);
data = cell2mat(data);
fclose(fid);

Start1stGals=13989;
Finish1stGals=16130.05;
Finish2ndGals=18549.85;
% Time1stGals=find(data(:, 1)>Start1stGals & data(:, 1)<Finish1stGals);
Time2ndGals=find(data(:, 1)>=Finish1stGals & data(:, 1)<Finish2ndGals);

% have chousen 2nd gals

% plot(data(:, 1), data(:, 5));
% plot(data(Time1stGals,1), data(Time1stGals, 5));
% plot(data(Time2ndGals, 1), data(Time2ndGals, 5));

anomaly1gals    = [data(Time2ndGals, 5),data(Time2ndGals, 6),data(Time2ndGals, 7)] ;

anomaly1gals    =  anomaly1gals-anomaly1gals(1, :); %1st record eq to zero

new_time_anomal=linspace(data(Time2ndGals(1), 1), data(Time2ndGals(end), 1), 237500)';
InterpAnomal = interp1(data(Time2ndGals, 1), anomaly1gals, new_time_anomal, 'spline');

D               = importdata('../Input_data/Trajectory.txt');    
arr             = D.data;
clear D;
Length          = size(arr,1);
AnomalyAll=zeros(Length, 3);
% 
cfg;

NumGals=AllGals;

for i = 1:NumGals

    StartImu=650+i*2375;
    FinishImu=3025+i*2375;
    
    TimeImu=find(arr(:, 1)>StartImu & arr(:, 1) <= FinishImu);
    
    if i/2-floor(i/2)==0
        AnomalyAll(TimeImu, :)    = flipud(InterpAnomal);
    %     inverse array 
    else
        AnomalyAll(TimeImu, :)    = InterpAnomal;
    end
end

plot(arr(:, 1), AnomalyAll(:, 1));

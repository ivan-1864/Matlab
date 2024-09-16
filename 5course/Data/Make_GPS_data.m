%описание программы
% программа берет данные из файла с идеальными данными (время, коорд, скорости),
% обрезает и прореживает данные. далее добавляет шум,
% создает и записывает показания в новый файл.

%конфигурациионные данные
cfg

% инициализация
D = importdata('../Input_data/Trajectory.txt');    
arr = D.data;
clear D;
[Length,Width]=size(arr);
time1=0; %время начала отрезка
time2=arr(Length, 1); %время конца отрезка
delta_time=0.1; %желаемый шаг времени

% прореживание
tsin = timeseries(arr(:,2:Width),arr(:,1)); 
tsout = resample(tsin,(time1:delta_time:time2-0.1)); 
new_arr =[tsout.Time-time1,tsout.Data(:,:)];
New_len=(size(new_arr,1));


% -------------------------------------------------------------------------

%---------------добавление шумов-------------------------------------------
SKO_vel=Add_Noise*0.08; %СКО  скорости
rand1=SKO_vel*randn(1,New_len+2);
rand2=SKO_vel*randn(1,New_len+2);
rand3=SKO_vel*randn(1,New_len+2);
diffrand1=(rand1(3:end)-rand1(1:end-2))/(2*delta_time);
diffrand2=(rand2(3:end)-rand2(1:end-2))/(2*delta_time);
diffrand3=(rand3(3:end)-rand3(1:end-2))/(2*delta_time);
% -------------------------------------------------------------------------
GPS_data=zeros(New_len,7);
GPS_data(:,1)=new_arr(:,1); %время
GPS_data(:,2)=new_arr(:,2); %Lon
GPS_data(:,3)=new_arr(:,3); %Lat
GPS_data(:,4)=new_arr(:,4); %Hei
GPS_data(:,5)=new_arr(:,5)+SKO_vel*diffrand1'; %Ve
GPS_data(:,6)=new_arr(:,6)+SKO_vel*diffrand2'; %Vn
GPS_data(:,7)=new_arr(:,7)+SKO_vel*diffrand3'; %Vup
% STD = 0.045 m/s

% запись
StartNumGPS=find(new_arr(:, 1)==StartTime);
if AllGals==7
    EndNumGPS=find(new_arr(:, 1)==EndTime7Gals);
elseif AllGals==5
    EndNumGPS=find(new_arr(:, 1)==EndTime5Gals);
elseif AllGals==3
    EndNumGPS=find(new_arr(:, 1)==EndTime3Gals);
elseif AllGals==2
    EndNumGPS=find(new_arr(:, 1)==EndTime2Gals);    
elseif AllGals==1
    EndNumGPS=find(new_arr(:, 1)==EndTime1Gals);
else 
    disp('Error: check "AllGals" parametr')
end
    
F = fopen('../Output_data/GPS_data.txt', 'w');
fprintf(F,'%s\n','    Time[s]  Lon[d]       Lat[d]        Hei[m]      Ve[m/s]        Vn[m/s]     Vup[m/s]');
for i=StartNumGPS:EndNumGPS
    fprintf(F,'%10.4f %10.11f %10.11f %10.5f %10.5f %10.5f %10.5f\n', GPS_data(i,:));
end
fclose(F);


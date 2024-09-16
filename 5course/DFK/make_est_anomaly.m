D=importdata('../../Data/Output_data/anomaly_data.txt');
true_anomaly=D.data;
clear D


D=importdata('../Output_data/DFK_estimation.txt');
X=D.data;
clear D

D=importdata('../Output_data/Log.txt');
koef=D.data(1);
num_spline=D.data(2);
clear D


D=importdata('../Output_data/spline_coef.txt');
x=D;
clear D

Length=length(X);

approx_anomaly=zeros(Length, 3);
for i = 1:Length
    time_spline=(X(i, 1)-X(1,1))/koef;
    B=zeros(3, 3*num_spline);
    B1=B_spline(time_spline-floor(time_spline)+3)*eye(3);
    B2=B_spline(time_spline-floor(time_spline)+2)*eye(3);
    B3=B_spline(time_spline-floor(time_spline)+1)*eye(3);
    B4=B_spline(time_spline-floor(time_spline))*eye(3);
    B14=[B1, B2, B3, B4];
    B(:,3*floor(time_spline)+1:3*floor(time_spline)+12)=B14;
    approx_anomaly(i, :)=(B*x)';
end

est_anomaly_all=zeros(Length, 3);
est_anomaly=zeros(2375*10, 3);

  NumGals=5;
  for i = 1:NumGals-1
    
        StartImu=650+(i-1)*2375;
        FinishImu=3025+(i-1)*2375;
        
        TimeGals=find(X(:, 1)>StartImu & X(:, 1) <= FinishImu);
        
        if i/2-floor(i/2)==0
            est_anomaly    =  est_anomaly+flipud(approx_anomaly(TimeGals, :));     
        %     inverse array 
        else
            est_anomaly    = est_anomaly+approx_anomaly(TimeGals, :);
        end
   end
   est_anomaly=est_anomaly/NumGals;

for i = 1:NumGals-1
    
        StartImu=650+(i-1)*2375;
        FinishImu=3025+(i-1)*2375;
        
        TimeGals=find(X(:, 1)>StartImu & X(:, 1) <= FinishImu);
        
        if i/2-floor(i/2)==0
            est_anomaly_all(TimeGals, :)    =  flipud(est_anomaly);     
        %     inverse array 
        else
            est_anomaly_all(TimeGals, :)    = est_anomaly;
        end
end

% dg1
h=figure(4);
hold
plot(X(:, 1), 10^5*approx_anomaly(:,1));
plot(true_anomaly(:,1), 10^5*true_anomaly(:,2));
plot(X(:, 1), 10^5*est_anomaly_all(:, 1))
legend('approx \delta g_1','true \delta g_1', 'est_anomaly');
xlabel('time, [s]')
ylabel('dg1, [mGal]')
% saveas(h,'anomaly_1','fig');
% hold off

% dg2
h=figure(5);
hold
plot(X(:, 1), 10^5*approx_anomaly(:,2));
plot(true_anomaly(:,1), 10^5*true_anomaly(:,3));
plot(X(:, 1), 10^5*est_anomaly_all(:, 2))
legend('approx \delta g_2','true \delta g_2', 'est_anomaly_all');
xlabel('time, [s]')
ylabel('dg2, [mGal]')
% saveas(h,'anomaly_2','fig');
% hold off
% 
% dg3
h=figure(6);
hold
plot(X(:, 1), 10^5*approx_anomaly(:,3));
plot(true_anomaly(:,1), 10^5*true_anomaly(:,4));
plot(X(:, 1), 10^5*est_anomaly_all(:, 3))
legend('approx \delta g_3','true \delta g_3', 'est_anomaly_all');
xlabel('time, [s]')
ylabel('dg3, [mGal]')
% saveas(h,'anomaly_3','

% spl_gals=(num_spline-12)/3; % должно быть целым числом
% 
% c=zeros(3*spl_gals, 7);
% for i = 1:7
%     if floor(i/2)-i/2==0
%         c(:, i) = flipud(x(13+3*(i-1)*spl_gals:13+3*i*spl_gals-1));
%     else
%         c(:, i) = x(13+3*(i-1)*spl_gals:13+3*i*spl_gals-1);
%     end
% end
% c_est=zeros(3*spl_gals, 1);
% for i =1:3*spl_gals
%     c_est(i)=mean(c(i, :));
% end
% time_gals=time1(8001:29000);
% est_anomaly=zeros(210000, 3);
% for i = 1:210000-1
%     time_spline=(time_gals(i)-time_gals(1))/koef;
%     B=zeros(3, 3*spl_gals);
%     B1=B_spline(time_spline-floor(time_spline)+3)*eye(3);
%     B2=B_spline(time_spline-floor(time_spline)+2)*eye(3);
%     B3=B_spline(time_spline-floor(time_spline)+1)*eye(3);
%     B4=B_spline(time_spline-floor(time_spline))*eye(3);
%     B14=[B1, B2, B3, B4];
%     B(:,3*floor(time_spline)+1:3*floor(time_spline)+12)=B14;
%     est_anomaly(i, :)=(B*c_est)';
% end
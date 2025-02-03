GRS_a        =  6378.137;     % m, major semiaxis, kilometers, WGS-84 

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

D = importdata('../../Data/Output_data/GPS_data.txt');
% lon = D.data(:, 2);
lat = D.data(:, 3);
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
    if abs(time_spline-num_spline+3)>10^(-7) % костыль
        B14=[B1, B2, B3, B4];
        B(:,3*floor(time_spline)+1:3*floor(time_spline)+12)=B14;
    else
        B14=[B1, B2, B3];
        B(:,3*floor(time_spline)+1:3*floor(time_spline)+9)=B14;
    end
    approx_anomaly(i, :)=(B*x)';
end

est_anomaly_all=zeros(Length, 3);
est_anomaly=zeros(2375*10, 3);

  NumGals=3;
  for i = 1:NumGals
    
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

for i = 1:NumGals
    
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

TimeEnd = 2900;
TimeStart = 800;

TimeGals10Hz=TimeStart*10+1:TimeEnd*10;
TimeGals100Hz=TimeStart*100+1:10:TimeEnd*100;
sGals=deg2rad(lat(TimeStart*10+1:TimeEnd*10))*GRS_a;

% dg1
h=figure(4);
hold
% plot(X(TimeGals10Hz, 1), 10^5*approx_anomaly(TimeGals10Hz,1), 'LineWidth', 1.25);
plot(sGals, 10^5*true_anomaly(TimeStart*100+1:10:TimeEnd*100,2), 'LineWidth', 1.25);
plot(sGals, 10^5*est_anomaly_all(TimeStart*10+1:TimeEnd*10,1), 'LineWidth', 1.25)
legend('true \delta g_1', 'approx \delta g_1');
title('УОЛ (восточная компонента)')
xlabel('s, [km]')
ylabel('dg1, [mGal]')



% dg2
h=figure(5);
hold
% plot(X(TimeGals10Hz, 1), 10^5*approx_anomaly(TimeGals10Hz,2), 'LineWidth', 1.25);
plot(sGals, 10^5*true_anomaly(TimeStart*100+1:10:TimeEnd*100,3), 'LineWidth', 1.25);
plot(sGals, 10^5*est_anomaly_all(TimeStart*10+1:TimeEnd*10,2), 'LineWidth', 1.25)
title('УОЛ (северная компонента)')
legend('true \delta g_2','approx \delta g_2');
xlabel('s, [km]')
ylabel('dg2, [mGal]')

% dg3
h=figure(6);
hold
plot(X(:, 1), 10^5*approx_anomaly(:,3));
plot(true_anomaly(:,1), 10^5*true_anomaly(:,4));
plot(X(:, 1), 10^5*est_anomaly_all(:, 3))
legend('approx \delta g_3','true \delta g_3', 'est_anomaly_all');
xlabel('s, [km]')
ylabel('dg3, [mGal]')




figure(7);
hold
plot(sGals, 10^5*(est_anomaly_all(TimeStart*10+1:TimeEnd*10, 1)- true_anomaly(TimeStart*100+1:10:TimeEnd*100,2)), 'LineWidth', 1.25);
plot(sGals, 10^5*(est_anomaly_all(TimeStart*10+1:TimeEnd*10, 2)- true_anomaly(TimeStart*100+1:10:TimeEnd*100,3)), 'LineWidth', 1.25);
title('Ошибка определения УОЛ');
legend('dg1', 'dg2')
xlabel('s, [km]')
ylabel('|error dg|, [mGal]')


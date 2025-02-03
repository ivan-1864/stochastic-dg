D=importdata('../../Data/Output_data/anomaly_data.txt');
true_anomaly=D.data;
clear D

% D=importdata('../../Data/Output_data/delta_anomaly_data.txt');
% delta_anomaly=D.data;
% clear D


D=importdata('../Output_data/DFK_estimation_QR.txt');
X=D.data;
clear D

D=importdata('../Output_data/Log.txt');
koef=D.data(1);
num_spline=D.data(2);
clear D

Length=length(X);

D=importdata('../Output_data/spline_coef.txt');
x=D;
clear D
    % >
% time1=time(1:end-1, :);

% spl_gals=(num_spline)/7; % должно быть целым чисом
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

% x_est=[zeros(1, 12), c_est',c_est',c_est',c_est',c_est',c_est',c_est']';

%beta



h=figure(11);
plot(X(:, 1),rad2deg(X(:, 5:7))*3600, 'LineWidth', 1.25)
legend('\beta_1','\beta_2','\beta_3');
xlabel('time, [s]')
ylabel('beta, [arg.sec]')
title('Ошибки определения ориентации')
% saveas(h,'beta','fig')

 
% delta f
h=figure(12);
hold
plot(X(:, 1),10^(5)*X(:, 8:10), 'LineWidth', 1.25)
% plot(X(:, 1), -40*ones(size(:)),'--', 'LineWidth', 1.25)
% plot(X(:, 1), 30*ones(size(:)),'--', 'LineWidth', 1.25)
% legend('\Delta{f_1}','\Delta{f_2}','\Delta{f_3}', '\Delta{f_1} model', '\Delta{f_2} model');
title('Смещение нуля АКС')
xlabel('time, [s]')
ylabel('\Delta f^0, [mGal]')
% saveas(h,'delta_f','fig')

% nu
h=figure(13);
hold
plot(X(:, 1),rad2deg(X(:, 11:13))*3600, 'LineWidth', 1.25)
% plot(X(:, 1), -3*10^(-3)*ones(size(Time)),'--', 'LineWidth', 1.25)
% plot(X(:, 1), 3*10^(-3)*ones(size(Time)),'--', 'LineWidth', 1.25)
% plot(X(:, 1), 1*10^(-3)*ones(size(Time)),'--', 'LineWidth', 1.25)
title('Дрейф ДУС')
% legend('\nu_1','\nu_2','\nu_3', '\nu_1 model','\nu_2 model', '\nu_3 model');
xlabel('time, [s]')
ylabel('\nu^0, [deg/h]')
% saveas(h,'nu','fig')


% 
% approx_anomaly=zeros(Length, 3);
% est_anomaly=zeros(Length-1, 3);
% for i = 1:Length
%     time_spline=(X(i, 1)-X(1,1))/koef;
%     B=zeros(3, 3*num_spline);
%     B1=B_spline(time_spline-floor(time_spline)+3)*eye(3);
%     B2=B_spline(time_spline-floor(time_spline)+2)*eye(3);
%     B3=B_spline(time_spline-floor(time_spline)+1)*eye(3);
%     B4=B_spline(time_spline-floor(time_spline))*eye(3);
%     B14=[B1, B2, B3, B4];
%     B(:,3*floor(time_spline)+1:3*floor(time_spline)+12)=B14;
%     approx_anomaly(i, :)=(B*x)';
% 
% end
% % 
% % dg1
% h=figure(14);
% hold
% plot(X(:, 1), 10^5*approx_anomaly(:,1));
% plot(true_anomaly(:,1), 10^5*true_anomaly(:,2));
% legend('approx \delta g_1','true \delta g_1');
% xlabel('time, [s]')
% ylabel('dg1, [mGal]')
% % saveas(h,'anomaly_1','fig');
% % hold off
% 
% % dg2
% h=figure(15);
% hold
% plot(X(:, 1), 10^5*approx_anomaly(:,2));
% plot(true_anomaly(:,1), 10^5*true_anomaly(:,3));
% legend('approx \delta g_2','true \delta g_2');
% xlabel('time, [s]')
% ylabel('dg2, [mGal]')
% % saveas(h,'anomaly_2','fig');
% % hold off
% % 
% % dg3
% h=figure(16);
% hold
% plot(X(:, 1), 10^5*approx_anomaly(:,3));
% plot(true_anomaly(:,1), 10^5*true_anomaly(:,4));
% % plot(delta_anomaly(:, 1), 10^5*delta_anomaly(:, 4))
% legend('approx \delta g_3','true \delta g_3');
% xlabel('time, [s]')
% ylabel('dg3, [mGal]')
% % saveas(h,'anomaly_3','fig');
% % 
% % 
% % Fs = 100; 
% % L = length(X(:, 4:6));
% % INPUT=rad2deg(X(:, 5:7))*3600; 
% % [Pxx, freq] = pwelch(INPUT-mean(INPUT),[],0,[],Fs);
% % % h=figure(7);
% % % hold
% % figure(7)
% % plot(freq, Pxx(:, 1));
% % 
% % figure(8)
% % plot(freq, Pxx(:, 2));
% % 
% % figure(9)
% % plot(freq, Pxx(:, 3));
% 
% % xlabel('time, [s]')
% % ylabel('error dg1, [mGal]')
% % saveas(h,'err_anomaly_1','fig');
% % 
% figure(8);
% hold
% plot(X(:, 1), 10^5*(approx_anomaly(:,1)- true_anomaly(1:10:end-1,2)));
% plot(X(:, 1), 10^5*(approx_anomaly(:,2)- true_anomaly(1:10:end-1,3)));
% xlabel('time, [s]')
% ylabel('error dg, [mGal]')
% 
% 
% 

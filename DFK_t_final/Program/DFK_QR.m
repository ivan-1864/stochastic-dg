    % close all
clear

timetime = clock;

D=importdata('../../Data/Output_data/anomaly_data.txt');
true_anomaly=D.data;
clear D

% --------------download data----------------------------------------------
D               = importdata('../../INS/Output_data/output_INS_data.txt');  
time            = D.data(:, 1);
dVx_array       = D.data(:, 2:4);
Vx_array        = D.data(:, 5:7);
Omega_x_array   = D.data(:, 8:10);
u_E_array       = D.data(:, 23:25);
L_zx_array      = [D.data(:, 11), D.data(:, 12), D.data(:, 13),...
                   D.data(:, 14), D.data(:, 15), D.data(:, 16),...
                   D.data(:, 17), D.data(:, 18), D.data(:, 19)];
g0_array        = D.data(:, 20:22); 
omega_x_array   = Omega_x_array + u_E_array;
clear D



cfg; % download configuration param.
% ---------------global variables--------------------------------------
Length          = size(time, 1);
time_step       = mean(diff(time));
% -------------------------------------------------------------------------


x_val           = time;
x_min           = min(x_val);
x_max           = max(x_val);

if TimeEnd == 0
    TimeEnd = Length;
end

% ----------creating matrix P, Q, R (for KF)-------------------------------

d_dv            = Std_dV * ones(1, 3);
d_beta          = deg2rad(Std_beta / 60) * ones(1, 3);
d_dg            = Std_dg * 10^(-5) * ones(1, 3);
d_p             = Std_p * 10^(-5) * ones(1, 3);
d_f             = [Std_f * 10^(-5) * ones(1, 2), 10^(-5)];
d_nu            = deg2rad(Std_nu) / 3600 * ones(1, 3);
P_last_diag     = [d_dv, d_beta, d_f, d_nu, d_dg, d_p]; % diagonal of sqrt of cov. matrix at time t=0
P_last          = diag(P_last_diag);
% -------------------------------------------------------------------------
d_ax            = (10^(-5) * Std_AX) * ones(1, 3);
d_dus           = (deg2rad(Std_DUS) / 3600) * ones(1, 3);
d_q             = (10^(-5) * Std_q) * ones(1, 3);
Q_sqrt_diag     = [d_ax, d_dus, d_dg];
Q_sqrt          = diag(Q_sqrt_diag);
% -------------------------------------------------------------------------
d_vel           = (Std_GPS) * ones(1, 3);
R_sqrt_diag     = d_vel;
R_sqrt          = diag(R_sqrt_diag);
R_sqrt_inv      = diag(1 ./ R_sqrt_diag);

% --------------------------state vector--------------------
dimX            = 18;  %  x = (dlt V, dlt F, nu, beta, dg, p)
Y_last          = zeros(dimX, 1);
%----------------- Declare arrays for results ---------------------- 
Y_forw          = zeros(TimeEnd, dimX);
P_forw          = cell(1, TimeEnd);
P_forw{1}       = P_last;

% ------------------KF algorithm-------------------------------------------
disp('начало работы алгоритма')

for i = 1:TimeEnd
%    --------------initialization-------------------------------------------
    g_0         = Make_Cross_Matrix(g0_array(i, :));
    Vx          = Make_Cross_Matrix(Vx_array(i, :));
    Z_t         = dVx_array(i, :)';
    omega_x     = Make_Cross_Matrix(omega_x_array(i, :));
    omega_pl_u  = Make_Cross_Matrix(Omega_x_array(i, :)) + 2 * Make_Cross_Matrix(u_E_array(i, :));
%     ---------------------------------------------------------------------
    L_zx        = [L_zx_array(i,1), L_zx_array(i,2), L_zx_array(i,3);
                   L_zx_array(i,4), L_zx_array(i,5), L_zx_array(i,6);
                   L_zx_array(i,7), L_zx_array(i,8), L_zx_array(i,9)];
    
%     ---------------------------------------------------------------------
    A           = [omega_pl_u, g_0, L_zx' , Vx * L_zx', -eye(3), zeros(3);
                   zeros(3), omega_x, zeros(3), L_zx', zeros(3, 6);
                   zeros(6, 18);
                   zeros(3, 15), eye(3);
                   zeros(3, 18)];
       
    H_t         = [eye(3), -Vx, zeros(3, 12)];
    F_t         =  eye(size(A)) + A * time_step;
%     ---------------------------------------------------------------------
    J_t         = [L_zx', Vx * L_zx', zeros(3);
                  zeros(3), L_zx', zeros(3);         
                  zeros(9);
                  zeros(3, 6), eye(3);];
    
    [Y_pred,P_pred,Tmp1,Tmp2,Resid] = KF_forward( ...
        F_t,J_t,Q_sqrt,Z_t,H_t,R_sqrt,R_sqrt_inv,Y_last,P_last);
    Y_last = Y_pred; 
    P_last = P_pred;  
    
    % Copy into arrays
    Y_forw(i+1,:) = Y_pred';% 1 x Nx
    P_forw{i+1}   = P_pred;

end
%----------Estimates from Forward KF ------------------------------------
X_j=zeros(dimX, TimeEnd);

for j = 1 : TimeEnd       
    Y_current = Y_forw(j, :)';    
    P_current = P_forw{j};    
    
    P_j       = P_current * P_current';  % Covariance matrix    
    X_j(:, j) = P_current * Y_current;   % State vector estimate  
end
% ---------------------END of DFK--------------------------

% ---------------------record into the file--------------------------------

X_j = X_j';

% dv
figure(1)
hold
title('dV')
plot(X_j(:, 1:3))
% 
% beta
figure(2)
hold
title('beta')
plot(X_j(:, 4:6))
% 
% 
% df
figure(3)
hold
title('Delta f')
plot(X_j(:, 7:9))


% nu
figure(4)
hold
title('nu')
plot(X_j(:, 10:12))

time_end = 5250;
windowSize = 80*10; 
sft = windowSize/2;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
dg_sm = filter(b, a, X_j(:, 13:15));

% dg
figure(5)
hold
plot(time(1:time_end*10), 10^5*dg_sm(:, 3))
plot(true_anomaly(1+sft*10:10:time_end*100+sft*10,1),10^5*true_anomaly(1:10:time_end*100,4));
legend('comp', 'true')

figure(6)
hold
plot(time(1:time_end*10), 10^5*dg_sm(:, 1))
plot(true_anomaly(1+sft*10:10:time_end*100+sft*10,1), 10^5*true_anomaly(1:10:time_end*100,2));
legend('comp', 'true')

figure(7)
hold
plot(time(1:time_end*10), 10^5*dg_sm(:, 1))
plot(true_anomaly(1+sft*10:10:time_end*100+sft*10,1), 10^5*true_anomaly(1:10:time_end*100,3));
legend('comp', 'true')

disp(['время работы программы: ', num2str(etime(clock, timetime)), ' секунд'])
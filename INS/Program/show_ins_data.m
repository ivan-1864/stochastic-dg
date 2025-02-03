D=importdata('../Output_data/output_INS_data.txt');
time=D.data(:, 1);
dVx_array = D.data(:, 2:4);
Vx_array = D.data(:, 5:7);
omega_x_array = D.data(:, 8:10);
L_zx_array = [D.data(:, 11),D.data(:, 12),D.data(:, 13),...
    D.data(:, 14),D.data(:, 15),D.data(:, 16),...
    D.data(:, 17),D.data(:, 18),D.data(:, 19)];
g0_array=D.data(:, 20:22);
clear D

figure(1)
plot(time, Vx_array);
legend('V_{x_1}', 'V_{x_2}', 'V_{x_3}')

figure(22)
plot(time, dVx_array);
legend('\Delta V_{x_1}', '\Delta V_{x_2}', '\Delta V_{x_3}')

figure(3)
plot(time, omega_x_array);
legend('\omega_{x_1}', '\omega_{x_2}', '\omega_{x_3}')

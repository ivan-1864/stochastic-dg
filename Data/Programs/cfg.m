TrajectoryFile = '../Input_data/Trajectory.txt';
FlightLinesFile = '../Input_data/FlightLines.txt';
AnomalyFile = '../Input_data/anomaly/XGM2019e_2159_2100.dat';
OuputAnomalyFile = '../Output_data/anomaly_data.txt';
IMUFile = '../Output_data/IMU_data.txt';
GPSFile = '../Output_data/GPS_data.txt';

AllGals         = 6; %кол-во нужных галсов (1, 3, 5, 7)
Add_Noise       = 1; %добавить  шум
Add_Bias        = 1; %добавить смещения и дрейфы
Add_Anomal      = 1; %добавить аномалию

StartTime       = 0; %время начала записи в секунах
% EndTime1Gals    = 3025; % время конца записи для 1о галсового полета
% EndTime2Gals    = 5250; % время конца записи для 2х галсового полета 
% EndTime3Gals    = 7775; % время конца записи для 3х галсового полета
% EndTime4Gals    = 10300; % время конца записи для 4х галсового полета
% EndTime5Gals    = 12525;% время конца записи для 5и галсового полета 
% EndTime6Gals    = 15050;% время конца записи для 6и галсового полета 
% EndTime7Gals    = 17125;% время конца записи для 7и галсового полета

AccBias         = [30, -40, 0] * 10^(-5); % [m/s/s]
GyroDrift       = [-3, 3, 1] * 10^(-3); % [deg/h]

StdVel10Hz      = 0.05; % [m/s] 10Hz
StdGyro100Hz    = 3; %[deg/h]  100Hz
StdAcc100Hz     = 30 * 10^(-5); % [m/s/s] 100Hz

TStep100Hz           = 0.01; %timestep




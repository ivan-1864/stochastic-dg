% input: ideal data from Trajectory.txt,  configuration from cfg, real grav grav dist-ce vector from XGM2019 (up to 2100 degree)
% cut data and add IMU noise and bias, add grav dist-ce vector
% create IMU_data.txt file and record data

close all
clear

% configuration data
cfg

% ---------------ititalization--------------------------------------------
D              = importdata(TrajectoryFile);
Time100Hz      = D.data(:, 1);
Data100Hz      = D.data(:, 2:end); % frequency  = 100 Hz
clear D


fid = fopen(FlightLinesFile);
C = textscan(fid,'%s %f %c %f','HeaderLines',3); 
fclose(fid);

TimeGalsBeg = C{2};
TimeGalsEnd = C{4};
clear fid C

DeltaTime = mean(diff(Time100Hz));
Freq = round(1/DeltaTime);
TimeGals = TimeGalsEnd(1) - TimeGalsBeg(1) + (TimeGalsBeg(2) - TimeGalsEnd(1));

if Add_Anomal==1
% -------download anomaly--------------------------------------------------

    fid  = fopen(AnomalyFile); 
    data = textscan(fid, '%f %f %f %f %f %f %f','HeaderLines', 34);
    data = cell2mat(data);
    fclose(fid);

    clear fid

% --------------choose number of lines-------------------------------------
% StartInter = TimeGalsBeg(1) - (TimeGalsBeg(2) - TimeGalsEnd(1)) / 2;
StartInter = 0;
switch AllGals
    case 7
        EndInter = TimeGalsEnd(7) + (TimeGalsBeg(7) - TimeGalsEnd(6)) / 2 ;
    case 6
        EndInter = (TimeGalsEnd(6) + TimeGalsBeg(7)) / 2; 
    case 5
        EndInter = (TimeGalsEnd(5) + TimeGalsBeg(6)) / 2; 
    case 4
        EndInter = (TimeGalsEnd(4) + TimeGalsBeg(5)) / 2; 
    case 3
        EndInter = (TimeGalsEnd(3) + TimeGalsBeg(4)) / 2; 
    case 2
        EndInter = (TimeGalsEnd(2) + TimeGalsBeg(3)) / 2; 
    case 1
        EndInter = (TimeGalsEnd(1) + TimeGalsBeg(2)) / 2;
    otherwise
        disp('Error: check "AllGals" parametr')
end

clear EndTime7Gals EndTime6Gals EndTime5Gals 
clear EndTime4Gals EndTime3Gals EndTime2Gals EndTime1Gals

IndexInterval = (Time100Hz(:, 1) >= StartInter) & (Time100Hz(:, 1) < EndInter);
IndexIMU = find(IndexInterval);

clear IndexInterval

LenInter = length(IndexIMU);


AnomalyAll=zeros(LenInter, 3);
AnomalyAll_x    = zeros(LenInter,4);

% choose the interval from 2 line (from graph) 
Finish1stGals = 16130.05;
Finish2ndGals = 18549.85;
Index2ndGals=find(data(:, 1) >= Finish1stGals & data(:, 1) < Finish2ndGals);

Anomaly1Gals    = [data(Index2ndGals, 6),data(Index2ndGals, 7),data(Index2ndGals, 5)] ; % xi [arcsec], eta [arcsec], \Delta g [mGal]
Anomaly1Gals    =  Anomaly1Gals-Anomaly1Gals(1, :); % 1st record equal to zero
NumRec1Gals = TimeGals * Freq; % number of records in 1 line 

NewIndexAnomal=linspace(data(Index2ndGals(1), 1), data(Index2ndGals(end), 1), NumRec1Gals)';
InterpAnomal = interp1(data(Index2ndGals, 1), Anomaly1Gals, NewIndexAnomal, 'spline');

% expand anomaly to all lines

    for i = 1:AllGals
    
        StartImu  = TimeGalsBeg(i) - (TimeGalsBeg(i + 1) - TimeGalsEnd(i)) / 2;
        FinishImu = TimeGalsEnd(i) + (TimeGalsBeg(i + 1) - TimeGalsEnd(i)) / 2;
        
        TimeImu = find(Time100Hz(:, 1) >= StartImu & Time100Hz(:, 1) < FinishImu);
        
        if i / 2 - floor(i / 2) == 0
            AnomalyAll(TimeImu, :) = flipud(InterpAnomal); %  inverse array 
        else
            AnomalyAll(TimeImu, :) = InterpAnomal;
        end
    end


    for i = 1:LenInter
       L_zx         = Make_L_zx_Matrix(deg2rad(Data100Hz(i, 9)), deg2rad(Data100Hz(i, 8)), deg2rad(Data100Hz(i, 7)));
       AnomalyAll(i, 1)   = -deg2rad(AnomalyAll(i, 1) / 3600) * Geodesy_NormalGravity(Data100Hz(i, 2), Data100Hz(i, 3)); %arcsec->m/s/s
       AnomalyAll(i, 2)   = -deg2rad(AnomalyAll(i, 2) / 3600) * Geodesy_NormalGravity(Data100Hz(i, 2), Data100Hz(i, 3)); %arcsec->m/s/s
       AnomalyAll(i, 3)   = 10^(-5) * AnomalyAll(i, 3); %mGal->m/s/s
       AnomalyAll_x(i, :) = [Time100Hz(i, 1), AnomalyAll(i, :)];
       AnomalyAll(i, :)   = AnomalyAll(i, :) * L_zx'; %Mx->Mz
    end

end

if Add_Noise == 1

    

    Rand1 = randn(1, LenInter + 1)';
    DiffRand1 = Rand1(2 : end) - Rand1(1 : end - 1);
    DiffRand1 = DiffRand1 ./ std(DiffRand1);
    WxNoise = StdGyro100Hz * DiffRand1; 
    
    Rand2 = randn(1, LenInter + 1)';
    DiffRand2 = Rand2(2 : end) - Rand2(1 : end - 1);
    DiffRand2 = DiffRand2 ./ std(DiffRand2);
    WyNoise = StdGyro100Hz * DiffRand2; %std=1*10^(-4) rad/s = 2.1 deg/h
    
    Rand3 = randn(1, LenInter + 1)';
    DiffRand3 = Rand3(2 : end) - Rand3(1 : end - 1);
    DiffRand3 = DiffRand3 ./ std(DiffRand3);
    WzNoise = StdGyro100Hz * DiffRand3; %model x=a\dot{q}, 100 Hz

    
    AxNoise = StdAcc100Hz * randn(1, LenInter)'; %std=30 mGal, model x=aq
    AyNoise = StdAcc100Hz * randn(1, LenInter)'; %std=30 mGal, model x=aq, 100 Hz
    AzNoise = StdAcc100Hz * randn(1, LenInter)'; %std=30 mGal, model x=aq, 

    clear Rand1 Rand3 Rand3 DiffRand1 DiffRand2 DiffRand3

end

% adding the errors in the data
IMU_data = zeros(LenInter, 10);
IMU_data(:, 1)   = Time100Hz(IndexIMU, 1); %time
IMU_data(:, 2)   = Data100Hz(IndexIMU, 7); %pitch 
IMU_data(:, 3)   = Data100Hz(IndexIMU, 8); %roll
IMU_data(:, 4)   = Data100Hz(IndexIMU, 9); %heading
IMU_data(:, 5)   = Data100Hz(IndexIMU, 10); %fs1 
IMU_data(:, 6)   = Data100Hz(IndexIMU, 11); %fs2
IMU_data(:, 7)   = Data100Hz(IndexIMU, 12); %fs3
IMU_data(:, 8)   = Data100Hz(IndexIMU, 13); %omega_s1
IMU_data(:, 9)   = Data100Hz(IndexIMU, 14); %omega_s2
IMU_data(:, 10)  = Data100Hz(IndexIMU, 15); %omega_s3

if Add_Bias == 1
    IMU_data(:, 5)  = IMU_data(:, 5) + AccBias(1); %fs1 
    IMU_data(:, 6)  = IMU_data(:, 6) + AccBias(2); %fs2
    IMU_data(:, 7)  = IMU_data(:, 7) + AccBias(3); %fs3
    IMU_data(:, 8)  = IMU_data(:, 8) - GyroDrift(1); %omega_s1
    IMU_data(:, 9)  = IMU_data(:, 9) - GyroDrift(2); %omega_s2
    IMU_data(:, 10) = IMU_data(:, 10) - GyroDrift(3); %omega_s3
end

if Add_Noise == 1
    IMU_data(:, 5)  = IMU_data(:, 5) + AxNoise; %fs1 
    IMU_data(:, 6)  = IMU_data(:, 6) + AyNoise; %fs2
    IMU_data(:, 7)  = IMU_data(:, 7) + AzNoise; %fs3
    IMU_data(:, 8)  = IMU_data(:, 8) - WxNoise; %omega_s1
    IMU_data(:, 9)  = IMU_data(:, 9) - WyNoise; %omega_s2
    IMU_data(:, 10) = IMU_data(:, 10) - WzNoise; %omega_s3
end

if Add_Anomal == 1
    IMU_data(:, 5) = IMU_data(:,5) - AnomalyAll(:, 1); %fs1 
    IMU_data(:, 6) = IMU_data(:,6) - AnomalyAll(:, 2); %fs2
    IMU_data(:, 7) = IMU_data(:,7) - AnomalyAll(:, 3); %fs3
end

% -------------recording---------------------------------------------------
if Add_Anomal == 1
    F = fopen(OuputAnomalyFile, 'w');
    fprintf(F,'%20s %20s %20s %20s\n ','time[s]', 'dg1[m/s^2]', 'dg2[m/s^2]', 'dg3[m/s^2]');
    for i = 1 : LenInter
        fprintf(F,'%20.15f %20.15f %20.15f %20.15f \n', AnomalyAll_x(i, :));
    end
    fclose(F);
end

    
F = fopen(IMUFile, 'w');
fprintf(F,'%20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n','Time[s]','Pitch[d]','Roll[d]','TrueHeading[d]','Fs1[m/s^2]','Fs2[m/s^2]','Fs3[m/s^2]','omega_s1[deg/h]','omega_s2[deg/h]','omega_s3[deg/h]');
for i = 1 : LenInter
    fprintf(F,'%20.4f %20.11f %20.11f %20.11f %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f\n', IMU_data(i,:));
end
fclose(F);





% input: ideal data from Trajectory.txt,  configuration from cfg
% resamle it to 10Hz, cut data and add noise (to the velocity)
% create GPS_data.txt file and record data

close all
clear

% configuration data
cfg

%---------------ititalization--------------------------------------------
D              = importdata(TrajectoryFile);
Time100Hz      = D.data(:, 1);
Data100Hz      = D.data(:, 2:end); % frequency  = 100 Hz
clear D;

Len100Hz       = length(Data100Hz);
TimeStart      = 0; 
TimeEnd        = Data100Hz(Len100Hz, 1); 

% resempling
Time10Hz  = Time100Hz(1:10:end);
Data10Hz  = resample(Data100Hz, 1, 10);
Len10Hz   = size(Data10Hz, 1);
DeltaTime = mean(diff(Data10Hz(:, 1)));

clear Data100Hz Len100Hz 

fid = fopen('../Input_data/FlightLines.txt');
C = textscan(fid,'%s %f %c %f','HeaderLines',3); 
fclose(fid);

TimeGalsBeg = C{2};
TimeGalsEnd = C{4};
clear fid C

% --------------choose number of lines-------------------------------------

StartInter = StartTime;
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

    IndexInterval = (Time10Hz(:, 1) >= StartInter) & (Time10Hz(:, 1) < EndInter);
    IndexGPS = find(IndexInterval);

    clear IndexInterval
    
    LenInter = length(IndexGPS);

%------------add noise -------------------------------------------
if Add_Noise == 1
    Rand1 = randn(1, LenInter + 2);
    Rand2 = randn(1, LenInter + 2);
    Rand3 = randn(1, LenInter + 2);

    Diffrand1 = (Rand1(3:end) - Rand1(1:end - 2)) / (2 * DeltaTime); % central difference
    Diffrand2 = (Rand2(3:end) - Rand2(1:end - 2)) / (2 * DeltaTime);
    Diffrand3 = (Rand3(3:end) - Rand3(1:end - 2)) / (2 * DeltaTime);

    Diffrand1 = Diffrand1 ./ std(Diffrand1); % normalization
    Diffrand2 = Diffrand2 ./ std(Diffrand2);
    Diffrand3 = Diffrand3 ./ std(Diffrand3);

    clear Rand1 Rand3 Rand3
end

GPS_data = zeros(LenInter, 7);

GPS_data(:, 1) = Time10Hz(IndexGPS, 1); %time
GPS_data(:, 2) = Data10Hz(IndexGPS, 1); %Lon
GPS_data(:, 3) = Data10Hz(IndexGPS, 2); %Lat
GPS_data(:, 4) = Data10Hz(IndexGPS, 3); %Hei

if Add_Noise == 1
    GPS_data(:, 5) = Data10Hz(IndexGPS, 4) + StdVel10Hz * Diffrand1'; %Ve
    GPS_data(:, 6) = Data10Hz(IndexGPS, 5) + StdVel10Hz * Diffrand2'; %Vn
    GPS_data(:, 7) = Data10Hz(IndexGPS, 6) + StdVel10Hz * Diffrand3'; %Vup
else
    GPS_data(:, 5) = Data10Hz(IndexGPS, 4); %Ve
    GPS_data(:, 6) = Data10Hz(IndexGPS, 5); %Vn
    GPS_data(:, 7) = Data10Hz(IndexGPS, 6); %Vup
end

%--------------recording---------------------------------------------------
    
F = fopen(GPSFile, 'w');
fprintf(F,'%s\n','    Time[s]  Lon[d]       Lat[d]        Hei[m]      Ve[m/s]        Vn[m/s]     Vup[m/s]');
for i = 1:LenInter
    fprintf(F,'%10.4f %10.11f %10.11f %10.5f %10.5f %10.5f %10.5f\n', GPS_data(i, :));
end
fclose(F);

% -----------plot----------------------------------------------------------

figure(1)
hold
plot(GPS_data(:, 1), GPS_data(:, 2:3))
title('coordinates from GPS')
legend('Lon', 'Lan')
xlabel('time, s')
ylabel('ang, deg')


figure(2)
hold
plot(GPS_data(:, 1), GPS_data(:, 4))
title('height from GPS')
xlabel('time, s')
ylabel('hei, m')


figure(3)
hold
plot(GPS_data(:, 1), GPS_data(:, 5:7))
title('velocity from GPS')
legend('Ve', 'Vn', 'Vup')
xlabel('time, s')
ylabel('vel, m/s')


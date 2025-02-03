% -------------------------------------------------------------------------
num_knt = 25*2; %суммарное кол-во узлов на отрезке (включая границы)
TimeBeg = 0; % начало счета(в итерациях)
TimeEnd = 52500; %окончание счета (в итерациях) , 0 -> to the end of rec

% -----------------начальные условия---------------------------------------

% -----------для матрицы P-------------------------------------------------
Std_dV      = 0.1; %СКО ошибок определения скоростей [м/с], 10Hz 0.1
Std_beta    = 0.2; %СКО ошибок определения углов [угл. мин], 10Hz 0.1
Std_dg      = 10; %СКО оценки ВСТ
Std_p       = 10; %СКО оценки производной ВСТ
Std_f       = 20; %СКО смещения нуля АКС [мГал], 10Hz 20
Std_nu      = 0.005; %СКО дрейфов ДУС [град/час], 10Hz 0.005


% -----------для матрицы Q-------------------------------------------------
Std_AX      = 10; %СКО акселерометров [мГал], 1Hz 3 = 10 10 Hz
Std_DUS     = 0.01; %СКО ДУС [град/ч], 1Hz 0.003 = 0.01 10 Hz
Std_q       = 10; %СКО
% -----------для матрицы R-------------------------------------------------
Std_GPS=0.03; %СКО приемника ГНСС [м/с], 10Hz 0.03
% -------------------------------------------------------------------------

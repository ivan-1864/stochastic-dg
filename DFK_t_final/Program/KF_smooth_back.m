function [Y_sm,I_sm] = KF_smooth_back(Y,Y_prev,Y_sm_prev,I_sm_prev,Tmp1,Tmp2,Tmp3)
%--------------------------------------------------------------------------
% Inputs:       Y is an n x 1 vector containing output of the Forward-Time 
%               Kalam Filter at the time instant t_k  
%
%               Y_prev is an n x 1 vector containing output of the  
%               Forward-Time Kalam Filter at the time instant t_{k+1} 
%
%               Y_sm_prev is an n x 1 vector containing output of the  
%               Kalam Smoother at the time instant t_{k+1} 
%
%               I_sm_prev is an n x n matrix containing output of the  
%               Kalam Smoother at the time instant t_{k+1} 
%
%               Tmp1, Tmp2, Tmp3 are auxiliary matrices
%
% Outputs:      Y_sm is an n x 1 vector = the information vector of  
%               the state vector estimate at the time instant t_{k} 
%
%               I_sm is an n x n matrix  = the information matrix of  
%               the state vector estimate error at the time instant t_{k} 
%
% Description:  Kalman Smoothing at one iteration using Square-Root
%               Information Form (RTS algorithm),
%               see T. Kailath "Linear Estimation" (2001) for description
%
% Author:       Vadim Vyazmin, NavLab, MSU
%
% Date:         December 3, 2016
%--------------------------------------------------------------------------



%--------------------- Backward-time iteration ----------------------------
% Information vector for smoothed estimate:
Y_sm = Y + Tmp1*Tmp3  + Tmp2*(Y_sm_prev - Y_prev);

% Information matrix for the smoothed estimate error: 
I_sm = eye(size(Tmp1,1)) - Tmp1*Tmp1' - Tmp2*Tmp2' + Tmp2*I_sm_prev*Tmp2';

end
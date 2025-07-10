%% 
clc;

% フロー速度などから雑にIを推定
ne = 1.5e20;
q = 1.6e-19;
me = 9.1e-28;
k = 1.38e-23;    % [J/K]
Te = 4*q/k; % [K]
Ti = 4*q/k; 
% mi = 1.672e-27; % [kg] Hの場合
mi = (18 / 6.02e23 - me) * 1e-3; % [kg] Arの場合
ud = 10000;
E = ud/sqrt(2*k*Ti/mi);
j = ne*q*sqrt(k*Ti/(2*pi*mi))*(exp(-E^2)+E*sqrt(pi)*(erf(E)+1));

uth = sqrt(k*Ti/mi);
cs = sqrt(2*k*Te/mi);

J = ne*q*exp(-0.5)*sqrt(k*Te/mi)*exp((uth+ud)^2/(cs^2));
S = 0.001 * 0.005;
I = j * S;
II = J * S;
disp(strcat("I = ", num2str(I)))
disp(strcat("II = ", num2str(II)))

%%
% 電流比(Iratio)とフロー速度(flow)の関係
Iratio = 10;

K = 2;
M = log(Iratio)/K;
k = 1.38e-23;    % [J/K]
% T_e：電子温度→青山さんのp.35を参照 5~8[eV]
q = 1.60e-19;  % [C]
Te = 4*q/k; % [K]
mi = 1.672e-27;   % [kg]
cs = sqrt(2*k*Te/mi);
ud = cs*M;
flow = log(Iratio)/K * cs *1e-3;  % [km/s]
disp(flow)
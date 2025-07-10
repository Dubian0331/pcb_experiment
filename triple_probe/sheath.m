% シースの厚さ:Dを計算
% https://www.jstage.jst.go.jp/article/jspf1958/68/6/68_6_550/_pdf
% https://www.jstage.jst.go.jp/article/jvsj2/59/7/59_16-LC-007/_pdf


m = 9.109*10^(-31); % e-の質量[kg]
M_H = 1.673*10^(-27); % H+の質量[kg]
M_Ar = 6.63352*10^(-26); % Ar+の質量[kg]
T_e = 5; % [eV]
n_e = 10^19; % [m^-3]
lambda_D = 0.74*(T_e/n_e)^(1/2); % デバイ長[m]

% D_H = (sqrt(2)/3)*exp(1/4)*(M_H/(2.3*m))^(4/3)*lambda_D % [m]
% D_Ar = (sqrt(2)/3)*exp(1/4)*(M_Ar/(2.3*m))^(4/3)*lambda_D % [m]

D = 50*lambda_D % [m]
% D_H_um = D_H / 10^(-6)
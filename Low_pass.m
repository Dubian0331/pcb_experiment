clear all

% Low-pass filter
% f_c = 1/(2*pi*R*C)

R = 100;
C = 1*10^(-9);
% f_c = 81*10^6;

% R = 1/(2*pi*C*f_c)
% C = 1/(2*pi*R*f_c)
f_c = 1/(2*pi*R*C)

T = 1/f_c

if T < 1*10^(-6)
    disp('T is under 1 [us]')
else
    disp('Check your calc')
end

G = tf([1], [R*C, 1]);
T = feedback(G, 1);

bode(G), grid
clc;clear;close all
%load('D:\Assignment_Data_SC42145.mat')
load('E:\TU DELFT\Q2\ROBUST\PART1\Assignment_Data_SC42145.mat')

%%
% open loop Bode plot AND pole-zero map
SS=ss(A,B,C,D);
TFs=tf(SS);
G1=TFs(1,1)

 figure()
 pzplot(G1)
% title('Pole-zero map SISO system')
% opts = bodeoptions('cstprefs');opts.FreqUnits = 'Hz';

figure()
opt=bodeoptions;
opt.PhaseWrapping='on';
bode(G1,opt); grid on;
title('Bode diagram SISO OL system (with phase wrapping)')

% figure()
 sisoolstepinfo=stepinfo(G1)
% step(G1)
% title('Step response SISO system')

figure()
%rlocus(G1)
nichols(G1)
grid on;
%% damping ratio, sensitivity, PM
zeta=(0.203-0.189)/(0.209-0.189);
ts=168.7029;
bdf=bandwidth(G1);
bd=4.05/ts;
PM=100*zeta;

%% filter omega1=0.2r/s omega2=3.29r/s
s=tf('s');
omega_1=0.2;
omega_2=3.29;
zeta_n1=0.1*(db2mag(-24.2)-db2mag(-15.3))/db2mag(-15.3);
zeta_d1=0.1;
zeta_n2=(db2mag(-43)-db2mag(-33))/db2mag(-33);
zeta_d2=1;
%cn1=(s^2+2*zeta_d1*omega_1*s+omega_1^2)/(s^2+2*zeta_n1*omega_1*s+omega_1^2);
%cn2=(s^2+2*zeta_d2*omega_2*s+omega_2^2)/(s^2+2*zeta_n2*omega_2*s+omega_2^2);
cn1=tf([25 0.39 1],[25 0.12 1]);
cn2=tf([0.09 0.037 1],[0.09 0.015 1]);
figure()
opt=bodeoptions;
opt.PhaseWrapping='on';
Ln=cn1*cn2*G1;
bode(Ln,opt); grid on;
hold on
bode(cn1);
bode(cn2);
title('Bode diagram SISO w/filter system (w/phase wrapping)')
figure()
Ciso=  -12.207*(s+1)^2*(s^2 + 13.28*s + 64)^2/(s*(s+100)^2);
stepinfo(Ciso*Ln/(1+Ciso*Ln))
step(Ciso*Ln/(1+Ciso*Ln));
title('Step response SISO CL')
figure()
bode(Ciso*Ln/(1+Ciso*Ln));
title('Bode diagram SISO CL')
margin(Ciso*Ln/(1+Ciso*Ln))
grid on;
figure()
bode(1/(1+Ciso*Ln));
title('Bode diagram SISO CL S')
figure()
bode(1-1/(1+Ciso*Ln));
title('Bode diagram SISO CL T')



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
% figure()
% opt=bodeoptions;
% opt.PhaseWrapping='on';
% Ln=cn1*cn2*G1;
% bode(G1,opt)
% bode(Ln,opt); grid on;
% hold on
% bode(cn1);
% bode(cn2);
% title('Bode diagram SISO w/filter system (w/phase wrapping)')
figure()
Ciso=-7031.2*(s+1)^2*(s^2 + 13.28*s + 64)/(s*(s+100)^2*(s+10));
stepinfo(Ciso*Ln/(1+Ciso*Ln))
step(Ciso*Ln/(1+Ciso*Ln));
title('Step response SISO CL')
figure()
bode(Ciso*Ln/(1+Ciso*Ln));
title('Bode diagram SISO CL')
margin(Ciso*Ln/(1+Ciso*Ln))
bandwidth(Ciso*Ln/(1+Ciso*Ln))
grid on;
figure()
bode(1/(1+Ciso*Ln));
title('Bode diagram SISO CL S')
figure()
bode(1-1/(1+Ciso*Ln));
title('Bode diagram SISO CL T')

minreal(Ciso*Ln/(1+Ciso*Ln))
%% 1.4 disturbance
Gd=TFs(1,3);

%% 2.1
SS=ss(A,B,C,D);
TFs=tf(SS);
G11=TFs(1,1)
G12=TFs(1,2)
G21=TFs(2,1)
G22=TFs(2,2)
G=[G11 G12;G21 G22]
w1=0;w2=0.8*pi;
format short
val1 = evalfr(G,1j*w1)
val2 = evalfr(G,1j*w2)
hw1=val1.*transpose(inv(val1))
hw2=val2.*transpose(inv(val2))

%% 2.2
G=minreal(G);
P22=pole(G)
Z22=tzero(G)

%% 2.7
%use robust control toolbox (code from Lec 2 Page 41/26)
s=tf('s');
systemnames ='G Wp Wu'
Wu=[0.01 0;0 (5*10^-3*s^2+7*10^-4*s+5*10^-5)/(s^2+14*10^-4*s+10^-6)];
Wp11=(s/1.8+0.8*pi)/(s+8*10^-5*pi);
Wp=[Wp11 0;0 0.2];
Wt=[];
[K,CL,GAM,INFO]=mixsyn(G,Wp,Wu,Wt);
P=[Wp -Wp*G;zeros(2) Wu;eye(2) -G]
inputvar ='[w(2);u(2)]';
input_to_G='[u]';
input_to_Wu='[u]';
input_to_Wp='[w-G]';
outputvar ='[Wp;Wu;w-G]';
sysoutname='P';
sysic;
P=minreal(P)
[K2,CL2,GAM2,INFO2]=hinfsyn(P,2,2);
K2=minreal(K2);
P=minreal(P);
size(P)
size(G)
size(K2)
%L22=K2*P;
%% 2.7  2.8 plot
%nyquist(K2*G) %(not for MIMO)
L22=tf(K2)*G;
figure()
det11=L22(1,1)+1;
det12=L22(1,2);
det21=L22(2,1);
det22=L22(2,2)+1;
L2=det11*det22-det12*det21;
nyquist(L2)

G13=TFs(1,3)
G23=TFs(2,3)
Gd=[G13;G23];

%% 3 d
load('E:\TU DELFT\Q2\ROBUST\robust\td.mat')
load('E:\TU DELFT\Q2\ROBUST\robust\tvari.mat')
d=1;
%tscol = tscollection(Wind_Data);
%tscol=addts(tscol,Wind_Data,'Fit1')
%plot(tscol.Fit1,'Color','r')
tsine=2*sin(pi/500*td);
d3=tvari-tsine;
figure()
plot(td,tvari);
hold on
plot(td,tsine);
figure()
plot(td,d3);

%% 3.2 3.3 3.4 3.5
s=tf('s');
SS=ss(A,B,C,D);
TFs=tf(SS);
G11=TFs(1,1);
G12=TFs(1,2);
G21=TFs(2,1);
G22=TFs(2,2);
G=[G11 G12;G21 G22];
G13=TFs(1,3);
G23=TFs(2,3);
Gd=[G13;G23];

omegalpf=2*pi/1000;
omegahpf=10;
LPF=omegalpf/(s+omegalpf);
HPF=s/(s+omegahpf);
Wu3=[LPF 0;0 HPF];
Wp3=(s/2+0.8*pi)/(s+8*10^-5*pi);
Wt=[];

P3=[Wp*Gd Wp*G;[1;0] Wu;-Gd -G];
inputvar ='[d(1);u(2)]';
input_to_G='[u]';
input_to_Wu='[u]';
input_to_Wp='[Gd-G]';
outputvar ='[Wp;Wu;Gd-G]';
sysoutname='P3';
sysic;
P3=minreal(P3);
[K3,CL3,GAM3,INFO3]=hinfsyn(P3,2,2);
K3=minreal(K3);
P3=minreal(P3);
size(P3)
size(G)
size(K3)
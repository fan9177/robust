clc;clear
load('D:\Assignment_Data_SC42145.mat')
%%
% open loop Bode plot AND pole-zero map
SS=ss(A,B,C,D);
TFs=tf(SS);
[num_g11,den_g11]=tfdata(TFs(1,1),'v');
g11=tf(num_g11,den_g11);

g_siso=g11;
%pzmap(g_siso);
%opt=bodeoptions;
%opt.PhaseWrapping='on';
%bode(g_siso,opt); grid on;

%% PID
s=tf('s');
% ku=5.25;
% tu=0;
% kp=0.7*ku;
% ki=2.5/tu;
% kd=tu/6.7;
kp=-0.0001;
% cpid=kp*(1+kd*s+ki/s);
figure(5)
% step(G1*cpid/(1+G1*cpid));
step(g_siso*kp/(1+g_siso*kp));
%turnedsisoolstepinfo=stepinfo(G1);

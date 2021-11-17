clc;clear
load('D:\BaiduNetdiskWorkspace\代尔夫特\课件\Q2  15\SC42145 Robust Control\Assignments\Assignment_Data_SC42145.mat')
%%
% open loop Bode plot AND pole-zero map
SS=ss(A,B,C,D);
TFs=tf(SS);
[num_g11,den_g11]=tfdata(TFs(1,1),'v');
g11=tf(num_g11,den_g11);

g_siso=g11;
pzmap(g_siso);
opt=bodeoptions;
opt.PhaseWrapping='on';
bode(g_siso,opt); grid on;


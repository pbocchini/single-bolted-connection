clc;
clear all;
close all;

%Input Variables 
Pslip=6.4512*0.153*2;
width_brace=2;
thickness_brace=5/16;
width_leg=4;
thickness_leg=5/16;
Fy_brace=36;
Fu_brace=60;
Fy_leg=50;
Bolt_diameter=0.625;
E=29000;
v=0.26;
Hole_diameter=0.6875;
e1=1;
%%%%%%%%%%%%%%%%%%Initiate the model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,U,K,Ftrial,Utrial,Ktrial,Keb, Leff, Kci, Kti, Rc, Rt, REPc, REPt,u_s0c, u_b0c,u_s0t,u_b0t,f_s0c, f_b0c,f_s0t,f_b0t,...
    u_bc, u_bt, u_sc, u_st,f_bc, f_bt, f_sc, f_st,ubkp,fbkp,kbkp,Kslip,u_inc,status,dir,u_current,f_current,u_transition, f_transition,clearance]...
    =SingBoltAngleJoint(Pslip,width_brace, thickness_brace, thickness_leg, Fy_brace,Fu_brace, Fy_leg, Bolt_diameter, E, v, Hole_diameter, e1);

%% Start the Loop
%Upeak=[0 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.13 -0.13 0.18 -0.18 0.3 -0.18 0.4 -0.18 ];  %H1_4
%Upeak=[0 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.13 -0.13 0.18 -0.18 0.3 -0.3 0.5 -0.3 ];   %H2_4
%Upeak=[0 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.13 -0.13 0.18 -0.18 0.23 -0.23 0.28 -0.28 ];  %H3_4
%H4_4 H25_1 H25_2 H25_3 H33_1 H33_2 H33_3
Upeak=[0 0.04 -0.04 0.06 -0.06 0.08 -0.08 0.13 -0.13 0.18 -0.18 0.23 -0.23 0.28 0 ];
%Upeak=[0  0.04 -0.04 0.06 -0.06 0.08 -0.08 -0.3 -0.1 -0.4];  


Urecord = history2(1000,Upeak);
%%%%%%%%%%%%%%%Calculate the break points based on current point%%%%%%%%%%
for i=1:size(Urecord,2)
[ubkp, fbkp, kbkp]=CalcBreakpoints(dir, status, Pslip,  u_bc, f_bc, u_bt, f_bt, u_b0c,f_b0c,u_b0t,f_b0t,u_s0c, f_s0c, u_s0t, f_s0t, u_current,f_current, Keb, Kti, Kci, REPc, REPt, Rc, Rt,Kslip, u_transition, f_transition);
%%%%%%%%%%%%%%%%%%%%%%%%%input the new displacement%%%%%%%%%%%%%%%%%
u_trial=Urecord(i);
%%%%%%%Based on break points compute the trial stress and trialstatus and
%%%%%%%trial dir%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f_trial,status_trial,dir_trial,u_transition_trial,f_transition_trial,u_bc_trial,f_bc_trial,u_bt_trial,f_bt_trial,u_b0t_trial,f_b0t_trial,u_b0c_trial,f_b0c_trial]...
    =CalcForceStatusDir(ubkp, fbkp, kbkp,u_trial,status,dir,u_current,f_current,u_st, f_st, u_sc, f_sc, Kslip, Kti, Kci,u_transition,f_transition,...
    u_bc,f_bc,u_bt,f_bt,REPt,REPc,Rt,Rc,u_b0t,f_b0t,u_b0c,f_b0c,Keb,Pslip);
%Commit the updates 
u_current=u_trial;
f_current=f_trial;
status=status_trial;
dir=dir_trial;
u_transition=u_transition_trial;
f_transition=f_transition_trial;
u_bc=u_bc_trial;
f_bc=f_bc_trial;
u_bt=u_bt_trial;
f_bt=f_bt_trial;
u_b0t=u_b0t_trial;
f_b0t=f_b0t_trial;
u_b0c=u_b0c_trial;
f_b0c=f_b0c_trial;
%Record Force
Frecord(i)=f_trial;


end


figure;
plot(Urecord(1:i)*2.54,Frecord*4.4482216);
hold on
%plot([ u_b0c u_s0c 0 u_s0t u_b0t],[f_b0c f_s0c 0 f_s0t f_b0t]);
load('H4_4.mat','U','F');
plot(U*2.54,F*4.4482216,'--','LineWidth',1,'Color','r');

xlabel('displacement (cm)','FontSize', 12);
ylabel('force (kN)','FontSize', 12);
lgd=legend('Zero-length element model','Brick elements model','Location','southeast');
lgd.FontSize=12;
title('L 2\times2\times5/16','FontSize', 12)
% xlim([-0.18 0.18])
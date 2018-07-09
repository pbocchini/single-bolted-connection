function [F,U,K,Ftrial,Utrial,Ktrial,Keb, Leff, Kci, Kti, Rc, Rt, REPc, REPt,u_s0c, u_b0c,u_s0t,u_b0t,f_s0c, f_b0c,f_s0t,f_b0t,...
    u_bc, u_bt, u_sc, u_st,f_bc, f_bt, f_sc, f_st,ubkp,fbkp,kbkp,Kslip,u_inc,status,dir,u_current,f_current,u_transition, f_transition,clearance]...
    =SingBoltAngleJoint(Pslip,width_brace, thickness_brace, ...
    thickness_leg, Fy_brace,Fu_brace, Fy_leg, Bolt_diameter, E, v, Hole_diameter, e1)

F=0;
U=0;
K=0;
Ftrial=0;
Utrial=0;
Ktrial=0;

%%%%%%%%%%%%%%%%%% Fixed Properties%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearance=Hole_diameter-Bolt_diameter;
G=E/(2*(1+v));
%Effective length calulated using euro code 3 with 30 degrees force angle
Leff=(width_brace+(width_brace-thickness_brace)/2)*3^0.5;
%Elastic stiffness of the brace plate using Leff
An=width_brace*thickness_brace+(width_brace-thickness_brace)*thickness_brace;
Aeff=An/2;
Af=width_brace*thickness_brace;
Anet=(width_brace-Bolt_diameter)*thickness_brace;
Keb=Aeff*E/Leff;
%Initial Compressive stiffness
Kbr=120*Fy_brace*thickness_brace*Bolt_diameter;  %braceplate bearing stiffness
Knet=Keb;
Abolt=pi*(Bolt_diameter/2)^2;
Lbolt=thickness_brace+thickness_leg;
Ibolt=pi/4*(Bolt_diameter/2)^4;
Kboltshear=Abolt*G/Lbolt*32/37;
Kboltbending=3*E*Ibolt/(Lbolt^3);	
KbrLeg=120*Fy_leg*thickness_leg*Bolt_diameter;
Kci=1/(1/Kbr+1/Knet+1/Kboltshear+1/Kboltbending+1/KbrLeg);
%Initial Tensile stiffness
Kbben=32*E*thickness_brace*(e1/Bolt_diameter-1/2)^3;%brace plate bending stiffness
Kbv=6.67*G*thickness_brace*(e1/Bolt_diameter-1/2);%brace plate shear stiffness
Kti=1/(1/Kci+1/Kbben+1/Kbv);
%Compressive Capacity
Rc=min(2.4*Bolt_diameter*thickness_brace*Fu_brace,Af*Fy_brace);
%Tensile Capacity
Rt=min(Anet*Fu_brace,0.7*Fu_brace*2*thickness_brace*(e1-Hole_diameter/2)/cos(pi/6));
%Richard Equation Parameter:
REPc=[7.28899865552825   %Compression
-0.00711376636993643
2.77499874147550
0.330198149367655
]; 
REPt=[4.56831690798003 %Tension
0.0137414954357986
1.04617240539503
0.493239575034980
];   

%Assign very small stiffness to slippage to avoid numerical problems 
Kslip=Keb/10000;
%u increment used for tangent evaluation 
u_inc=10^-7;


%%%%%%%%%%%%%%%%%%%%%% Variables subject to change%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Compressive Bearing Strain
u_s0c=-Pslip/Keb;                     %compression slipping displacement 
f_s0c=-Pslip;
u_b0c=-clearance+u_s0c;               % start of Richard Equation   
f_b0c=-clearance*Kslip+f_s0c;
%Initial Tensile Bearing Strain
u_s0t=Pslip/Keb;                      %tension slipping displacement 
f_s0t=Pslip;
u_b0t=clearance+u_s0t;                % start of Richard Equation  
f_b0t=clearance*Kslip+f_s0t;



%internal Variables subjected to change 
u_bc=u_b0c;                           % start of bearing 
f_bc=f_b0c; 
u_bt=u_b0t;                           % start of bearing 
f_bt=f_b0t;
u_sc= u_s0c;
f_sc= f_s0c;
u_st=u_s0t;
f_st=f_s0t;

%Initiate the breakpoints array
ubkp=[u_b0c; u_b0c; u_s0c; u_s0t; u_b0t; u_b0t];
fbkp=[f_b0c; f_b0c; f_s0c; f_s0t; f_b0t; f_b0t];
%for u_current<=ubkp[1] and u_current>=ubkp[6] we use Richard Equation so
%there is no linear stiffness
kbkp=[Kci; Kslip; Keb; Kslip; Kti];

%status 1 tension bearing nonstick -1 compression bearing nonstick 
%status 4 tension bearing stick -4 compressiong bearing stick
%status 0 before bearing 
status=0; 
%direction 1 loading -1 reversal loading 0 for start
dir=0;

u_current=0;
f_current=0;
u_transition=0;
f_transition=0;










end
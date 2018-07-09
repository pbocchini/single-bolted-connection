clc;
clear all;
close all;

%%
%establish values of ui
n=1000;

%Case 1  Applied force less than critical value
%Upeak=[0 0.001 -0.001];

%Case 2 Slipping before bearing from -0.068 to 0.068
%Upeak=[0 0.01 -0.01  0.02 -0.02 0.03 -0.03 0.04 -0.04 0.05 -0.05 0.06 -0.06 0]

%Case 3 Elastic Bearing from -0.08 to 0.08
%Upeak=[0  0.04 -0.04 0.05 -0.05 0.06 -0.06 0.07 -0.07 0.075 -0.075  0.078 -0.078 0 ];

%Case 4 Plastic bearing part 1   from -0.11 to 0.14
%Upeak=[0  0.04 -0.04 0.075 -0.075  0.078 -0.078 0.08 -0.08 0.1 -0.08 0.1 -0.1 0.1 ];

%Case 5 Positive Plastic bearing part 2   from 0.14  to 0.3
%Upeak=[0  0.04 -0.04 0.075 -0.075  0.078 -0.078 0.08 -0.08 0.1 -0.08 0.1 -0.1 0.14 0 0.15 0 0.2 0 0.25 0 0.23 0.2 0.24 0.2 0.25 0.2 0.3 0];

%Case 6 buckling with negligible effect on stiffness  from -0.15999 to 0.3
%Upeak=[0  -0.1 0 -0.11 0 -0.12 0 -0.13 0 -0.14 0 -0.15 0 -0.155 0 -0.157 0 -0.15999 0 ];

%Case 7 buckling with effect on stiffness from -0.3 to 0.3
%Upeak=[0 -0.06 0 -0.08 0  -0.1 0 -0.12 0 -0.15 0 -0.18 0 -0.21 0 -0.22 0 -0.23 0 -0.24 0 -0.27 0 -0.3 0];

%Load History 5
Upeak=[0 0.03 -0.03 0.06 -0.06 0.09 -0.09 0.12 -0.12  0.15 -0.15 0.18 -0.18 0.21 -0.21  0.24 -0.24 0.27 -0.27  0.3 0   ];
%   
%if position=0.8 0.126*(0.8-0.5)=0.0378 so it moves to left by 0.0378
%Upeak=[0 0.03 -0.03 0.06 -0.06 0.09 -0.09 0.12 -0.12 0.15 -0.15 0.18 -0.18 0.21 -0.21  0.24 -0.24 0.27 -0.27-0.0378   0.3-0.0378  0];
%   
%Load History 6
%Upeak=[0 0.06 0 0.08 0 0.1 0 0.12 0 0.15 0 0.18 0 0.21 0 0.24 0 0.27 0 0.3 0];
%   
%Load History 12
%Upeak=[0 -0.06 0 -0.08 0  -0.1 0 -0.12 0 -0.15 0 -0.18 0 -0.21 0 -0.24 0 -0.27 0 -0.3 0];

Urecord = history2(n,Upeak);


figure;
plot((0:size(Urecord,2)-1),Urecord);


%Urecord=[0.2;0];
%%

%The simulation results from abaqus
load('loadhistory3.mat');
Uhis3=U;
Fhis3=F;
load('loadhistory4.mat');
Uhis4=U;
Fhis4=F;
figure 
plot(Uhis3,Fhis3,'r');
hold on
plot(Uhis4,Fhis4,'r');
hist3=horzcat(Uhis3,Fhis3);
hist4=horzcat(Uhis4,Fhis4);

%Figure out the input parameters
u(1)=0.005; 
f(1)=2;
u(2)=0.068;
f(2)=2.02;
u(3)=0.08;
f(3)=5;
u(4)=0.14;
f(4)=7.5;
u(5)=0.3;
kpend=3;
f(5)=f(4)+(u(5)-u(4))*kpend;   %figure out the kpend
u(6)=0.25;
u(7)=0.17;
u(8)=0.27;


un(1)=-u(1); 
fn(1)=-f(1);
un(2)=-u(2);
fn(2)=-f(2);
un(3)=-0.08;
fn(3)=-5;
un(4)=-0.11;
fn(4)=-8;
un(5)=-0.16;
fn(5)=-6.0;
un(6)=-0.27;
knend=-10;
fn(6)=fn(5)+(un(6)-un(5))*knend;   %figure out the knend
un(7)=-0.215;
un(8)=-0.2;
un(9)=-0.14;

position=0.5;

%% input 
% k1=400;
% criticalfr=2;
% initial_clearance=0.125;
% k2=0.0001;
% position=0.9;  %in the middle
% [u(1),f(1)]=Intersect(k1,0,0, k2,0,criticalfr);
% [un(1),fn(1)]=Intersect(k1,0,0, k2,0,-criticalfr);
% u(2)=u(1)+initial_clearance*(1-position);
% f(2)=f(1)+k2*(u(2)-u(1));
% un(2)=un(1)-initial_clearance*position;
% fn(2)=fn(1)+k2*(un(2)-un(1));
%%
%Initialize the spring

[Fsbs,Usbs,Ksbs,Fsbst,Usbst,Ksbst,envlpPosForce,envlpNegForce,envlpPosDeform,envlpNegDeform,envlpTangent,Clearance,Umin,FUmin,Status] ...
   = SbsInit(u(1),f(1),u(2),f(2),u(3),f(3),u(4),f(4),u(5),u(6),u(7),u(8),kpend,...
   un(1),fn(1),un(2),fn(2),un(3),fn(3),un(4),fn(4),un(5),fn(5),un(6),un(7),un(8),un(9),knend,position);

envlpPosForceOriginal=envlpPosForce;
envlpNegForceOriginal=envlpNegForce;
envlpPosDeformOriginal=envlpPosDeform;
envlpNegDeformOriginal=envlpNegDeform;

%Envlope of the spring
array=[1 2 4 2 3 5 3 7 6 4 1 -5 -1 4 -1 -2 -5 -6 -3 -2 -3 -4 -8 -9 -7 -4 -7 -6 -5 -1 -11 -10 -11 6 7 9 8 1];

[ x,y ] = DrawEnvlp( array,envlpPosDeform,envlpPosForce,envlpNegDeform,envlpNegForce);
figure;
line (x,y);

load('loadhistory5.mat');
Uhis5=U;
Fhis5=F;
hold on
plot(Uhis5,Fhis5,'r');
hold on 
plot(Uhis3,Fhis3,'g');
hold on
plot(Uhis4,Fhis4,'g');

%%

%Compute F
for i=1:size(Urecord,2)
    ureocrd(i)=Urecord(i);  %for check purpose
    %[ Fsbst,Ksbst,envlpPosDeformt,envlpPosForcet,envlpNegDeformt,envlpNegForcet]  = SbsResp(Urecord(i),Fsbs,Usbs,Ksbs,envlpPosForce,envlpNegForce,envlpPosDeform,envlpNegDeform,envlpTangent,...
    %   Clearance,Umax,Umin,Fmax,Fmin);
    if Status==0  %no buckling 
    [ Path ] = CreatePath( Usbs,Fsbs,envlpPosForce,envlpNegForce,envlpPosDeform,envlpNegDeform,envlpTangent );
    [ Fsbst,Ksbst ]=Path2FsbstKsbst( Urecord(i),Path );
    [ envlpPosDeformtrial,envlpPosForcetrial,envlpNegDeformtrial,envlpNegForcetrial,envlpTangentrial,Clearancetrial,Statustrial,Umintrial,FUmintrial ] =...
        UpdateEnvlp_K_Clearance( Urecord(i),Fsbst,envlpPosForce,envlpNegForce,envlpPosDeform,envlpNegDeform,envlpTangent,Clearance,Status,Umin,FUmin );
    else %buckling 
    %1 compute the load and unload stiffness considering degradation with respect to minimal deformation
    [ kNegUnload,kNegLoad,kPosUnload,kPosLoad ] = BucklingStiffness( Umin,FUmin,envlpTangent,envlpNegDeform,envlpNegForce );
     %2 compute Fsbst Ksbst and update envelope, status, clearance,umin,and fumin
    [ Fsbst,Ksbst,envlpNegDeformtrial,envlpNegForcetrial,envlpPosDeformtrial,envlpPosForcetrial,Statustrial,Clearancetrial,Umintrial,FUmintrial,Path]...
     =ForceandUpdate_Buckling( Urecord(i),Usbs,Fsbs,kNegUnload,kNegLoad,kPosUnload,kPosLoad,Clearance,Status...
    ,envlpTangent,envlpNegDeform,envlpNegForce,envlpPosDeform,envlpPosForce,Umin,FUmin );         
    end       
    %Update  
    Fsbs=Fsbst;
    Usbs=Urecord(i);
    Ksbs=Ksbst;
    
    envlpPosDeform=envlpPosDeformtrial;
    envlpPosForce=envlpPosForcetrial;
    envlpNegDeform=envlpNegDeformtrial;
    envlpNegForce=envlpNegForcetrial;
    envlpTangent=envlpTangentrial;
    Clearance=Clearancetrial;
    Status=Statustrial;
    Umin=Umintrial;
    FUmin=FUmintrial;
    %Record Force
    Frecord(i)=Fsbs;
    Crecord(i)=Clearance;
    
 
    
end



figure;
% subplot(2,1,1);
% plot((0:size(Urecord,2)-1),Urecord);
% ylabel('displacement (in)');
% xlabel('timestep');
% subplot(2,1,2);
plot(Urecord,Frecord);
xlabel('displacement (in)');
ylabel('kips');
% 
load('loadhistory5.mat');
Uhis=U;
Fhis=F;
hold on
plot(Uhis,Fhis,'r');


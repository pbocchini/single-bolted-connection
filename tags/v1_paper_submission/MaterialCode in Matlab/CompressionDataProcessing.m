clc;
clear all;
close all;

%% initial point

y0=-6.4512*0.153*2;
x0=-0.0625+y0/2055.28*2;

%%

% 
% 
% load('C1_1.mat');
% u_hist1=U;
% f_hist1=F;
% 
% load('C1_2.mat');
% u_hist2=U;
% f_hist2=F;
% 
% load('C1_3.mat');
% u_hist3=U;
% f_hist3=F;
% 
% load('C1_4.mat');
% u_hist4=U;
% f_hist4=F;
% 
% load('C1_5.mat');
% u_hist5=U;
% f_hist5=F;
% 
% load('C1_4H.mat');
% u_hist6=U;
% f_hist6=F;
% 
% figure;
% plot(u_hist1,f_hist1);
% hold on
% plot(u_hist2,f_hist2);
% hold on
% plot(u_hist3,f_hist3);
% hold on
% plot(u_hist4,f_hist4);
% hold on
% plot(u_hist5,f_hist5);
% hold on
% plot(u_hist6,f_hist6);
% legend('C1_1','C1_2','C1_3','C1_4','C1_5');



%% Plot the original deformation vs force
load('C1_4.mat');
u_hist14=U;
f_hist14=F;
load('C2_4.mat');
u_hist24=U;
f_hist24=F;
load('C3_4.mat');
u_hist34=U;
f_hist34=F;
load('C4_4.mat');
u_hist44=U;
f_hist44=F;

load('C25_1.mat');
u_hist251=U;
f_hist251=F;
load('C25_2.mat');
u_hist252=U;
f_hist252=F;
load('C25_3.mat');
u_hist253=U;
f_hist253=F;

load('C33_1.mat');
u_hist331=U;
f_hist331=F;
load('C33_2.mat');
u_hist332=U;
f_hist332=F;
load('C33_3.mat');
u_hist333=U;
f_hist333=F;

figure;
scatter(u_hist14,f_hist14);
hold on
scatter(u_hist24,f_hist24);
scatter(u_hist34,f_hist34);
scatter(u_hist44,f_hist44);
scatter(u_hist251,f_hist251);
scatter(u_hist252,f_hist252);
scatter(u_hist253,f_hist253);
scatter(u_hist331,f_hist331);
scatter(u_hist332,f_hist332);
scatter(u_hist333,f_hist333);
title('origianl deformation vs force ');
legend('L2*1/8','L2*3/16','L2*1/4','L2*5/16','L2.5*3/16','L2.5*1/4','L2.5*5/16','L3*3/16','L3*1/4','L3*5/16');
%% Initial Stiffness and ultimate strength 
%from C14 C24 C34 C44 C331 C332 C333 C251 C252 C253
Ki=[208.2151499
284.4874084
328.3832774
354.132422
293.8252426
338.0600826
363.3630234
289.9089857
332.7073127
359.5201947
];
Rn=[
9
13.5
18
22.5
16.875
22.5
28.125
16.875
22.5
28.125
]+y0;
RnAbaqus=[-8.266
-13.95
-20.23
-25.11
-17.1
-23.91
-28.28
-15.69
-21.94
-26.75
]*-1+y0;


%% Exclude the part before bearing 
Rn_afterbearing=Rn+y0;

u14=-(u_hist14-x0);
u24=-(u_hist24-x0);
u34=-(u_hist34-x0);
u44=-(u_hist44-x0);
u331=-(u_hist331-x0);
u332=-(u_hist332-x0);
u333=-(u_hist333-x0);
u251=-(u_hist251-x0);
u252=-(u_hist252-x0);
u253=-(u_hist253-x0);
f14=-(f_hist14-y0);
f24=-(f_hist24-y0);
f34=-(f_hist34-y0);
f44=-(f_hist44-y0);
f331=-(f_hist331-y0);
f332=-(f_hist332-y0);
f333=-(f_hist333-y0);
f251=-(f_hist251-y0);
f252=-(f_hist252-y0);
f253=-(f_hist253-y0);
%% Normalize using analytical load capacity 
f14n=f14/Rn(1);
f24n=f24/Rn(2);
f34n=f34/Rn(3);
f44n=f44/Rn(4);
f331n=f331/Rn(5);
f332n=f332/Rn(6);
f333n=f333/Rn(7);
f251n=f251/Rn(8);
f252n=f252/Rn(9);
f253n=f253/Rn(10);

u14n=u14/(Rn(1)/Ki(1));
u24n=u24/(Rn(2)/Ki(2));
u34n=u34/(Rn(3)/Ki(3));
u44n=u44/(Rn(4)/Ki(4));
u331n=u331/(Rn(5)/Ki(5));
u332n=u332/(Rn(6)/Ki(6));
u333n=u333/(Rn(7)/Ki(7));
u251n=u251/(Rn(8)/Ki(8));
u252n=u252/(Rn(9)/Ki(9));
u253n=u253/(Rn(10)/Ki(10));

figure;
scatter(u14n,f14n,'y');
hold on
scatter(u14n,f14n,'m');
scatter(u24n,f24n,'c');
scatter(u34n,f34n,'r');
scatter(u44n,f44n,'g');
scatter(u331n,f331n,'b');
scatter(u332n,f332n,'k');
scatter(u333n,f333n);
scatter(u251n,f251n);
scatter(u252n,f252n);
scatter(u253n,f253n);

title('Normalized deformation vs force, normalized by analytical capacity ');
legend('L2*1/8','L2*3/16','L2*1/4','L2*5/16','L2.5*3/16','L2.5*1/4','L2.5*5/16','L3*3/16','L3*1/4','L3*5/16');
xlim([0,inf]);
ylim([0,1]);
%% Normalize and exclude the values after buckling and below 0 for curve fitting 
%C1_4 
[M,I] = max(f14n);
f14_cf=f14n(1:I); %exclude the values after buckling 
u14_cf=u14n(1:I); %exclude the values after buckling 

u14_cf=u14_cf(u14_cf>0); %exclude all the values below 0
Np=size(u14_cf,1); %number of points
f14_cf=f14_cf(end-Np+1:end);%exclude all the values below 0

figure;
scatter(u14_cf,f14_cf);
hold on
%C2_4 
[M,I] = max(f24n);
f24_cf=f24n(1:I); %exclude the values after buckling 
u24_cf=u24n(1:I); %exclude the values after buckling 

u24_cf=u24_cf(u24_cf>0); %exclude all the values below 0
Np=size(u24_cf,1); %number of points
f24_cf=f24_cf(end-Np+1:end);%exclude all the values below 0

scatter(u24_cf,f24_cf);
%C3_4 
[M,I] = max(f34n);
f34_cf=f34n(1:I); %exclude the values after buckling 
u34_cf=u34n(1:I); %exclude the values after buckling 

u34_cf=u34_cf(u34_cf>0); %exclude all the values below 0
Np=size(u34_cf,1); %number of points
f34_cf=f34_cf(end-Np+1:end);%exclude all the values below 0

scatter(u34_cf,f34_cf);
%C4_4 
[M,I] = max(f44n);
f44_cf=f44n(1:I); %exclude the values after buckling 
u44_cf=u44n(1:I); %exclude the values after buckling 

u44_cf=u44_cf(u44_cf>0); %exclude all the values below 0
Np=size(u44_cf,1); %number of points
f44_cf=f44_cf(end-Np+1:end);%exclude all the values below 0

scatter(u44_cf,f44_cf);
%C25_1 
[M,I] = max(f251n);
f251_cf=f251n(1:I); %exclude the values after buckling 
u251_cf=u251n(1:I); %exclude the values after buckling 

u251_cf=u251_cf(u251_cf>0); %exclude all the values below 0
Np=size(u251_cf,1); %number of points
f251_cf=f251_cf(end-Np+1:end);%exclude all the values below 0

scatter(u251_cf,f251_cf);
%C25_2 
[M,I] = max(f252n);
f252_cf=f252n(1:I); %exclude the values after buckling 
u252_cf=u252n(1:I); %exclude the values after buckling 

u252_cf=u252_cf(u252_cf>0); %exclude all the values below 0
Np=size(u252_cf,1); %number of points
f252_cf=f252_cf(end-Np+1:end);%exclude all the values below 0

scatter(u252_cf,f252_cf);
%C25_3 
[M,I] = max(f253n);
f253_cf=f253n(1:I); %exclude the values after buckling 
u253_cf=u253n(1:I); %exclude the values after buckling 

u253_cf=u253_cf(u253_cf>0); %exclude all the values below 0
Np=size(u253_cf,1); %number of points
f253_cf=f253_cf(end-Np+1:end);%exclude all the values below 0

scatter(u253_cf,f253_cf);
%C33_1 
[M,I] = max(f331n);
f331_cf=f331n(1:I); %exclude the values after buckling 
u331_cf=u331n(1:I); %exclude the values after buckling 

u331_cf=u331_cf(u331_cf>0); %exclude all the values below 0
Np=size(u331_cf,1); %number of points
f331_cf=f331_cf(end-Np+1:end);%exclude all the values below 0

scatter(u331_cf,f331_cf);
%C33_2 
[M,I] = max(f332n);
f332_cf=f332n(1:I); %exclude the values after buckling 
u332_cf=u332n(1:I); %exclude the values after buckling 

u332_cf=u332_cf(u332_cf>0); %exclude all the values below 0
Np=size(u332_cf,1); %number of points
f332_cf=f332_cf(end-Np+1:end);%exclude all the values below 0

scatter(u332_cf,f332_cf);
%C33_3 
[M,I] = max(f333n);
f333_cf=f333n(1:I); %exclude the values after buckling 
u333_cf=u333n(1:I); %exclude the values after buckling 

u333_cf=u333_cf(u333_cf>0); %exclude all the values below 0
Np=size(u333_cf,1); %number of points
f333_cf=f333_cf(end-Np+1:end);%exclude all the values below 0

scatter(u333_cf,f333_cf);
legend('L2*1/8','L2*3/16','L2*1/4','L2*5/16','L2.5*3/16','L2.5*1/4','L2.5*5/16','L3*3/16','L3*1/4','L3*5/16');
%% Normalize using abaqus maximum 

f14n=f14/RnAbaqus(1);
f24n=f24/RnAbaqus(2);
f34n=f34/RnAbaqus(3);
f44n=f44/RnAbaqus(4);
f331n=f331/RnAbaqus(5);
f332n=f332/RnAbaqus(6);
f333n=f333/RnAbaqus(7);
f251n=f251/RnAbaqus(8);
f252n=f252/RnAbaqus(9);
f253n=f253/RnAbaqus(10);

u14n=u14/(RnAbaqus(1)/Ki(1));
u24n=u24/(RnAbaqus(2)/Ki(2));
u34n=u34/(RnAbaqus(3)/Ki(3));
u44n=u44/(RnAbaqus(4)/Ki(4));
u331n=u331/(RnAbaqus(5)/Ki(5));
u332n=u332/(RnAbaqus(6)/Ki(6));
u333n=u333/(RnAbaqus(7)/Ki(7));
u251n=u251/(RnAbaqus(8)/Ki(8));
u252n=u252/(RnAbaqus(9)/Ki(9));
u253n=u253/(RnAbaqus(10)/Ki(10));

figure;
scatter(u14n,f14n,'y');
hold on
scatter(u14n,f14n,'m');
scatter(u24n,f24n,'c');
scatter(u34n,f34n,'r');
scatter(u44n,f44n,'g');
scatter(u331n,f331n,'b');
scatter(u332n,f332n,'k');
scatter(u333n,f333n);
scatter(u251n,f251n);
scatter(u252n,f252n);
scatter(u253n,f253n);

title('Normalized deformation vs force,normalized by abaqus capacity ');
legend('L2*1/8','L2*3/16','L2*1/4','L2*5/16','L2.5*3/16','L2.5*1/4','L2.5*5/16','L3*3/16','L3*1/4','L3*5/16');
xlim([0,inf]);
ylim([0,1]);
%% Normalize and exclude the values after buckling and below 0 for curve fitting 
%C1_4 
[M,I] = max(f14n);
f14_cf=f14n(1:I); %exclude the values after buckling 
u14_cf=u14n(1:I); %exclude the values after buckling 

u14_cf=u14_cf(u14_cf>0); %exclude all the values below 0
Np=size(u14_cf,1); %number of points
f14_cf=f14_cf(end-Np+1:end);%exclude all the values below 0

figure;
scatter(u14_cf,f14_cf);
hold on
%C2_4 
[M,I] = max(f24n);
f24_cf=f24n(1:I); %exclude the values after buckling 
u24_cf=u24n(1:I); %exclude the values after buckling 

u24_cf=u24_cf(u24_cf>0); %exclude all the values below 0
Np=size(u24_cf,1); %number of points
f24_cf=f24_cf(end-Np+1:end);%exclude all the values below 0

scatter(u24_cf,f24_cf);
%C3_4 
[M,I] = max(f34n);
f34_cf=f34n(1:I); %exclude the values after buckling 
u34_cf=u34n(1:I); %exclude the values after buckling 

u34_cf=u34_cf(u34_cf>0); %exclude all the values below 0
Np=size(u34_cf,1); %number of points
f34_cf=f34_cf(end-Np+1:end);%exclude all the values below 0

scatter(u34_cf,f34_cf);
%C4_4 
[M,I] = max(f44n);
f44_cf=f44n(1:I); %exclude the values after buckling 
u44_cf=u44n(1:I); %exclude the values after buckling 

u44_cf=u44_cf(u44_cf>0); %exclude all the values below 0
Np=size(u44_cf,1); %number of points
f44_cf=f44_cf(end-Np+1:end);%exclude all the values below 0

scatter(u44_cf,f44_cf);
%C25_1 
[M,I] = max(f251n);
f251_cf=f251n(1:I); %exclude the values after buckling 
u251_cf=u251n(1:I); %exclude the values after buckling 

u251_cf=u251_cf(u251_cf>0); %exclude all the values below 0
Np=size(u251_cf,1); %number of points
f251_cf=f251_cf(end-Np+1:end);%exclude all the values below 0

scatter(u251_cf,f251_cf);
%C25_2 
[M,I] = max(f252n);
f252_cf=f252n(1:I); %exclude the values after buckling 
u252_cf=u252n(1:I); %exclude the values after buckling 

u252_cf=u252_cf(u252_cf>0); %exclude all the values below 0
Np=size(u252_cf,1); %number of points
f252_cf=f252_cf(end-Np+1:end);%exclude all the values below 0

scatter(u252_cf,f252_cf);
%C25_3 
[M,I] = max(f253n);
f253_cf=f253n(1:I); %exclude the values after buckling 
u253_cf=u253n(1:I); %exclude the values after buckling 

u253_cf=u253_cf(u253_cf>0); %exclude all the values below 0
Np=size(u253_cf,1); %number of points
f253_cf=f253_cf(end-Np+1:end);%exclude all the values below 0

scatter(u253_cf,f253_cf);
%C33_1 
[M,I] = max(f331n);
f331_cf=f331n(1:I); %exclude the values after buckling 
u331_cf=u331n(1:I); %exclude the values after buckling 

u331_cf=u331_cf(u331_cf>0); %exclude all the values below 0
Np=size(u331_cf,1); %number of points
f331_cf=f331_cf(end-Np+1:end);%exclude all the values below 0

scatter(u331_cf,f331_cf);
%C33_2 
[M,I] = max(f332n);
f332_cf=f332n(1:I); %exclude the values after buckling 
u332_cf=u332n(1:I); %exclude the values after buckling 

u332_cf=u332_cf(u332_cf>0); %exclude all the values below 0
Np=size(u332_cf,1); %number of points
f332_cf=f332_cf(end-Np+1:end);%exclude all the values below 0

scatter(u332_cf,f332_cf);
%C33_3 
[M,I] = max(f333n);
f333_cf=f333n(1:I); %exclude the values after buckling 
u333_cf=u333n(1:I); %exclude the values after buckling 

u333_cf=u333_cf(u333_cf>0); %exclude all the values below 0
Np=size(u333_cf,1); %number of points
f333_cf=f333_cf(end-Np+1:end);%exclude all the values below 0

scatter(u333_cf,f333_cf);

legend('L2*1/8','L2*3/16','L2*1/4','L2*5/16','L2.5*3/16','L2.5*1/4','L2.5*5/16','L3*3/16','L3*1/4','L3*5/16');
%% Save
save('CurefitCompression','u14_cf','f14_cf','u24_cf','f24_cf',...
    'u34_cf','f34_cf','u44_cf','f44_cf','u251_cf','f251_cf','u252_cf','f252_cf',...
'u253_cf','f253_cf','u331_cf','f331_cf','u332_cf','f332_cf','u333_cf','f333_cf');








% %% save data for curve fitting  
% 
% %C33_1 is chosen to do the curve fitting 
% [M,I] = max(f252n);
% f252_cf=f252n(1:I); %exclude the values after buckling 
% u252_cf=u252n(1:I); %exclude the values after buckling 
% 
% u252_cf=u252_cf(u252_cf>0); %exclude all the values below 0
% Np=size(u252_cf,1); %number of points
% f252_cf=f252_cf(end-Np+1:end);%exclude all the values below 0
% 
% figure;
% scatter(u252_cf,f252_cf);
% hold on
% % scatter(u252n,f252n);
% save('CurefitCompression','u252_cf','f252_cf');

% 
% figure;
% plot(u14,f14);
% hold on
% plot(u24,f24);
% hold on
% plot(u34,f34);
% hold on
% plot(u44,f44);
% xlim([0,0.7]);
% xpoint1=0.01;
% delta1=(0:0.0001:xpoint1);
% delta2=(0:0.0001:xpoint1*fu24/fu14);
% delta3=(0:0.0001:xpoint1*fu34/fu14);
% delta4=(0:0.0001:xpoint1*fu44/fu14);
% y1=delta1*k14;
% y2=delta2*k24;
% y3=delta3*k34;
% y4=delta4*k44;
% hold on
% plot(delta1,y1);
% hold on
% plot(delta2,y2);
% hold on
% plot(delta3,y3);
% hold on
% plot(delta4,y4);
% %%
% %Normalize
% % y1max=6.229;
% % y2max=11.98;
% % y3max=18.25;
% % y4max=23.19;
% y1max=fu14;
% y2max=fu24;
% y3max=fu34;
% y4max=fu44;
% 
% x1max=0.1334;
% x2max=0.2683;
% x3max=0.4164;
% x4max=0.5294;
% 
% x1factor=y1max/k14;
% x2factor=y2max/k24;
% x3factor=y3max/k34;
% x4factor=y4max/k44;
% 
% 
% figure;
% scatter(u14/x1factor,f14/y1max);
% hold on
% scatter(u24/x2factor,f24/y2max);
% hold on
% scatter(u34/x3factor,f34/y3max);
% hold on
% scatter(u44/x4factor,f44/y4max);
% 
% u14n=u14/x1factor;f14n=f14/y1max; %nomalized disp and force
% u24n=u24/x2factor;f24n=f24/y2max;
% u34n=u34/x3factor;f34n=f34/y3max;
% u44n=u44/x4factor;f44n=f44/y4max;
% 
% save('C1_4n','u14n','f14n');
% save('C2_4n','u24n','f24n');
% save('C3_4n','u34n','f34n');
% save('C4_4n','u44n','f44n');
% 
% 
% %%
% %Richard Equation
% % def=(0:0.001:0.7);
% % ndef=def*0.3*k14/y1max;
% % %R1=y1max*((1.74*ndef)./(1+ndef.^0.5).^2-0.009*ndef);
% % R1=y1max*((3*ndef)./(1+ndef.^0.5).^2-0*ndef);
% % figure;
% % plot(u14,f14);
% % hold on
% % plot(u24,f24);
% % hold on
% % plot(u34,f34);
% % hold on
% % plot(u44,f44);
% % hold on
% % plot(def,R1);
% f44_f=f44/y4max; %factored value
% [M,I] = max(f44_f);
% f44_cf=f44_f(1:I); %exclude the values after buckling 
% u44_f=u44/x4factor; %factored value
% u44_cf=u44_f(1:I); %exclude the values after buckling 
% 
% u44_cf=u44_f(u44_cf>=0); %exclude all the values below 0
% Np=size(u44_cf,1); %number of points
% f44_cf=f44_cf(end-Np+1:end);%exclude all the values below 0
% 
% 
% 
% 
% 
% figure;
% scatter(u44_cf,f44_cf);
% save('curfitC4_4','u44_cf','f44_cf');
% 
% 
% 
% 
% 

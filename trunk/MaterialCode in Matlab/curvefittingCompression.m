clc;
close all;
clear all;
%% Step 1 nonlinear regression to get 4 parameters 
load('CurefitCompression.mat');
f = [f14_cf; f24_cf;f34_cf;f44_cf;f251_cf;f252_cf;f253_cf;f331_cf;f332_cf;f333_cf] ;
u = [u14_cf; u24_cf;u34_cf;u44_cf;u251_cf;u252_cf;u253_cf;u331_cf;u332_cf;u333_cf] ;

%K1---------a(1)
%Kp---------a(2)
%Ro---------a(3)
%n----------a(4)
mdl=@(a,x) (x.*a(1)./(1+(x.*a(1)./a(3)).^a(4)).^(1./a(4))+x.*a(2));
a0=[1.7;0;1.74;0.5];
[ahat,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(u,f,mdl,a0);

minX=min(u);
maxX=max(u);

figure;
gridX=linspace(minX,maxX,200);
plot(gridX, mdl(ahat,gridX),'r');
hold on 
scatter(u14_cf,f14_cf);
scatter(u24_cf,f24_cf);
scatter(u34_cf,f34_cf);
scatter(u44_cf,f44_cf);
scatter(u251_cf,f251_cf);
scatter(u252_cf,f252_cf);
scatter(u253_cf,f253_cf);
scatter(u331_cf,f331_cf);
scatter(u332_cf,f332_cf);
scatter(u333_cf,f333_cf);

%% Step 2 use the model to predict the behavior and compare to its original behavior in Abaqus
% initial point
y0=-6.4512*0.153*2;
x0=-0.0625+y0/2055.28*2;
%from C14 C24 C34 C44 C331 C332 C333 C251 C252 C253
Ki=[207.7922526
272.4905205
316.6087436
343.3311546
272.6907886
316.8866448
343.667171
272.6114264
316.7768887
343.5349261
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


% displacement 
u=(x0:-0.01:-0.7);
%C14
load('C1_4.mat');
u_hist14=U;
f_hist14=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(1)/Rn(1);
f14=-(mdl(ahat,def_norm)*Rn(1))+y0;

figure;
scatter(u_hist14,f_hist14);
hold on
plot(u,f14);

%C24
load('C2_4.mat');
u_hist24=U;
f_hist24=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(2)/Rn(2);
f24=-(mdl(ahat,def_norm)*Rn(2))+y0;

figure;
scatter(u_hist24,f_hist24);
hold on
plot(u,f24);

%C34
load('C3_4.mat');
u_hist34=U;
f_hist34=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(3)/Rn(3);
f34=-(mdl(ahat,def_norm)*Rn(3))+y0;

figure;
scatter(u_hist34,f_hist34);
hold on
plot(u,f34);

%C44
load('C4_4.mat');
u_hist44=U;
f_hist44=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(4)/Rn(4);
f44=-(mdl(ahat,def_norm)*Rn(4))+y0;

figure;
scatter(u_hist44,f_hist44);
hold on
plot(u,f44);

%C331
load('C33_1.mat');
u_hist331=U;
f_hist331=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(5)/Rn(5);
f331=-(mdl(ahat,def_norm)*Rn(5))+y0;

figure;
scatter(u_hist331,f_hist331);
hold on
plot(u,f331);

%C332
load('C33_2.mat');
u_hist332=U;
f_hist332=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(6)/Rn(6);
f332=-(mdl(ahat,def_norm)*Rn(6))+y0;

figure;
scatter(u_hist332,f_hist332);
hold on
plot(u,f332);

%C333
load('C33_3.mat');
u_hist333=U;
f_hist333=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(7)/Rn(7);
f333=-(mdl(ahat,def_norm)*Rn(7))+y0;

figure;
scatter(u_hist333,f_hist333);
hold on
plot(u,f333);

%C251
load('C25_1.mat');
u_hist251=U;
f_hist251=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(8)/Rn(8);
f251=-(mdl(ahat,def_norm)*Rn(8))+y0;

figure;
scatter(u_hist251,f_hist251);
hold on
plot(u,f251);

%C252
load('C25_2.mat');
u_hist252=U;
f_hist252=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(9)/Rn(9);
f252=-(mdl(ahat,def_norm)*Rn(9))+y0;

figure;
scatter(u_hist252,f_hist252);
hold on
plot(u,f252);

%C253
load('C25_3.mat');
u_hist253=U;
f_hist253=F;
%model predictted behavior 
def_norm=(-(u-x0))*Ki(10)/Rn(10);
f253=-(mdl(ahat,def_norm)*Rn(10))+y0;

figure;
scatter(u_hist253,f_hist253);
hold on
plot(u,f253);


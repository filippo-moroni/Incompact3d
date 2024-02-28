%% This is a code to plot statistics for BLs and channel flows

clc
close all
clear 

%% Input from the user
Re = 4200.0;              % Reynolds number (1/nu)
nu = 1/Re;

%% Latex interpreter
set(0,'defaulttextInterpreter','latex') 

%% External functions
%loadobj cmocean;
%loadobj contourfcmap;

%% Reading of file and variables
M = readtable('MEAN500.0.txt',NumHeaderLines=1);

% Averages of velocity components
mean_u = M{:,1};   % mean of u
mean_v = M{:,2};   % mean of v
mean_w = M{:,3};   % mean of w

% Variances of velocity components
var_u = M{:,4};    % variance of u
var_v = M{:,5};    % variance of v
var_w = M{:,6};    % variance of w

% Skewnesses of velocity components
skew_u = M{:,7};   % skewness of u
skew_v = M{:,8};   % skewness of v
skew_w = M{:,9};   % skewness of w

% Kurtoses of velocity components
kurt_u = M{:,10};  % kurtosis of u
kurt_v = M{:,11};  % kurtosis of v
kurt_w = M{:,12};  % kurtosis of w

% Reynolds stresses
mean_uv = M{:,13};  % <u'v'>
mean_uw = M{:,14};  % <u'w'>
mean_vw = M{:,15};  % <v'w'>

% Pressure and scalar fields
mean_p   = M{:,16};   % mean of p
var_p    = M{:,17};   % variance of p

mean_phi = M{:,18};   % mean of phi
var_phi  = M{:,19};   % variance of phi

% Mixed fluctuations (velocity and scalar)
mean_uphi = M{:,20};  % <u'phi'>
mean_vphi = M{:,21};  % <v'phi'>
mean_wphi = M{:,22};  % <w'phi'>

% Reading of grid points
G = readtable('yp.dat',NumHeaderLines=0);

y = G{:,1};           % y-coordinate at the faces of the cells  

%% Calculating shear velocity and viscous unit

% Mean gradient at the first face (shared by first 2 grid elements)
mean_gradient = mean_u(2)/y(2);    % partial U / partial y  

% Shear velocity
sh_vel = sqrt(nu*mean_gradient);

% Viscous unit
delta_nu = nu/sh_vel;

%% Rescaling variables through wall units
y = y/delta_nu;
mean_u = mean_u/sh_vel;
%% Plotting

% trial

h4 = figure;

%c = colorbar('FontSize',16);
%c.Label.String = '$U$';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize = 48;

plot(y,mean_u)

xlabel('$y^+$','FontSize',40)
ylabel('$U^+$','FontSize',40)
set(h4,'PaperSize',[40 18]);



















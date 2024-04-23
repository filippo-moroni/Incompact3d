%% This is a code to plot statistics for BLs and channel flows

clc
close all
clear 

%% Input from the user
Re = 4200.0;              % Reynolds number (1/nu)
nu = 1/Re;                % kinematic viscosity

ny = 65;                  % number of points in y-direction
nh = (ny-1)/2 + 1;        % half number of points in y-direction

Ly = 2.0;                 % total height of the channel

%% External functions
loadobj cloneYAxisFromLeftToRight;

%% Latex interpreter
set(0,'defaulttextInterpreter','latex') 

%% Set some useful colors
blue   = [57 106 177]./255;
red    = [204 37 41]./255;
black  = [83 81 84]./255;
green  = [62 150 81]./255;
brown  = [146 36 40]./255;
purple = [107 76 154]./255;

yellow = [0.9290 0.6940 0.1250];
orange = [0.8500 0.3250 0.0980];
lblue  = [0.3010 0.7450 0.9330];
grey   = [0.5 0.5 0.5];

%% Reading of file and variables

% CFR - default code & modified code
M1 = readtable('mean_stats400.0_default_extra_diss.txt',NumHeaderLines=1);
M2 = readtable('mean_stats.txt',NumHeaderLines=1);

% CPG
% M2 = readtable('mean_stats400.0_mycode_cpg.txt',NumHeaderLines=1);

% Vorticity - default code
M3 = readtable('vort_stats400.0.txt',NumHeaderLines=1);

%% Default code variables
mean_u  = M1{:,1};          % mean of u default code
mean_v  = M1{:,2};          % mean of v default code
var_u   = M1{:,4};          % variance of u
var_v   = M1{:,5};          % variance of v
mean_uv = M1{:,13};         % <u'v'>

%% Modified code variables
mean_u_mycode  = M2{:,1};   % mean of u my code
mean_v_mycode  = M2{:,2};   % mean of v my code
var_u_mycode   = M2{:,4};   % variance of u my code
var_v_mycode   = M2{:,5};   % variance of v my code
mean_uv_mycode = M2{:,13};  % <u'v'> my code

%% Vorticity
vort_x = M3{:,1};
vort_y = M3{:,2};
vort_z = M3{:,3};

%% Reading of grid points
G = readtable('yp.dat',NumHeaderLines=0);

y = G{:,1};           % y-coordinate at the faces of the cells 

%% Reference data by Lee & Moser
A1 = readtable('data_lee_retau180.txt',NumHeaderLines=72);

yplus_lee  = A1{3:end,2};
mean_u_lee = A1{3:end,3};

A2 = readtable('data_lee_fluct_retau180.txt',NumHeaderLines=75);

yplus_lee2  = A2{3:end,2};
var_u_lee   = A2{3:end,3};
var_v_lee   = A2{3:end,4};
mean_uv_lee = A2{3:end,6};

%% Calculations for default code

% Mean gradient at the first face (shared by first 2 grid elements)
mean_gradient = mean_u(2)/y(2);     % partial U / partial y 

% Shear velocity
sh_vel = sqrt(nu*mean_gradient);
 
% Viscous unit
delta_nu = nu/sh_vel;

% Viscous time 
t_nu = nu/(sh_vel^2);

%% Rescaling variables through wall units
y_plus_default = y/delta_nu;
mean_u  = mean_u/sh_vel;
var_u   = var_u/(sh_vel^2);
var_v   = var_v/(sh_vel^2);
mean_uv = mean_uv/(sh_vel^2);

vort_x  = vort_x*t_nu;
vort_y  = vort_y*t_nu;
vort_z  = vort_z*t_nu;

%% Calculations for modified code

% Mean gradient at the first face (shared by first 2 grid elements)
mean_gradient_mycode = mean_u_mycode(2)/y(2);     % partial U / partial y 

% Shear velocity
sh_vel_mycode = sqrt(nu*mean_gradient_mycode);
 
% Viscous unit
delta_nu_mycode = nu/sh_vel_mycode;

%% Rescaling variables through wall units
y_mycode = y/delta_nu_mycode;
mean_u_mycode  = mean_u_mycode/sh_vel_mycode;
var_u_mycode   = var_u_mycode/(sh_vel_mycode^2);
var_v_mycode   = var_v_mycode/(sh_vel_mycode^2);
mean_uv_mycode = mean_uv_mycode/(sh_vel_mycode^2);

%% Von Karman law 

% k = 0.41, B = 5.2 (Pope, "Turbulent flows")
% k = 0.37, B = 5.2 (Cimarelli, turb. lecture notes on channel flows)

B = 5.2;  
k = 0.37;

y_plus = linspace(5,180,175);
u_plus = (1/k)*log(y_plus)+B;

%% Viscous sub-layer

y_plus_vsl = linspace(1,15,15);
u_plus_vsl = y_plus_vsl;

%% Resizing arrays
y = y(1:nh);
y_mycode = y_mycode(1:nh);
y_plus_default = y_plus_default(1:nh);

mean_u  = mean_u (1:nh);
var_u   = var_u (1:nh);
var_v   = var_v (1:nh);
mean_uv = mean_uv(1:nh); 

mean_u_mycode  = mean_u_mycode (1:nh);
var_u_mycode   = var_u_mycode (1:nh);
var_v_mycode   = var_v_mycode (1:nh);
mean_uv_mycode = mean_uv_mycode(1:nh);

vort_x = vort_x (1:nh);
vort_y = vort_y (1:nh);
vort_z = vort_z (1:nh);

%% Mean velocity profile plot
h4 = figure;

scatter(y_plus_default,mean_u,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
hold on
semilogx(y_plus,u_plus,'Color',grey,'LineStyle', '--',LineWidth=1.5)
hold on
semilogx(yplus_lee,mean_u_lee,'Color',yellow,LineWidth=1.5)
scatter(y_mycode,mean_u_mycode,"MarkerEdgeColor",red,'Marker','o',LineWidth=1.5)
hold on
semilogx(y_plus_vsl,u_plus_vsl,'Color',grey,'LineStyle', '--',LineWidth=1.5)

legend({'Default Incompact3d', 'Viscous sublayer and log law', 'Lee and Moser (2015)','Modified Incompact3d'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',40)

yaxis_lim = 20;  % upper bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel('$U^+$','FontSize',40)
ylim([0,yaxis_lim]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
ylim([0,yaxis_lim]);

set(h4,'PaperSize',[22 12]);

caption = 'Log law with constants: k = 0.37, B = 5.2 (channel flow) \n(see lecture notes on turbulence prof. Cimarelli)';
text(0.05, -1, sprintf(caption), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');

%% Variance of u' plot
h4 = figure;

scatter(y_plus_default,var_u,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
hold on
semilogx(yplus_lee,var_u_lee,'Color',yellow,LineWidth=1.5)
scatter(y_mycode,var_u_mycode,"MarkerEdgeColor",red,'Marker','o',LineWidth=1.5)

legend({'Default Incompact3d', 'Lee and Moser (2015)','Modified Incompact3d'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',40)

yaxis_lim = 8;  % upper bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel("$\langle u'^2 \rangle/u_\tau^2$",'FontSize',40)
ylim([0,yaxis_lim]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
ylim([0,yaxis_lim]);

set(h4,'PaperSize',[22 12]);

%% Variance of v' plot
h4 = figure;

scatter(y_plus_default,var_v,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
hold on
semilogx(yplus_lee,var_v_lee,'Color',yellow,LineWidth=1.5)
scatter(y_mycode,var_v_mycode,"MarkerEdgeColor",red,'Marker','o',LineWidth=1.5)

legend({'Default Incompact3d', 'Lee and Moser (2015)','Modified Incompact3d'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',40)

yaxis_lim = 0.8;  % upper bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel("$\langle v'^2 \rangle/u_\tau^2$",'FontSize',40)
ylim([0,yaxis_lim]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
ylim([0,yaxis_lim]);

set(h4,'PaperSize',[22 12]);

%% Reynolds stresses <u'v'> plot
h4 = figure;

scatter(y_plus_default,mean_uv,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
hold on
semilogx(yplus_lee,mean_uv_lee,'Color',yellow,LineWidth=1.5)
scatter(y_mycode,mean_uv_mycode,"MarkerEdgeColor",red,'Marker','o',LineWidth=1.5)

legend({'Default Incompact3d', 'Lee and Moser (2015)','Modified Incompact3d'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',40)

yaxis_lim = -0.9;  % lower bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel("$\langle u'v' \rangle/u_\tau^2$",'FontSize',40)
ylim([yaxis_lim,0.1]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
ylim([yaxis_lim,0.1]);

set(h4,'PaperSize',[22 12]);

%% Vorticity plot
h4 = figure;

scatter(y_plus_default,vort_y,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)

legend({'Default Incompact3d', 'Lee and Moser (2015)','Modified Incompact3d'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',40)

%yaxis_lim = -0.9;  % lower bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel("$\langle \omega_i \rangle t_\nu$",'FontSize',40)
%ylim([yaxis_lim,0.1]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
%ylim([yaxis_lim,0.1]);

set(h4,'PaperSize',[22 12]);


























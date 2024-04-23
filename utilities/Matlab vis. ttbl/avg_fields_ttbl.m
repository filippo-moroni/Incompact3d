%% This is a code to plot statistics for temporal TBLs

clc
close all
clear 

%% Input from the user
Re = 500.0;               % Reynolds number (1/nu)
nu = 1/Re;                % kinematic viscosity

ny = 649;                 % number of points in y-direction
Ly = 24.0;                % total height of the domain

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

% Mean stats - modified code
M1 = readtable('mean_stats430.0.txt',NumHeaderLines=1);

% Vorticity - modified code
M3 = readtable('vort_stats430.0.txt',NumHeaderLines=1);

%% Default code variables
mean_u  = M1{:,1};          % mean of u default code
mean_v  = M1{:,2};          % mean of v default code
var_u   = M1{:,4};          % variance of u
var_v   = M1{:,5};          % variance of v
mean_uv = M1{:,13};         % <u'v'>

uwall = 1.0;
mean_u = uwall - mean_u;

%% Vorticity
vort_x = M3{:,1};
vort_y = M3{:,2};
vort_z = M3{:,3};

%% Reading of grid points
G = readtable('yp.dat',NumHeaderLines=0);

y = G{:,1};           % y-coordinate at the faces of the cells 

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
y_plus = y/delta_nu;
mean_u  = mean_u/sh_vel;
var_u   = var_u/(sh_vel^2);
var_v   = var_v/(sh_vel^2);
mean_uv = mean_uv/(sh_vel^2);

vort_x  = vort_x*t_nu;
vort_y  = vort_y*t_nu;
vort_z  = vort_z*t_nu;

%% Viscous sub-layer & Von Karman law 

% Viscous sublayer
y_plus_vsl = linspace(1,15,15);
u_plus_vsl = y_plus_vsl;

% Values from Kozul et al. (2016)
k = 0.384;
B = 4.173;  

% Log region
y_plus_k = linspace(5,180,175);  % k: Kozul et al.
u_plus_k = (1/k)*log(y_plus_k)+B;

%% Mean velocity profile plot
h4 = figure;

scatter(y_plus,mean_u,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
hold on
semilogx(y_plus_vsl,u_plus_vsl,'Color',grey,'LineStyle', '--',LineWidth=1.5)
hold on
semilogx(y_plus_k,u_plus_k,'Color',grey,'LineStyle', '--',LineWidth=1.5)

legend({'Present', 'Viscous sublayer and log law (Kozul et al. (2016))'}, 'Interpreter', 'latex',Location='northwest',FontSize=18);

grid on;
grid minor;

xlim([0,180]);
xticks([0 5 30 60 100 180])
set(gca,'xscale','log')
xlabel('$y^+$','FontSize',50)

yaxis_lim = 20;  % upper bound of y axes

yyaxis left
ax = gca;
ax.YColor = 'black'; 
ylabel('$U^+$','FontSize',50)
ylim([0,yaxis_lim]);
yyaxis right
ax = gca;
ax.YColor = 'black'; 
ylim([0,yaxis_lim]);

set(h4,'PaperSize',[22 12]);

caption = 'Log law with constants: k = 0.384, B = 4.173 (Kozul et al. (2016))';
text(0.32, -1, sprintf(caption), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 16, 'FontWeight', 'bold');


%% Vorticity plot
% 
% h4 = figure;
% 
% scatter(y_plus,vort_z,"MarkerEdgeColor",blue,'Marker','o',LineWidth=1.5)
% 
% legend({'Present'}, 'Interpreter', 'latex',Location='northwest',FontSize=12);
% 
% grid on;
% grid minor;
% 
% xlim([0,y_plus(ny)]);
% xticks([0 1 5 30 60 100 180, 300, y_plus(ny)])
% set(gca,'xscale','log')
% xlabel('$y^+$','FontSize',40)
% 
% yaxis_lim = 20;  % upper bound of y axes
% 
% yyaxis left
% ax = gca;
% ax.YColor = 'black'; 
% ylabel('$\langle \omega_z \rangle$','FontSize',40)
% %ylim([0,yaxis_lim]);
% yyaxis right
% ax = gca;
% ax.YColor = 'black'; 
% %ylim([0,yaxis_lim]);
% 
% set(h4,'PaperSize',[22 12]);





























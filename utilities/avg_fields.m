%% This is a code to plot statistics for BLs and channel flows

clc
close all
clear 

%% Input from the user
Re = 4200.0;              % Reynolds number (1/nu)

ny = 65;                  % number of points in y-direction
nh = (ny-1)/2 + 1;        % half number of points in y-direction

Ly = 2.0;                 % total height of the channel

x_vertical = 180;         % friction Reynolds number of the simulation

%% Latex interpreter
set(0,'defaulttextInterpreter','latex') 

%% External functions
%loadobj cmocean;
%loadobj contourfcmap;

%% Reading of file and variables

% CFR
M  = readtable('mean_stats400.0_default_extra_diss.txt',NumHeaderLines=1);
M2 = readtable('mean_stats400.0_mycode_extra_diss.txt',NumHeaderLines=1);

% CPG
% M = readtable('mean_stats400.0_mycode_cpg.txt',NumHeaderLines=1);

% Averages of velocity components
mean_u = M{:,1};         % mean of u default code
mean_u_mycode = M2{:,1};  % mean of u my code

% Reading of grid points
G = readtable('yp.dat',NumHeaderLines=0);

y = G{:,1};           % y-coordinate at the faces of the cells 

%% Calculations

% Bulk velocity calculation and rescaling
Ub = sum(mean_u(1:ny))/Ly;
Ub = sqrt(Ub);

mean_u = mean_u/Ub;

nu = 1/Re;

% Mean gradient at the first face (shared by first 2 grid elements)
mean_gradient = mean_u(2)/y(2);     % partial U / partial y 

% Shear velocity
sh_vel = sqrt(nu*mean_gradient);
 
% Viscous unit
delta_nu = nu/sh_vel;

%% Rescaling variables through wall units
y_plus_default = y/delta_nu;
mean_u = mean_u/sh_vel;


%% Bulk velocity calculation and rescaling
Ub = sum(mean_u_mycode(1:ny))/Ly;
Ub = sqrt(Ub);

mean_u_mycode = mean_u_mycode/Ub;

nu = 1/Re;

% Mean gradient at the first face (shared by first 2 grid elements)
mean_gradient_mycode = mean_u_mycode(2)/y(2);     % partial U / partial y 

% Shear velocity
sh_vel_mycode = sqrt(nu*mean_gradient_mycode);
 
% Viscous unit
delta_nu_mycode = nu/sh_vel_mycode;

%% Rescaling variables through wall units
y_mycode = y/delta_nu_mycode;
mean_u_mycode = mean_u_mycode/sh_vel_mycode;

%% Data by Lee & Moser
M = readtable('data_lee_retau180.txt',NumHeaderLines=72);

yplus_lee = M{3:end,2};
mean_u_lee = M{3:end,3};
%% Plotting

% k = 0.41, B = 5.2 data reported by Pope (Turbulent flows)
B = 5.2;  
k = 0.37; % for a channel flow (more specific constant, see turbulence notes 
          % by prof. Cimarelli)

y_plus = linspace(1,180,180);
u_plus = (1/k)*log(y_plus)+B;

h4 = figure;

y = y(1:nh);

y_mycode = y_mycode(1:nh);
y_plus_default = y_plus_default(1:nh);

mean_u = mean_u(1:nh);
mean_u_mycode = mean_u_mycode(1:nh);

%plot(y,mean_u)
%hold on
%plot(y,u_plus)

semilogx(y_plus_default,mean_u,LineWidth=1.5)
hold on
semilogx(y_plus,u_plus,LineWidth=1.5)
hold on
semilogx(yplus_lee,mean_u_lee,LineWidth=1.5)
semilogx(y_mycode,mean_u_mycode,LineWidth=1.5)
line([x_vertical, x_vertical], ylim, 'Color', 'black', 'LineStyle', '--','Linewidth',1.5); 
legend({'default code', 'log law', 'Lee and Moser','mycode','$Re_\tau$'}, 'Interpreter', 'latex',Location='northwest');

xlim([0,240]);
xticks([0 5 30 60 100 180])

grid on;
grid minor;

xlabel('$y^+$','FontSize',40)
ylabel('$U^+$','FontSize',40)
set(h4,'PaperSize',[16 16]);






























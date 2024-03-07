%% Setup parameters for temporal BLs

clc
close all
clear

% Derivative at the wall

nu = 0.000002;          % kinematic viscosity

D = 0.001;              % trip wire diameter
uwall = 1.0;            % fixed velocity of the wall (U_wall)
thetasl = 54*nu/uwall;  % thickness of the shear layer (theta_sl)


mg = - uwall/(4*thetasl)*(sech(D/2/thetasl))^2;  % mean gradient of the IC
sh_vel = sqrt(nu*abs(mg));  % shear velocity at the start of the simulation


ReD = uwall*D/nu;           % trip Reynolds number
delta_nu = nu/sh_vel;       % viscous unit
Re = 1/nu;                  % Reynolds used in Incompact3d
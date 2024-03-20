%% Setup parameters for temporal BLs

clc
close all
clear

% Inputs
nu = 0.002;                  % kinematic viscosity
D = 1.0;                     % trip wire diameter
uwall = 1.0;                 % fixed velocity of the wall (U_wall)
cf = 0.007;                  % maximum cf of the temporal BL (Cimarelli et al.(2024))


% Outputs
sh_vel = sqrt((cf/2))*uwall; % corresponding shear velocity 
ReD = uwall*D/nu;            % trip Reynolds number
thetasl = 54*nu/uwall;       % thickness of the shear layer (theta_sl)
delta_nu = nu/sh_vel;        % viscous unit
Re = 1/nu;                   % Reynolds used in Incompact3d


%% Derivative at the wall

% mg = - uwall/(4*thetasl)*(sech(D/2/thetasl))^2; % mean gradient of the IC
% sh_vel = sqrt(nu*abs(mg));                      % shear velocity at the start of the simulation
                                                  
% the real value is much higher (e.g. about 0.059 vs 0.0013)
% see Cimarelli et al. (2024)
% real reference values are thus needed!

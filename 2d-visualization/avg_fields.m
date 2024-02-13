%% This is a code to plot 2d averaged fields

clc
close all
clear 

%% Latex interpreter
set(0,'defaulttextInterpreter','latex') 

%% External functions
loadobj cmocean;
%loadobj contourfcmap;

%% Grey color for the airfoil
htmlGray = [200 200 200]/255;

%% Imporing the airfoil geometry, creation of polyshape and rotation
airfoil = readtable('naca_4412.txt',NumHeaderLines=0);
xa = airfoil{:,1};
ya = airfoil{:,2};

pgon = polyshape(xa,ya);
pgon_rot = rotate(pgon,-15);

%% Reading of the file and of variables
M = readtable('Time+Span_Average_tvmd131.txt',NumHeaderLines=1);
x = M{:,1};          % x-coordinate
y = M{:,2};          % y-coordinate

U = M{:,4};          % U velocity component
V = M{:,5};          % V velocity component

vortz = M{:,22};     % <omega_z> component
% Iturb = M{:,23};     % I_turb[%]

% Only for ILES
% vortz = M{:,13};     % <omega_z> component
% Iturb = M{:,14};     % I_turb[%]

% k = M{:,12};         % <k>
% uu = M{:,8};         % <u'u'>
% vv = M{:,9};         % <v'v'>
% ww = M{:,10};        % <w'w'>
% uv = M{:,11};        % <u'v'>
% 
tau11 = M{:,13};     % <tau_11>
tau22 = M{:,14};     % <tau_22>
tau33 = M{:,15};     % <tau_33>
tau12 = M{:,16};     % <tau_12>

esgs = M{:,19};      % eps_sgs
Esgs = M{:,20};      % E_sgs
esgs_p = M{:,21};    % eps_sgs'

% Preparing the additional grid

dx = 0.0005;
dy = 0.0005;
xg = -0.1:dx:1.1;
yg = -0.35:dy:0.25;
yg = yg';

% Plotting U and <omega_z>

% U velocity component
Ui = griddata(x,y,U,xg,yg);
h4 = figure;

ll = -1;
ul =  2;

contourf(xg,yg,Ui,ll:0.1:ul,'LineColor','none')
c = colorbar('FontSize',16);
c.Label.String = '$U$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.0)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);


% % <omega_z> 
% 
% dx = 0.0001;
% dy = 0.0001;
% xg = -0.05:dx:0.2;
% yg = -0.05:dy:0.1;
% yg = yg';
% 
% vortzi = griddata(x,y,vortz,xg,yg);
% h4 = figure;
% 
% ll = -400;
% ul =  400;
% 
% constrained_data = max(min(vortzi, ul), ll);
% axis equal
% contourf(xg,yg,constrained_data,ll:8:ul,'LineColor','none')
% xlim([-0.04 0.2])
% ylim([-0.03 0.1])
% c = colorbar('FontSize',16);
% c.Label.String = '$\langle\omega_z\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% clim([ll ul]);
% 
% hold on
% cmocean('bal','pivot',0.0)
% hold on
% 
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlim([-0.04 0.2])
% ylim([-0.03 0.1])
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);


%% Plotting variances, Reynolds stresses and <k>

% % <k>
% ki = griddata(x,y,k,xg,yg);
% h4 = figure;
% 
% ll = 0.0;
% ul = 0.4;
% 
% contourf(xg,yg,ki,ll:0.03:ul,'LineColor','none')
% c = colorbar('FontSize',16);
% c.Label.String = '$\langle k\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% clim([ll ul]);
% hold on
% cmocean('haline')
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);

% % <u'u'>
% uui = griddata(x,y,uu,xg,yg);
% h4 = figure;
% 
% ll = 0.0;
% ul = 0.6;
% 
% contourf(xg,yg,uui,ll:0.04:ul,'LineColor','none')
% c = colorbar('FontSize',16);
% c.Label.String = '$\langle u^{\prime 2}\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% clim([ll ul]);
% hold on
% cmocean('haline')
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);

%% <v'v'> and <w'w'> must be finalized (they are not needed probably)

% % <v'v'>
% vvi = griddata(x,y,vv,xg,yg);
% h4 = figure;
% 
% ll = 0.0;
% ul = 0.1;
% 
% contourf(xg,yg,vvi,ll:0.02:ul,'LineColor','none')
% c = colorbar('FontSize',16);
% c.Label.String = '$\langle v^{\prime 2}\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% clim([ll ul]);
% hold on
% cmocean('haline')
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);
% 
% % <w'w'>
% wwi = griddata(x,y,ww,xg,yg);
% figure
% contourf(xg,yg,wwi,0:0.04:0.4,'LineColor','none')
% c = colorbar;
% c.Label.String = '$\langle w^{\prime 2}\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 14;
% cmocean('haline')
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlabel('$x/c$')
% ylabel('$y/c$')

%%

% % <u'v'>
% uvi = griddata(x,y,uv,xg,yg);
% h4 = figure;
% 
% ll = -0.1;     % lower limit to show
% ul =  0.3;     % upper limit to show
% 
% constrained_data = max(min(uvi, ul), ll);
% contourf(xg,yg,constrained_data,ll:0.02:ul,'LineColor','none')
% 
% c = colorbar('FontSize',16);
% c.Label.String = '$\langle u^{\prime} v^{\prime}\rangle$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% clim([ll ul]);
% hold on
% cmocean('bal','pivot',0.0)
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);

%% Subgrid stresses

% <tau11>
tau11i = griddata(x,y,tau11,xg,yg);
h4=figure;

ll = -0.001;
ul =  0.001;

contourf(xg,yg,tau11i,ll:0.0001:ul,'LineColor','none')
c = colorbar('FontSize',16);
c.Label.String = '$\langle \tau_{11}\rangle$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48; 
clim([ll ul]);
hold on
cmocean('bal','pivot',0.0)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);


% <tau22>
tau22i = griddata(x,y,tau22,xg,yg);
h4=figure;

ll = -0.001;
ul =  0.001;

contourf(xg,yg,tau22i,ll:0.00005:ul,'LineColor','none')
c = colorbar('FontSize',16);
c.Label.String = '$\langle \tau_{22}\rangle$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.0)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);
% 
% % % <tau33>
% % tau33i = griddata(x,y,tau33,xg,yg);
% % h4=figure;
% % contourf(xg,yg,tau33i,-0.0001:0.00001:0.0001,'LineColor','none')
% % c = colorbar('FontSize',16);
% % c.Label.String = '$\langle \tau_{33}\rangle$';
% % c.Label.Interpreter = 'latex';
% % c.Label.FontSize = 48;
% % hold on
% % plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
% % xlabel('$x/c$','FontSize',40)
% % ylabel('$y/c$','FontSize',40)
% % set(h4,'PaperSize',[40 18]);
% % 
% <tau12>
tau12i = griddata(x,y,tau12,xg,yg);
h4=figure;

ll = -0.003;     % lower limit to show
ul =  0.003;     % upper limit to show

constrained_data = max(min(tau12i, ul), ll);
contourf(xg,yg,constrained_data,ll:0.00005:ul,'LineColor','none')

c = colorbar('FontSize',16);
c.Label.String = '$\langle \tau_{12}\rangle$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.0)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);

%% Streamlines

% x_stream =  -0.1:0.01:1.1;
% y_stream =  -0.3:0.01:0.1;
% y_stream = y_stream';
% 
% [startX,startY] = meshgrid(x_stream,y_stream);
% 
% Ustream = griddata(x,y,U,x_stream,y_stream);
% Vstream = griddata(x,y,V,x_stream,y_stream);
% 
% figure
% verts = stream2(x_stream,y_stream,Ustream,Vstream,startX,startY);
% lineobj = streamline(verts);
% hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)


% %% Plotting streamlines and airfoil
% figure
% % streamslice(x_stream,y_stream,Ustream,Vstream,10.0) 
% % hold on
% plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1,LineWidth=1.0)
% axis equal

%%
% figure
% quiver(x,y,U,V);

%% Subgrid dissipation rates

% eps_sgs
esgsi = griddata(x,y,esgs,xg,yg);
h4=figure;

ll = -0.7;     % lower limit to show
ul =  0.7;     % upper limit to show

constrained_data = max(min(esgsi, ul), ll);
contourf(xg,yg,constrained_data,ll:0.05:ul,'LineColor','none')

c = colorbar('FontSize',16);
c.Label.String = '$ \varepsilon_{sgs}$'; 
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.00001)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);

% eps_prime_sgs
esgs_pi = griddata(x,y,esgs_p,xg,yg);
h4=figure;

ll = -0.7;     % lower limit to show
ul =  0.4;     % upper limit to show

constrained_data = max(min(esgs_pi, ul), ll);
contourf(xg,yg,constrained_data,ll:0.05:ul,'LineColor','none')

c = colorbar('FontSize',16);
c.Label.String = '$ \varepsilon_{sgs}^\prime$'; 
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.00001)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);

% E_sgs
Esgsi = griddata(x,y,Esgs,xg,yg);
h4=figure;

ll = -0.5;     % lower limit to show
ul =  0.5;     % upper limit to show

constrained_data = max(min(Esgsi, ul), ll);
contourf(xg,yg,constrained_data,ll:0.05:ul,'LineColor','none')

c = colorbar('FontSize',16);
c.Label.String = '$E_{sgs}$'; 
c.Label.Interpreter = 'latex';
c.Label.FontSize = 48;
clim([ll ul]);
hold on
cmocean('bal','pivot',0.00001)
hold on
plot(pgon_rot,'FaceColor', htmlGray,FaceAlpha=1)
xlabel('$x/c$','FontSize',40)
ylabel('$y/c$','FontSize',40)
set(h4,'PaperSize',[40 18]);

%% Turbulence intensity

% % Preparing the additional grid, only for the turbulence intensity
% 
% dx = 0.0005;
% dy = 0.0005;
% xg = -0.1:dx:0.0;
% yg = -0.2:dy:0.2;
% yg = yg';
% 
% % I_turb[%]
% Iturbi = griddata(x,y,Iturb,xg,yg);
% h4 = figure;
% contourf(xg,yg,Iturbi,0:0.25:6,'LineColor','none')
% c = colorbar('FontSize',16);
% c.Label.String = '$I_\infty[\%]$';
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 48;
% axis equal
% cmocean('haline')
% xlabel('$x/c$','FontSize',40)
% ylabel('$y/c$','FontSize',40)
% set(h4,'PaperSize',[40 18]);

%% Calculation of the turbulence intensity on the monitoring plane

% dy = 0.0005;
% xg = -0.1;
% yg = -0.2:dy:0.2;
% yg = yg';
% 
% Iturbi = griddata(x,y,Iturb,xg,yg);
% I_turb_mean = mean(Iturbi);




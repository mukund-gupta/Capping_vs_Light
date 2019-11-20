clear all; clc; close all;

% Adding libraries
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library/AutoDiff_R2016b')
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_channel\library\AutoDiff_R2016b');
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_channel\library');
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library');

my_pwd = pwd();

% Constants
yr = 3600*360*24; % one year [s]
Mc = 12;
rho0 = 1030;

%% Model data

nyears = 100;
tt = 1:(12*nyears);

% Grid stuff
load('grid2d.mat')
XC = grid2d.XC;
YC = grid2d.YC;
RC = grid2d.RC;
DRF = grid2d.DRF;
DXG = grid2d.DXG;
DYC = grid2d.DYC;
ylat_min = -69.5; % I am not sure about this value
yc = mean(YC,1);
YC_lat = ylat_min + YC/(111*10^3);
yc_lat = ylat_min + yc/(111*10^3);
ny = length(yc);

av_idx = length(yc)/2;
idx_st = 1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation

my_runs = {'par0_ice0','seas_sens/par0_ice25seas','seas_sens/par0_ice50seas','seas_sens/par0_ice75seas','seas_sens/par0_ice100seas'};
alpha_vent = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
   
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_vent(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_vent_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light

my_runs = {'par0_ice0','seas_sens/par25seas_ice0','seas_sens/par50seas_ice0','seas_sens/par75seas_ice0','seas_sens/par100seas_ice0'};
alpha_light = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_light(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_light_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation and light

my_runs = {'par0_ice0','seas_sens/par25seas_ice25seas','seas_sens/par50seas_ice50seas','seas_sens/par75seas_ice75seas','seas_sens/par100seas_ice100seas'};
alpha_both = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_both(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_both_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;
    
end

%%
% Plotting

figure
width = 615;
height = 250;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);
mksz = 130;
clrs = [237, 141, 14; 12, 160, 22; 0,0,0]/255;

subplot(1,2,1)
sc1 = scatter(alpha_vent,Fint_vent,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
hold on
sc2 = scatter(alpha_light,Fint_light,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
sc3 = scatter(alpha_both,Fint_both,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
xlabel('Sea ice fraction')
ylabel('Outgassing carbon flux (gC.m^{-1}.yr^{-1})')
legend([sc1 sc2 sc3],{'Capping','Light Attenuation','Capping + Light Atten.'},'Location','NorthWest')
ylim([0 4e7])
title('(a) Ice zone')

subplot(1,2,2)
sc1 = scatter(alpha_vent,Fint_vent_tot,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
hold on
sc2 = scatter(alpha_light,Fint_light_tot,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
sc3 = scatter(alpha_both,Fint_both_tot,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
xlabel('Sea ice fraction')
ylabel('Integrated flux (gC.m^{-1}.yr^{-1})')
% legend([sc1 sc2 sc3],{'Capping','Light Attenuation','Capping + Light Atten.'},'Location','NorthWest')
title('(b) Whole domain')


%%
% Functions

function F = getFbio(alpha,ke,kb,V,H,L,Cin)
    F = -ke*Cin*V*H ./ (ke + kb*(1-alpha)) .* ( exp(-L/V/H*(ke+kb*(1-alpha))) -1 );
end

function F = getFblock(alpha,ke,kb,V,H,L,Cin)
    F = -ke*(1-alpha)*Cin*V*H ./ (ke*(1-alpha) + kb) .* ( exp(-L/V/H*(ke*(1-alpha) + kb)) -1 );
end

function F = getFtot(alpha,ke,kb,V,H,L,Cin)
    F = -ke*(1-alpha)*Cin*V*H ./ (ke*(1-alpha) + kb*(1-alpha)) .* ( exp(-L/V/H*(ke*(1-alpha) + kb*(1-alpha))) -1 );
end

function F = getFnorm(alpha,ke,kb,V,H,L,Cin)
    lambda_0 = (ke + kb)*L/V/H;
    lambda = lambda_0*(1-alpha);   
    F = (exp(-lambda) - 1)/(exp(-lambda_0) - 1);
end

function [F,C] = getFlat(alpha_e,alpha_b,ke,kb,V,H,L,Ci,y)
    lambda = (ke*(1-alpha_e) + kb*(1-alpha_b))*L/V/H;
    C = Ci*exp(-lambda*y/L);
    F = ke*(1-alpha_e)*C;
end


function [dFdke,dFdkb] = getGradfe(alpha,ke,kb,V,H,L,Cin)
    f = @(ke_v,kb_v) -ke_v*(1-alpha)*Cin*V*H ./ (ke_v*(1-alpha) + kb) .* ( exp(-L/V/H*(ke_v*(1-alpha) + kb)) -1 );
    [ke_v,kb_v] = ainit(ke,kb,3);   
    z = f(ke_v,kb_v);                               
    dFdke = z{1,0};
    dFdkb = z{0,1};
end

function [dFdke,dFdkb] = getGradfb(alpha,ke,kb,V,H,L,Cin)
    f = @(ke_v,kb_v) -ke_v*Cin*V*H ./ (ke_v + kb_v*(1-alpha)) .* ( exp(-L/V/H*(ke_v+kb_v*(1-alpha))) -1 );
    [ke_v,kb_v] = ainit(ke,kb,3);   
    z = f(ke_v,kb_v);                               
    dFdke = z{1,0};
    dFdkb = z{0,1};
end

function supertitle(mytitle,varargin)

axes('Units','Normal');
h = title(mytitle,varargin{1:length(varargin)});
set(gca,'visible','off')
set(h,'visible','on')

end
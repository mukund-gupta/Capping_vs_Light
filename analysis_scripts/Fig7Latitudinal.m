clear all; clc; close all;

my_pwd = pwd();

% Constants
yr = 360*3600*24;
Mc = 12;
rho0 = 1030;
SI_salt = 4/1000; % salt concentration in sea ice (kg salt / kg water)
pCO2a = 278*1e-6;
clrs = [255,0,0; 0,255,0; 0,0,255; 255,165,0]/255;

% Adding libraries
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_channel\library');
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library');

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

% % Control simulations
% cd('../ctrl/sim_1')
% load('DICPCO2.mat')
% load('DICCFLX.mat') % [mol/m2/sec]
% 
% load('SALT.mat')
% load('THETA.mat')
% load('VVELMASS.mat')
% load('GM_PsiY.mat')
% load('SI_Fract.mat')
% cd(my_pwd)
% 
% THETA_c = THETA;
% SALT_c = SALT;
% GM_PsiY_c = GM_PsiY;
% VVELMASS_c = VVELMASS;
% SI_Fract_c = SI_Fract;
% DICPCO2_c = DICPCO2;
% DICCFLX_c = DICCFLX;
% [~,~,psiR_c] = get_psi(VVELMASS_c,GM_PsiY_c,DRF,DXG);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plotting
% No ice effect on par 
my_runs = {'par0_ice0','par0_ice25','par0_ice50','par0_ice75','par0_ice100'};
my_labels = {'0%','25%','50%','75%','100%'};

figure (1)
width = 650;
height = 550;
set(gcf,'color','w')
set(gcf,'units','points','position',[10,10,width,height]) % Display size
plt_idx = 300;

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd('../bin_files')
    load('sea_ice.mat')
    cd(my_pwd)
    
    alpha_av_nopar(i) = mean(mean(sea_ice(idx_st:av_idx,:),2),1);
    Fint_nopar(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
    dpCO2 = (1-mean(sea_ice,2)) .* (mean(DICPCO2(:,tt),2) - pCO2a)*1e6;
    dpCO2_av_nopar(i) = mean(mean(dpCO2(idx_st:av_idx,:),2),1);
    
    % Plotting
    
    figure (1)
    subplot(3,2,1)
    fld = -mean(DICCFLX(:,tt),2)*yr*Mc;
    plot(yc_lat(1:plt_idx),fld(1:plt_idx),'linewidth',1.5,'Displayname',my_labels{i})
%     xlabel('latitude')
    ylabel('Carbon flux (gCm^{-2}yr^{-1})')
    title('Channel - Capping')
    hold on
    ylim([-20 55])
    
end

subplot(3,2,1)
legend show

% supertitle(['2D channel - Ice affects only blocking (not PAR)' newline newline],'FontSize',14)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No ice effect on blocking 

my_runs = {'par0_ice0','par25_ice0','par50_ice0','par75_ice0','par100_ice0'};
my_labels = {'0%','25%','50%','75%','100%'};

figure (1)
plt_idx = 300;

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])

    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd('../bin_files')
    load('BCs.mat')
    cd(my_pwd)
    
    alpha_av_noblock(i) = mean(mean(par_ice(idx_st:av_idx,:),2),1);
    Fint_noblock(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
    dpCO2 = (1-mean(par_ice,2)) .* (mean(DICPCO2(:,tt),2) - pCO2a)*1e6;
    dpCO2_av_noblock(i) = mean(mean(dpCO2(idx_st:av_idx,:),2),1);
    
    % Plotting
    figure (1)
    subplot(3,2,3)
    fld = -mean(DICCFLX(:,tt),2)*yr*Mc;
    plot(yc_lat(1:plt_idx),fld(1:plt_idx),'linewidth',1.5,'Displayname',my_labels{i})
%     xlabel('latitude')
    ylabel('Carbon flux (gCm^{-2}yr^{-1})')
    title('Channel - Light Attenuation')
    hold on
    ylim([-20 55])
    
    
end

% subplot(3,2,3)
% legend show
%grid on

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both effects on blocking and par

my_runs = {'par0_ice0','par25_ice25','par50_ice50','par75_ice75','par100_ice100'};
my_labels = {'0%','25%','50%','75%','100%'};

figure (1)
plt_idx = 300;

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd('../bin_files')
    load('BCs.mat')
    cd(my_pwd)
    
    alpha_av(i) = mean(mean(par_ice(idx_st:av_idx,:),2),1);
    Fint(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
    dpCO2 = (1-mean(par_ice,2)) .* (mean(DICPCO2(:,tt),2) - pCO2a)*1e6;
    dpCO2_av(i) = mean(mean(dpCO2(idx_st:av_idx,:),2),1);
    
    % Plotting
    figure (1)
    subplot(3,2,5)
    fld = -mean(DICCFLX(:,tt),2)*yr*Mc;
    plot(yc_lat(1:plt_idx),fld(1:plt_idx),'linewidth',1.5,'Displayname',my_labels{i})
    xlabel('latitude (°)')
    ylabel('Carbon flux (gCm^{-2}yr^{-1})')
    title('Channel - Capping and Light Atten.')
    hold on
    ylim([-20 55])
    
end
% 
% subplot(3,2,5)
% legend show
% grid on

supertitle(['2D channel - Ice affects blocking and PAR' newline newline],'FontSize',14)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D analytical Model

% Plotting parameters
xmin = 0;
xmax = 1650;

% Model parameters
load('fit_params.mat')

% Grid and timestep parameters
dy = 1*10^3; % meridional spacing [m]
N = floor(L/dy); % Number of grid cells
y_ana = linspace(dy/2,(N-1)*dy,N);
y_ana_lat = -65 + y_ana/(110*1e3);

yi = 150*1000;
[~,y_idx] = min( abs( y_ana-yi ) );

% Setting sea ice concentration in all simulations
alpha_all = [0 0.25 0.5 0.75 1];
Nsim = length(alpha_all);
my_labels = {'0%','25%','50%','75%','100%'};


% Flux only
C_all = zeros(N,Nsim);
Fe_all = zeros(N,Nsim);
Fav_flux = zeros(1,Nsim);

for i=1:Nsim
    % Initialization
   
    alpha = alpha_all(:,i);
    ke_mod = ke*(1-alpha);
    kb_mod = kb;
    lambda = (ke_mod + kb_mod)*L / (V*H); % flux only
    C = Cin*exp(-lambda*y_ana/L);
    Fe = -ke_mod*C; % flux only
    
    % Average flux calculation
%     Fav_flux(i) = -ke_mod*Cin/lambda*L*(exp(-lambda)-1);
    Fav_flux(i) = -trapz(y_ana,Fe);
    C_all(:,i) = C;
    Fe_all(:,i) = Fe;
    
    % Plotting
    figure (1)
    
    subplot(3,2,2)
    plot(y_ana/1000,-Fe,'linewidth',1.5)
    hold on
    ylim([0 110])

end

subplot(3,2,2)
% set(gca, 'XDir','reverse')
% xlabel('Meridional distance (km)')
% ylabel('Air-sea carbon flux (+ve positive)')
title('Analytical - Capping')
ylabel('Carbon flux (gCm^{-2}yr^{-1})')
xlim([xmin xmax])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bio only - analytical 

C_all = zeros(N,Nsim);
Fe_all = zeros(N,Nsim);
Fav_bio = zeros(1,Nsim);

for i=1:Nsim
    % Initialization
   
    alpha = alpha_all(:,i);
    ke_mod = ke;
    kb_mod = kb*(1-alpha);
    lambda = (ke_mod + kb_mod)*L / (V*H); % flux only
    C = Cin*exp(-lambda*y_ana/L);
    Fe = -ke_mod*C; % flux only
    
    % Average flux calculation
%     Fav_bio(i) = -ke_mod*Cin/lambda*L*(exp(-lambda)-1);
    Fav_bio(i) = -trapz(y_ana,Fe);
    C_all(:,i) = C;
    Fe_all(:,i) = Fe;
    
    % Plotting
    figure (1)
    
    subplot(3,2,4)
    plot(y_ana/1000,-Fe,'linewidth',1.5)
    hold on
    ylim([0 110])

end

subplot(3,2,4)
% set(gca, 'XDir','reverse')
% xlabel('Meridional distance (km)')
% ylabel('Air-sea carbon flux (+ve positive)')
title('Analytical - Light Attenuation')
ylabel('Carbon flux (gCm^{-2}yr^{-1})')
xlim([xmin xmax])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PAR and blocking - Analytical
C_all = zeros(N,Nsim);
Fe_all = zeros(N,Nsim);
Fav_both = zeros(1,Nsim);

for i=1:Nsim
    % Initialization
   
    alpha = alpha_all(:,i);
    ke_mod = ke*(1-alpha);
    kb_mod = kb*(1-alpha);
    lambda = (ke_mod + kb_mod)*L / (V*H); % flux only
    C = Cin*exp(-lambda*y_ana/L);
    Fe = -ke_mod*C; % flux only
    
    % Average flux calculation
%     Fav_both(i) = -ke_mod*Cin/lambda*L*(exp(-lambda)-1);
    Fav_both(i) = -trapz(y_ana,Fe);
    C_all(:,i) = C;
    Fe_all(:,i) = Fe;
    
    % Plotting
    figure (1)
    
    subplot(3,2,6)
    plot(y_ana/1000,-Fe,'linewidth',1.5)
    hold on
    ylim([0 110])
    

end

subplot(3,2,6)
% set(gca, 'XDir','reverse')
xlabel('Meridional distance (km)')
ylabel('Carbon flux (gCm^{-2}yr^{-1})')
title('Analytical - Capping and Light Atten.')
xlim([xmin xmax])

%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions

function map = bmap1(cmin,cmax,cint)
    n=floor(abs( (cmax-cmin)  /cint));
    r = linspace(1,0,n);
    g = linspace(1,0,n);
    b = linspace(1,0.5,n);
    map = [r(:), g(:), b(:)];
end

function map = bmap2(cmin,cmax,cint)
    n=floor(abs( (cmax-cmin)  /cint));
    r = linspace(1,1,n);
    g = linspace(0,0,n);
    b = linspace(0,0,n);
    map = [r(:), g(:), b(:)];
end

function [ map ] = bmap_mono(cmin,cmax,cint)
    n=floor(abs( (cmax-cmin)  /cint));
    r = linspace(255,0,n)/255;
    g = linspace(0,0,n)/255;
    b = linspace(0,255,n)/255;
    map = [r(:), g(:), b(:)];
end

function [ map ] = bluepalered(cmin,cmax,cint)
    n=floor(abs(cmin/cint));
    m=floor(abs(cmax/cint));
    r = [linspace(0,0.8,n) linspace(1,1,m)];
    g = [linspace(0.3,1,n) linspace(1,0.3,m)];
    b = [linspace(1,1,n) linspace(0.8,0,m)];
    map = [r(:), g(:), b(:)];
end

function [ map ] = bluepalered2(cmin,cmax,cint_c,cint_w)
    n=floor(abs(cmin/cint_c));
    m=floor(abs(cmax/cint_w));
    r = [linspace(0,0.8,n) linspace(1,1,m)];
    g = [linspace(0.3,1,n) linspace(1,0.3,m)];
    b = [linspace(1,1,n) linspace(0.8,0,m)];
    map = [r(:), g(:), b(:)];
end

function fld_mean = area_mean(fld,DRF,DXG,msk,nyears)
    fld_mean = zeros(1,nyears);
    tot_area = 0;

    for i=1:length(DXG)
        for k=1:length(DRF)
            tot_area = tot_area + DRF(k)*DXG(i)*msk(i,k);
            for t_idx =1:nyears
                fld_mean(t_idx) = fld_mean(t_idx) + fld(i,k,t_idx)*DRF(k)*DXG(i)*msk(i,k);
            end
        end
    end
    
    fld_mean = fld_mean/tot_area;
    
    
end

function [psiE,psiB,psiR] = get_psi(VVELMASS,GM_PsiY,DRF,DXG)

    % Getting psi
    nt = size(VVELMASS,3);
    ny = size(VVELMASS,1);
    nz = size(VVELMASS,2);
    for tt=1:nt

        % Velocities
        vvel = VVELMASS(:,:,tt);
        gm_psi = zeros(ny,nz+1);
        gm_psi(:,1:nz) = GM_PsiY(:,:,tt);
        vvel_bol = zeros(ny,nz);

        for k=1:nz
           vvel_bol(:,k) = ( gm_psi(:,k+1) - gm_psi(:,k) ) / DRF(k);
        end

        % Integrating to get streamfunction
        varE = zeros(ny,nz+1);
        varB = zeros(ny,nz+1);

        for k=nz:-1:1
            varE(:,k) = varE(:,k+1) - DRF(k)*DXG'.*vvel(:,k);
            varB(:,k) = varB(:,k+1) - DRF(k)*DXG'.*vvel_bol(:,k);
        end
        psiE(:,:,tt) = varE(:,1:nz);
        psiB(:,:,tt) = varB(:,1:nz);
        psiR(:,:,tt) = (psiB(:,:,tt) + psiE(:,:,tt)); % Calculating residual
    end

end


clear all; clc; close all;

my_pwd = pwd();

% Library
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library')
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_channel\library')

% Constants
yr = 360*3600*24;
Mc = 12;
rho0 = 1030;
SI_salt = 4/1000; % salt concentration in sea ice (kg salt / kg water)
pCO2a = 278*1e-6;
clrs = [255,0,0; 0,255,0; 0,0,255; 255,165,0]/255;

%% Channel model

% Grid stuff
load('grid2d.mat')
XC = grid2d.XC;
YC = grid2d.YC;
RC = grid2d.RC;
DRF = grid2d.DRF;
DXG = grid2d.DXG;
DYC = grid2d.DYC;
dyc = DYC(1);
ylat_min = -69.5; % I am not sure about this value
yc = mean(YC,1);
YC_lat = ylat_min + YC/(111*10^3);
lat_chan = ylat_min + yc/(111*10^3);
ny = length(yc);

% Loading control data
cd('../ctrl_both/sim_1')
load('SI_Fract.mat')
load('DICCFLX.mat')
cd(my_pwd)

% Getting seasonality
% flx_clim_chan = -squeeze(MonthlyMean(DICCFLX))*yr*Mc; % positive out of ocean
[flx_clim_chan,~,~,flx_std_chan] = get_clim_stats(-DICCFLX*yr*Mc);


flx_chan_DJF = (flx_clim_chan(:,1) + flx_clim_chan(:,2) + flx_clim_chan(:,12))/3;
flx_chan_MAM = mean(flx_clim_chan(:,3:5),2);
flx_chan_JJA = mean(flx_clim_chan(:,6:8),2);
flx_chan_SON = mean(flx_clim_chan(:,9:11),2);
flx_chan_AM = mean(flx_clim_chan,2);

std_chan_DJF = nanmean(flx_std_chan(:,[1:2 12]),2);
std_chan_MAM = nanmean(flx_std_chan(:,3:5),2);
std_chan_JJA = nanmean(flx_std_chan(:,6:8),2);
std_chan_SON = nanmean(flx_std_chan(:,9:11),2);
std_chan_AM = nanmean(flx_std_chan,2);

alpha_chan = squeeze(MonthlyMean(SI_Fract));
alpha_chan_DJF = (alpha_chan(:,1) + alpha_chan(:,2) + alpha_chan(:,12))/3;
alpha_chan_MAM = mean(alpha_chan(:,3:5),2);
alpha_chan_JJA = mean(alpha_chan(:,6:8),2);
alpha_chan_SON = mean(alpha_chan(:,9:11),2);
alpha_chan_AM = mean(alpha_chan,2);

%% Landschutzer

cd('../../Landschutzer_data')
fname = 'spco2_1982-2015_MPI_SOM-FFN_v2016.nc';
flx = ncread(fname,'fgco2_smoothed')*Mc;
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
seamask = double(ncread(fname,'seamask'));
cd(my_pwd)

flx(flx>1e10) = NaN;
[flx_clim,flx_max,flx_min,flx_std] = get_clim_stats(flx);
flx_clim2 = flx_clim.*repmat(seamask,1,1,12);

flx_obs_zon = squeeze(nanmean(flx_clim,1));
flx_obs_DJF = nanmean(flx_obs_zon(:,[1:2 12]),2);
flx_obs_MAM = nanmean(flx_obs_zon(:,3:5),2);
flx_obs_JJA = nanmean(flx_obs_zon(:,6:8),2);
flx_obs_SON = nanmean(flx_obs_zon(:,9:11),2);
flx_obs_AM = nanmean(flx_obs_zon,2);

std_obs_zon = squeeze(nanmean(flx_std,1));
std_obs_DJF = nanmean(std_obs_zon(:,[1:2 12]),2);
std_obs_MAM = nanmean(std_obs_zon(:,3:5),2);
std_obs_JJA = nanmean(std_obs_zon(:,6:8),2);
std_obs_SON = nanmean(std_obs_zon(:,9:11),2);
std_obs_AM = nanmean(std_obs_zon,2);

%% Takahashi
% Getting Takahashi data
cd('../../Takahashi_data')
M = csvread('takahashi_data.csv',1,0);
cd(my_pwd)

lon_r = M(:,2);
lat_r = M(:,1);
mths_r = M(:,3);
area_r = M(:,16); % is the area of each box.  Unit is millions of squarekilometers, or 10^6 km^2 
ice_r = M(:,17); 
flx_r = M(:,19); % Monthly CO2 flux in Moles Carbon per square meter. Unit is Mole Carbon m^-2 month^-1.
N = length(lon_r);

lon_tak = 2.5:5:(360-2.5);
lat_tak = -76:4:4;
mths_tak = 1:1:12;

nlon_tak = length(lon_tak);
nlat_tak = length(lat_tak);
nmths_tak = length(mths_tak);

flx_tak = NaN*zeros(nlon_tak,nlat_tak,nmths_tak);
area_tak = NaN*zeros(nlon_tak,nlat_tak,nmths_tak);
ice_tak = NaN*zeros(nlon_tak,nlat_tak,nmths_tak);

for i=1:N
    ii = find(lon_tak==lon_r(i));
    jj = find(lat_tak==lat_r(i));
    kk = find(mths_tak==mths_r(i));
    flx_tak(ii,jj,kk) = flx_r(i);
    area_tak(ii,jj,kk) = area_r(i);
    ice_tak(ii,jj,kk) = ice_r(i);
end

flx_zon_tak = squeeze(nanmean(flx_tak,1))*12; % Mole Carbon m^-2 yr^-1
flx_tak_DJF = (flx_zon_tak(:,1) + flx_zon_tak(:,2) + flx_zon_tak(:,12))/3;
flx_tak_MAM = mean(flx_zon_tak(:,3:5),2);
flx_tak_JJA = mean(flx_zon_tak(:,6:8),2);
flx_tak_SON = mean(flx_zon_tak(:,9:11),2);
flx_tak_AM = nanmean(flx_zon_tak,2);

%% Getting sea ice fraction observations
cd('../../../SO_eddies/observational_data')
fname = 'seas_clim_2011_2015.nc';
alpha_obs = ncread(fname,'alpha');
lon_obs = ncread(fname,'lon');
lat_obs = ncread(fname,'lat');
cd(my_pwd)

alpha_obs = squeeze(nanmean(alpha_obs,1));
alpha_obs_DJF = nanmean([alpha_obs(:,1) alpha_obs(:,2) alpha_obs(:,12)],2);
alpha_obs_MAM = nanmean(alpha_obs(:,3:5),2);
alpha_obs_JJA = nanmean(alpha_obs(:,6:8),2);
alpha_obs_SON = nanmean(alpha_obs(:,9:11),2);
alpha_obs_AM = nanmean(alpha_obs,2);


%% Plotting

x0=50;
y0=50;
width=410;
height=520;
fnt_size=11;
fnt_name = 'Times New Roman';
ln_width = 0.5;

    % params for flux plot (channel)
ax_xLow_flx_chan = 55;
ax_yLow_flx_chan = 45;
ax_xLength_flx_chan = 265;
ax_yLength_flx_chan = 200;

    % params for ice plots (channel)
ax_xLow_ice_chan = ax_xLow_flx_chan+1;
ax_yLow_ice_chan = ax_yLow_flx_chan + 5;
ax_xLength_ice_chan = ax_xLength_flx_chan-2;
ax_yLength_ice_chan = 10;
yoffs_ice_chan = ax_yLength_ice_chan + 7;

    % params for flux plot (obs)
ax_xLow_flx_obs = 55;
ax_yLow_flx_obs = ax_yLow_flx_chan + ax_yLength_flx_chan + 45;
ax_xLength_flx_obs = 265;
ax_yLength_flx_obs = 200;

    % params for ice plots (obs)
nz_ice = 10;
z_ice = 1:nz_ice;
ax_xLow_ice_obs = ax_xLow_flx_obs+1;
ax_yLow_ice_obs = ax_yLow_flx_obs + 5;
ax_xLength_ice_obs = ax_xLength_flx_obs-2;
ax_yLength_ice_obs = 10;
yoffs_ice_obs = ax_yLength_ice_obs + 7;


figure
xmin = -68.5;
xmax = -41.5;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

% Observations
    % fluxes
    
err_idx = 5 + 1:5:180;

ax(5)=axes;
set(gca,'units','points','position',[ax_xLow_flx_obs,ax_yLow_flx_obs,ax_xLength_flx_obs,ax_yLength_flx_obs])

plot(lat_tak,flx_tak_AM,'k:','linewidth',3,'DisplayName','Annual Mean (Tak)')
hold on
plot(lat,flx_obs_AM,'k','linewidth',3,'DisplayName','Annual Mean (Lan)')
hold on
plt3 = plot(lat,flx_obs_DJF,'color',clrs(1,:),'linewidth',1.5,'DisplayName','DJF: Summer');
plt4 = plot(lat,flx_obs_MAM,'color',clrs(2,:),'linewidth',1.5,'DisplayName','MAM: Fall');
plt1 = plot(lat,flx_obs_JJA,'color',clrs(3,:),'linewidth',1.5,'DisplayName','JJA: Winter');
plt2 = plot(lat,flx_obs_SON,'color',clrs(4,:),'linewidth',1.5,'DisplayName','SON: Spring');

plt_er0 = errorbar(lat(err_idx),flx_obs_AM(err_idx),std_obs_AM(err_idx),'k','LineStyle','none','CapSize',10,'lineWidth',1.5);
plt_er1 = errorbar(lat(err_idx),flx_obs_DJF(err_idx),std_obs_DJF(err_idx),'color',clrs(1,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
plt_er2 = errorbar(lat(err_idx),flx_obs_MAM(err_idx),std_obs_MAM(err_idx),'color',clrs(2,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
plt_er3 = errorbar(lat(err_idx),flx_obs_JJA(err_idx),std_obs_JJA(err_idx),'color',clrs(3,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
plt_er4 = errorbar(lat(err_idx),flx_obs_SON(err_idx),std_obs_SON(err_idx),'color',clrs(4,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);

plot(lat,0*lat,'k--')
% xlabel('Latitude (°)')
ylabel('CO_2 flux [gCm^{-2}yr^{-1}]')
legend('Annual (Tak.)','Annual (Lan.)','DJF: Summer','MAM: Fall','JJA: Winter','SON: Spring')
title('Observed','FontSize',12)
xlim([xmin xmax])
ylim([-30 30])

annotation('textbox','String','[+ve outgassing]','units','points','linestyle','none','Position',[ax_xLow_flx_obs+10,ax_yLow_flx_obs-10,ax_xLength_flx_obs,ax_yLength_flx_obs],'FontSize',11)

    % SON 
ax(4)=axes;
fld = repmat(alpha_obs_SON,1,nz_ice); fld(fld<0.05) = NaN; 
[C1,h1] = contourf(lat_obs,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(4),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_obs,ax_yLow_ice_obs,ax_xLength_ice_obs,ax_yLength_ice_obs])
xlim([xmin xmax])
annotation('textbox','String','SON','units','points','linestyle','none','Position',[ax_xLow_ice_obs+110,ax_yLow_ice_obs + 5,ax_xLength_ice_obs,ax_yLength_ice_obs])

    % JJA 
ax(3)=axes;
fld = repmat(alpha_obs_JJA,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_obs,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(3),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_obs,ax_yLow_ice_obs+ yoffs_ice_obs,ax_xLength_ice_obs,ax_yLength_ice_obs])
xlim([xmin xmax])
annotation('textbox','String','JJA','units','points','linestyle','none','Position',[ax_xLow_ice_obs+110,ax_yLow_ice_obs+ yoffs_ice_obs + 5,ax_xLength_ice_obs,ax_yLength_ice_obs])

    % MAM 
ax(2)=axes;
fld = repmat(alpha_obs_MAM,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_obs,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(2),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_obs,ax_yLow_ice_obs + 2*yoffs_ice_obs,ax_xLength_ice_obs,ax_yLength_ice_obs])
xlim([xmin xmax])
annotation('textbox','String','MAM','units','points','linestyle','none','Position',[ax_xLow_ice_obs+67,ax_yLow_ice_obs + 2*yoffs_ice_obs + 4,ax_xLength_ice_obs,ax_yLength_ice_obs])

    % DJF 
ax(1)=axes;
fld = repmat(alpha_obs_DJF,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_obs,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(1),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_obs,ax_yLow_ice_obs + 3*yoffs_ice_obs,ax_xLength_ice_obs,ax_yLength_ice_obs])
xlim([xmin xmax])
annotation('textbox','String','DJF','units','points','linestyle','none','Position',[ax_xLow_ice_obs+55,ax_yLow_ice_obs + 3*yoffs_ice_obs + 5,ax_xLength_ice_obs,ax_yLength_ice_obs])



% % Channel model fluxes

err_idx = 50 + 1:60:320;
ax(10)=axes;
% smth_sp = 30;
smth_sp = 1;
set(gca,'units','points','position',[ax_xLow_flx_chan,ax_yLow_flx_chan,ax_xLength_flx_chan,ax_yLength_flx_chan])
plot(lat_chan,flx_chan_AM,'k','linewidth',3,'DisplayName','Annual Mean')
hold on
plt3 = plot(lat_chan,smooth(flx_chan_DJF,smth_sp),'color',clrs(1,:),'linewidth',1.5,'DisplayName','DJF: Summer');
plt4 = plot(lat_chan,smooth(flx_chan_MAM,smth_sp),'color',clrs(2,:),'linewidth',1.5,'DisplayName','MAM: Fall');
plt1 = plot(lat_chan,smooth(flx_chan_JJA,smth_sp),'color',clrs(3,:),'linewidth',1.5,'DisplayName','JJA: Winter');
plt2 = plot(lat_chan,smooth(flx_chan_SON,smth_sp),'color',clrs(4,:),'linewidth',1.5,'DisplayName','SON: Spring');


% plt_er0 = errorbar(lat_chan(err_idx),flx_chan_AM(err_idx),std_chan_AM(err_idx),'k','LineStyle','none','CapSize',10,'lineWidth',1.5);
% plt_er1 = errorbar(lat_chan(err_idx),flx_chan_DJF(err_idx),std_chan_DJF(err_idx),'color',clrs(1,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
% plt_er2 = errorbar(lat_chan(err_idx),flx_chan_MAM(err_idx),std_chan_MAM(err_idx),'color',clrs(2,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
% plt_er3 = errorbar(lat_chan(err_idx),flx_chan_JJA(err_idx),std_chan_JJA(err_idx),'color',clrs(3,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);
% plt_er4 = errorbar(lat_chan(err_idx),flx_chan_SON(err_idx),std_chan_SON(err_idx),'color',clrs(4,:),'LineStyle','none','CapSize',10,'lineWidth',1.5);


plot(lat_chan,0*lat_chan,'k--')
xlabel('Latitude (°)')
ylabel('CO_2 flux [gCm^{-2}yr^{-1}]')
% legend('Annual mean','DJF: Summer','MAM: Fall','JJA: Winter','SON: Spring')
title('Simulated','FontSize',12)
xlim([xmin xmax])
ylim([-50 50])

    % SON 
ax(11)=axes;
fld = repmat(alpha_chan_SON,1,nz_ice); fld(fld<0.05) = NaN; 
[C1,h1] = contourf(lat_chan,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(11),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_chan,ax_yLow_ice_chan,ax_xLength_ice_chan,ax_yLength_ice_chan])
xlim([xmin xmax])
annotation('textbox','String','SON','units','points','linestyle','none','Position',[ax_xLow_ice_chan+120,ax_yLow_ice_chan + 3,ax_xLength_ice_chan,ax_yLength_ice_chan])

    % JJA 
ax(12)=axes;
fld = repmat(alpha_chan_JJA,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_chan,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(12),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_chan,ax_yLow_ice_chan+ yoffs_ice_chan,ax_xLength_ice_chan,ax_yLength_ice_chan])
xlim([xmin xmax])
annotation('textbox','String','JJA','units','points','linestyle','none','Position',[ax_xLow_ice_chan+120,ax_yLow_ice_chan+ yoffs_ice_chan + 4,ax_xLength_ice_chan,ax_yLength_ice_chan])

    % MAM 
ax(13)=axes;
fld = repmat(alpha_chan_MAM,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_chan,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(13),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_chan,ax_yLow_ice_chan + 2*yoffs_ice_chan,ax_xLength_ice_chan,ax_yLength_ice_chan])
xlim([xmin xmax])
annotation('textbox','String','MAM','units','points','linestyle','none','Position',[ax_xLow_ice_chan+75,ax_yLow_ice_chan + 2*yoffs_ice_chan + 3,ax_xLength_ice_chan,ax_yLength_ice_chan])

    % DJF 
ax(14)=axes;
fld = repmat(alpha_chan_DJF,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(lat_chan,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(14),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow_ice_chan,ax_yLow_ice_chan + 3*yoffs_ice_chan,ax_xLength_ice_chan,ax_yLength_ice_chan])
xlim([xmin xmax])
annotation('textbox','String','DJF','units','points','linestyle','none','Position',[ax_xLow_ice_chan+75,ax_yLow_ice_chan + 3*yoffs_ice_chan + 5,ax_xLength_ice_chan,ax_yLength_ice_chan])


% colorbar
hcb = colorbar;
set(hcb,'units','points','position',[ax_xLow_flx_obs + ax_xLength_flx_obs + 15,ax_yLow_flx_chan + ax_xLength_flx_chan/2,20,ax_yLength_flx_obs])
ylabel(hcb,'Sea ice fraction')




% figure
% plot(yc_lat',flx_chan_AM,'k','linewidth',2.5)
% hold on
% plot(yc_lat',flx_chan_DJF,'color',clrs(1,:),'linewidth',1.5)
% hold on
% plot(yc_lat,flx_chan_MAM,'color',clrs(2,:),'linewidth',1.5)
% hold on
% plot(yc_lat,flx_chan_JJA,'color',clrs(3,:),'linewidth',1.5)
% hold on
% plot(yc_lat,flx_chan_SON,'color',clrs(4,:),'linewidth',1.5)
% % set(gca,'units','points','position',[ax_xLow,ax_yLow,ax_xLength,ax_yLength])
% xlabel('Latitude (°)')
% ylabel('(b) CO_2 flux (gCm^{-2}yr^{-1})')
% title('2d channel flux (+ve out of ocean)')
% legend('Annual mean','DJF: Summer','MAM: Fall','JJA: Winter','SON: Spring')
% xlim([-70 -41.5])
% ylim([-70 40])


%% Functions

function [ map ] = ice_map(cmin,cmax,cint)
    n=floor(abs(cmax-cmin)/cint)/2;
%     rgb_st = [239, 247, 249]/255;
    rgb_st = [255, 255, 255]/255;
    rgb_mid = [72, 194, 242]/255;
    rgb_end = [18, 6, 155]/255;
    
    r = [linspace(rgb_st(1),rgb_mid(1),n) linspace(rgb_mid(1),rgb_end(1),n)];
    g = [linspace(rgb_st(2),rgb_mid(2),n) linspace(rgb_mid(2),rgb_end(2),n)];
    b = [linspace(rgb_st(3),rgb_mid(3),n) linspace(rgb_mid(3),rgb_end(3),n)];
    map = [r(:), g(:), b(:)];
end

function [v_mean,v_max,v_min,v_std] = get_clim_stats(v)

    sl = 12; % length of slice (num of months in a year)
    dim = length(size(v))+1; % dimension of input array
    
    m = size(v);
    m2 = [m(1:end-1) sl m(end)/sl];
    v2 = reshape(v,m2);
    
    % Getting stats
    v_mean = squeeze(mean(v2,dim));
    v_max = squeeze(max(v2,[],dim));
    v_min = squeeze(min(v2,[],dim));
    v_std = squeeze(std(v2,0,dim));

end




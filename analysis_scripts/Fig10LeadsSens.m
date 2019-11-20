clear all; clc; close all;

my_pwd = pwd();

% Library
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_ctrlnel/library')
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_ctrlnel\library')

% Constants
yr = 360*3600*24;
Mc = 12;
rho0 = 1030;
SI_salt = 4/1000; % salt concentration in sea ice (kg salt / kg water)
pCO2a = 278*1e-6;
clrs = [255,0,0; 0,255,0; 0,0,255; 255,165,0]/255;

clrs2 = [242, 160, 7; 87, 199, 237]/255;

%% Control simulation

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
lat = ylat_min + yc/(111*10^3);
ny = length(yc);

idx_st = 1;
idx_en = 160;

% Loading control data
cd('../par0_ice0/sim_1')
load('DICCFLX.mat')
cd(my_pwd)

flx = -squeeze(MonthlyMean(DICCFLX))*yr*Mc; % positive out of ocean
flx_ctrl = mean(flx,2);
Fi_ctrl = trapz(yc(idx_st:idx_en),flx_ctrl(idx_st:idx_en))

%% par0_ice75

cd('../par0_ice75/sim_1')
load('DICCFLX.mat')
cd(my_pwd)

flx = -squeeze(MonthlyMean(DICCFLX))*yr*Mc; % positive out of ocean
flx_75 = mean(flx,2);
Fi_75 = trapz(yc(idx_st:idx_en),flx_75(idx_st:idx_en))



%% par0_ice75fr

cd('../par0_ice75fr/sim_1')
load('DICCFLX.mat')
cd(my_pwd)

flx = -squeeze(MonthlyMean(DICCFLX))*yr*Mc; % positive out of ocean
flx_75fr = mean(flx,2);
Fi_75fr = trapz(yc(idx_st:idx_en),flx_75fr(idx_st:idx_en))


%% Plotting

figure
width = 280;
height = 225;
ax_xLow = 47;
ax_yLow = 37;
ax_xLength = 210;
ax_yLength = 170;


set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);
set(gca,'units','points','position',[ax_xLow,ax_yLow,ax_xLength,ax_yLength])
plt0 = plot(lat,flx_ctrl,'k','linewidth',3.5,'DisplayName','Annual mean');
hold on
plt3 = plot(lat,flx_75,'color',[11, 132, 48]/255,'linewidth',2.5);
plt1 = plot(lat,flx_75fr,':','color',[11, 132, 48]/255,'linewidth',2.5);

plot(lat,lat*0,'k--')

xlim([-69 -42])
ylim([-20 50])
xlabel('Latitude (°)')
ylabel('CO_2 flux [gCm^{-2}yr^{-1}]')
legend('No ice','75% uniform ice','75% intermittent ice')
annotation('textbox','String','[+ve outgassing]','units','points','linestyle','none','Position',[ax_xLow+110,ax_yLow-50,ax_xLength,ax_yLength],'FontSize',11)


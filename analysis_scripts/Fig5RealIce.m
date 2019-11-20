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

% Model data

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

% Getting PAR
cd('../ctrl_par/bin_files')
fname = 'par_monthly.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
cd(my_pwd)
PAR = reshape(var,ny,12); 
PAR_DJF = mean(PAR(:,[12 1:2]),2);
PAR_MAM = mean(PAR(:,3:5),2);
PAR_JJA = mean(PAR(:,6:8),2);
PAR_SON = mean(PAR(:,9:11),2);

% Ventilation only
cd('../ctrl_ice/sim_1')
%load('DICPCO2.mat')
load('DICCFLX.mat') % [mol/m2/sec]
load('SI_Fract.mat')
cd(my_pwd)
flx_vent = -mean(DICCFLX,2)*yr*Mc;

alpha_chan = squeeze(MonthlyMean(SI_Fract));
alpha_chan_DJF = (alpha_chan(:,1) + alpha_chan(:,2) + alpha_chan(:,12))/3;
alpha_chan_MAM = mean(alpha_chan(:,3:5),2);
alpha_chan_JJA = mean(alpha_chan(:,6:8),2);
alpha_chan_SON = mean(alpha_chan(:,9:11),2);
alpha_chan_AM = mean(alpha_chan,2);

% light only 
cd('../ctrl_par/sim_1')
%load('DICPCO2.mat')
load('DICCFLX.mat') % [mol/m2/sec]
cd(my_pwd)
flx_par = -mean(DICCFLX,2)*yr*Mc;


% Both 
cd('../ctrl_both/sim_1')
%load('DICPCO2.mat')
load('DICCFLX.mat') % [mol/m2/sec]
cd(my_pwd)
flx_both = -mean(DICCFLX,2)*yr*Mc;

% No ice
cd('../par0_ice0/sim_1')
%load('DICPCO2.mat')
load('DICCFLX.mat') % [mol/m2/sec]
cd(my_pwd)
flx_no_ice = -mean(DICCFLX,2)*yr*Mc;


%% Plotting
figure
width = 440;
height = 280;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

nz_ice = 10;
z_ice = 1:nz_ice;
ax_x = 47;
ax_y = 37;
ax_xLength = 260;
ax_yLength = 230;


% Bar properties
xpos = ax_x+1;
bar_yLength = 5;
yoffs_par = 18;

xmin = -68.5;
xmax = -41.5;
ymin = -50;
ymax = 40;

clrs = [237, 141, 14; 12, 160, 22; 0,0,0]/255;

set(gca,'units','points','position',[ax_x,ax_y,ax_xLength,ax_yLength])
plt1 = plot(yc_lat,flx_vent,'color',clrs(1,:),'linewidth',2.5,'DisplayName','DJF: Summer');
hold on
plt2 = plot(yc_lat,flx_par,'color',clrs(2,:),'linewidth',2.5,'DisplayName','MAM: Fall');
plt3 = plot(yc_lat,flx_both,'color',clrs(3,:),'linewidth',2.5,'DisplayName','JJA: Winter');
plt4 = plot(yc_lat,flx_no_ice,'k:','linewidth',3);
plot(yc_lat,0*yc_lat,'k--')
xlabel('Latitude (°)')
ylabel('CO_2 flux [gCm^{-2}yr^{-1}]')
% lgd= legend([plt1 plt2 plt3],'Capping','Light Attenuation','Capping + Light Atten.');
lgd= legend([plt1 plt2 plt3 plt4],'Capping','Light Attenuation','Capping + Light Atten.','No ice');
lgd.FontSize = 10;
xlim([xmin xmax])
ylim([ymin ymax])


% myxticks = ((-65:5:-45)-xmin)/(xmax-xmin);
% set(gca,'xtick',myxticks)
% set(gca,'xticklabel',{'-65','-60','-55','-50','-45',})
set(gca,'ytick',-20:10:ymax)
% set(gca,'yticklabel',[])
xlabel('Latitude (°)')
set(gca,'units','points','position',[ax_x,ax_y,ax_xLength,ax_yLength])

ax(15)=axes;
fld = repmat(PAR_SON,1,nz_ice); 
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 80; cint = 0.5;
colormap(ax(15),par_map(cmin,cmax,cint))
caxis([cmin cmax])
ypos_par = ax_y+8;
set(gca,'units','points','position',[xpos,ypos_par,ax_xLength,bar_yLength])
xlim([xmin xmax])
set(gca,'xtick',[])
set(gca,'xticklabel',[])

    % SON 
ax(11)=axes;
ypos_ice = ypos_par+bar_yLength+2;
fld = repmat(alpha_chan_SON,1,nz_ice); fld(fld<0.05) = NaN; 
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(11),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_ice,ax_xLength,bar_yLength])
xlim([xmin xmax])
annotation('textbox','String','SON','units','points','linestyle','none','Position',[ax_x+130,ypos_ice + 6,ax_xLength,bar_yLength],'FontSize',9)
set(gca,'xtick',[])
set(gca,'xticklabel',[])

    % JJA 
ax(16)=axes;
ypos_par = ypos_par + yoffs_par;
fld = repmat(PAR_JJA,1,nz_ice); 
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 80; cint = 0.5;
colormap(ax(16),par_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_par,ax_xLength,bar_yLength])
xlim([xmin xmax])
set(gca,'xtick',[])
set(gca,'xticklabel',[])  
    
    
ax(12)=axes;
ypos_ice = ypos_par+bar_yLength+2;
fld = repmat(alpha_chan_JJA,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(12),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_ice,ax_xLength,bar_yLength])
xlim([xmin xmax])
annotation('textbox','String','JJA','units','points','linestyle','none','Position',[ax_x+130,ypos_ice+5,ax_xLength,bar_yLength],'FontSize',9)
set(gca,'xtick',[])
set(gca,'xticklabel',[])

    % MAM 


ax(16)=axes;
ypos_par = ypos_par + yoffs_par;
fld = repmat(PAR_MAM,1,nz_ice); 
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 80; cint = 0.5;
colormap(ax(16),par_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_par,ax_xLength,bar_yLength])
xlim([xmin xmax])
set(gca,'xtick',[])
set(gca,'xticklabel',[])

ax(13)=axes;
ypos_ice = ypos_par+bar_yLength+2;
fld = repmat(alpha_chan_MAM,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(13),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_ice,ax_xLength,bar_yLength])
xlim([xmin xmax])
annotation('textbox','String','MAM','units','points','linestyle','none','Position',[xpos+75,ypos_ice + 6,ax_xLength,bar_yLength],'FontSize',9)
set(gca,'xtick',[])
set(gca,'xticklabel',[])



    % DJF 
    
ax(17)=axes;
ypos_par = ypos_par + yoffs_par;
fld = repmat(PAR_DJF,1,nz_ice); 
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 80; cint = 0.5;
colormap(ax(17),par_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_par,ax_xLength,bar_yLength])
xlim([xmin xmax])
set(gca,'xtick',[])
set(gca,'xticklabel',[])  
hcb_par = colorbar;
    
ax(14)=axes;
ypos_ice = ypos_par+bar_yLength+2;
fld = repmat(alpha_chan_DJF,1,nz_ice); fld(fld<0.05) = NaN;
[C1,h1] = contourf(yc_lat,z_ice,fld');
set(h1,'LineColor','none')
set(gca, 'visible', 'off')
cmin = 0; cmax = 1; cint = 0.05;
colormap(ax(14),ice_map(cmin,cmax,cint))
caxis([cmin cmax])
set(gca,'units','points','position',[xpos,ypos_ice,ax_xLength,bar_yLength])
xlim([xmin xmax])
annotation('textbox','String','DJF','units','points','linestyle','none','Position',[xpos+90,ypos_ice + 6,ax_xLength,bar_yLength],'FontSize',9)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
hcb_ice = colorbar;

% ice colorbar
% ax(20)=axes;
xpos_cbar = ax_x+ax_xLength+15;
set(hcb_ice,'units','points','position',[xpos_cbar, ax_y+20, 10, 180])
ylabel(hcb_ice,'Ice fraction','FontSize',10)

% par colorbar
% ax(21)=axes;
xpos_cbar = xpos_cbar+55;
set(hcb_par,'units','points','position',[xpos_cbar, ax_y+20, 10, 180])
ylabel(hcb_par,'PAR (Wm^{-2})','FontSize',10)

annotation('textbox','String','[+ve outgassing]','units','points','linestyle','none','Position',[ax_x+165,ax_y-70,ax_xLength,ax_yLength],'FontSize',10.5)

 


% % Testing PAR stuff
% 
%     % DJF 
% figure
% fld = repmat(PAR_DJF,1,nz_ice); 
% [C1,h1] = contourf(yc_lat,z_ice,fld');
% set(h1,'LineColor','none')
% set(gca, 'visible', 'off')
% cmin = 0; cmax = 80; cint = 0.5;
% colormap(par_map(cmin,cmax,cint))
% caxis([cmin cmax])
% % set(gca,'units','points','position',[ax_xLow,ax_yLow + 3*yoffs,ax_xLength,ax_yLength])
% xlim([xmin xmax])
% % annotation('textbox','String','DJF','units','points','linestyle','none','Position',[ax_xLow+75,ax_yLow + 3*yoffs + 3,ax_xLength,ax_yLength])
% % colorbar
% % set(gca,'xtick',[])
% % set(gca,'xticklabel',[])

% 
%     % MAM 
% figure
% ax(14)=axes;
% fld = repmat(PAR_MAM,1,nz_ice); 
% [C1,h1] = contourf(yc_lat,z_ice,fld');
% set(h1,'LineColor','none')
% set(gca, 'visible', 'off')
% cmin = 0; cmax = 80; cint = 0.5;
% colormap(ax(14),par_map(cmin,cmax,cint))
% caxis([cmin cmax])
% set(gca,'units','points','position',[ax_xLow,ax_yLow + 3*yoffs,ax_xLength,ax_yLength])
% xlim([xmin xmax])
% annotation('textbox','String','MAM','units','points','linestyle','none','Position',[ax_xLow+75,ax_yLow + 3*yoffs + 3,ax_xLength,ax_yLength])
% colorbar
% 
% 
%     % JJA 
% figure
% ax(14)=axes;
% fld = repmat(PAR_JJA,1,nz_ice); 
% [C1,h1] = contourf(yc_lat,z_ice,fld');
% set(h1,'LineColor','none')
% set(gca, 'visible', 'off')
% cmin = 0; cmax = 80; cint = 0.5;
% colormap(ax(14),par_map(cmin,cmax,cint))
% caxis([cmin cmax])
% set(gca,'units','points','position',[ax_xLow,ax_yLow + 3*yoffs,ax_xLength,ax_yLength])
% xlim([xmin xmax])
% annotation('textbox','String','JJA','units','points','linestyle','none','Position',[ax_xLow+75,ax_yLow + 3*yoffs + 3,ax_xLength,ax_yLength])
% colorbar
% 
% 
%     % SON 
% figure
% ax(14)=axes;
% fld = repmat(PAR_SON,1,nz_ice); 
% [C1,h1] = contourf(yc_lat,z_ice,fld');
% set(h1,'LineColor','none')
% set(gca, 'visible', 'off')
% cmin = 0; cmax = 80; cint = 0.5;
% colormap(ax(14),par_map(cmin,cmax,cint))
% caxis([cmin cmax])
% set(gca,'units','points','position',[ax_xLow,ax_yLow + 3*yoffs,ax_xLength,ax_yLength])
% xlim([xmin xmax])
% annotation('textbox','String','SON','units','points','linestyle','none','Position',[ax_xLow+75,ax_yLow + 3*yoffs + 3,ax_xLength,ax_yLength])
% colorbar



% figure
% plot(yc_lat,PAR_DJF)
% hold on
% plot(yc_lat,PAR_MAM)
% plot(yc_lat,PAR_JJA)
% plot(yc_lat,PAR_SON)
% legend('DJF','MAM','JJA','SON')
% 
% figure
% contourf(1:12,yc_lat,PAR)
% colorbar




% Functions 

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

function [ map ] = par_map(cmin,cmax,cint)
    n=floor(abs(cmax-cmin)/cint)/3;
    rgb_st = [250, 248, 217]/255;
    rgb_mid1 = [247, 210, 114]/255;
    rgb_mid2 = [245, 147, 20]/255;
    rgb_end = [79, 37, 9]/255;
    
    r = [linspace(rgb_st(1),rgb_mid1(1),n) linspace(rgb_mid1(1),rgb_mid2(1),n) linspace(rgb_mid2(1),rgb_end(1),n)];
    g = [linspace(rgb_st(2),rgb_mid1(2),n) linspace(rgb_mid1(2),rgb_mid2(2),n) linspace(rgb_mid2(2),rgb_end(2),n)];
    b = [linspace(rgb_st(3),rgb_mid1(3),n) linspace(rgb_mid1(3),rgb_mid2(3),n) linspace(rgb_mid2(3),rgb_end(3),n)];
    map = [r(:), g(:), b(:)];
end



clear all; clc; close all;

my_pwd = pwd();

% Constants
rhoO = 1030;
g = 9.81;
yr = 360*3600*24;

% Loading libraries
addpath('C:\Users\Mukund\Dropbox (MIT)\MIT\research\work\glacial_interglacial\2d_channel\library');
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library');
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library/seawater_ver3_3.1')

% Grid stuff
load('grid2d.mat')
% msk = DIC(:,:,end);% Better to use HFac for the mask. But shouldn't be an issue really
% msk(msk ~= 0) = 1;
DXG = grid2d.DXG;
XC = grid2d.XC;
YC = grid2d.YC;
RC = grid2d.RC;
DRF = grid2d.DRF;
DYC = grid2d.DYC;
ylat_min = -69.5; % I am not sure about this value
yc = mean(YC,1);
YC_lat = ylat_min + YC/(111*10^3);
yc_lat = ylat_min + yc/(111*10^3);
ny = length(yc);
dxg = mean(DXG(:));
nz = length(RC);

% Loading data
cd('../ctrl_both/sim_1')
load('THETA.mat')
load('SALT.mat')
load('SI_Fract.mat')
load('GM_PsiY.mat')
load('VVELMASS.mat')
load('DIC.mat')
cd(my_pwd)

% cd('../ctrl_both/sim_laud')
% surfDiag = rdmnc('surfDiag*');
% MXLDEPTH = squeeze(surfDiag.MXLDEPTH);
% cd(my_pwd)

% Taking time means
SALT = mean(SALT,3);
THETA = mean(THETA,3);
DIC = mean(DIC,3);
GM_PsiY = mean(GM_PsiY,3);
VVELMASS = mean(VVELMASS,3);
% MXLDEPTH = mean(MXLDEPTH,2);
THETA(THETA==0) = NaN;

% Getting psi
vvel = VVELMASS;
gm_psi = zeros(ny,nz+1);
gm_psi(:,1:nz) = GM_PsiY;
vvel_bol = zeros(ny,nz);

for k=1:nz
   vvel_bol(:,k) = ( gm_psi(:,k+1) - gm_psi(:,k) ) / DRF(k);
end

varE = zeros(ny,nz+1);
varB = zeros(ny,nz+1);

for k=nz:-1:1
    varE(:,k) = varE(:,k+1) - DRF(k)*DXG'.*vvel(:,k);
    varB(:,k) = varB(:,k+1) - DRF(k)*DXG'.*vvel_bol(:,k);
end
psiE(:,:) = varE(:,1:nz)*dxg;
psiB(:,:) = varB(:,1:nz)*dxg;
psiR = psiB + psiE; % Calculating residual

% Calculating isopycnals
Patm = 10.13; % Atmospheric pressure [dbar]
Pres = -RC*rhoO*g*1e-4 + Patm; % Pressure [dbar]
Pref_sigma = 0*rhoO*g*1e-4 + Patm; % Pressure [dbar]
% Calculating in situ temperature

for j=1:length(YC)
    for k=1:length(RC)
        T(j,k) = sw_temp(SALT(j,k),THETA(j,k),Pres(k),Patm);
        pden(j,k) = sw_pden(SALT(j,k),T(j,k),Pres(k),Pref_sigma);
    end
end
sigma = pden - 1000;

%% Ice edge calculations

SI_Fract = SI_Fract(:,end-11:end);
SI_min_idx = 61;
SI_max_idx = 153;
SI_min = zeros(ny,1);
SI_max = zeros(ny,1);

SI_min(1:SI_min_idx) = 1; SI_min(SI_min<1) = NaN; 
SI_max(1:SI_max_idx) = 1; SI_max(SI_max<1) = NaN; 


%% Plotting 

fig = figure;
fig.Units = 'points';
x0=50;
y0=50;
width=560;
height=525;
fnt_size=11;
fnt_name = 'Times New Roman';
ln_width = 0.5;
ax_xLength = 400;
ax_yLength = 200;
ax_xLow = 55;
ax_yLow = 45;
ax_yHigh = ax_yLow + ax_yLength + 40;
ax_xHigh = ax_xLow + ax_xLength + 45;

nz_ice = 10;
z_ice = 1:nz_ice;

set(gcf,'units','points','position',[x0,y0,width,height])
set(gcf,'color','white')

% DIC 
ax(1)=axes;
fld = DIC'; fld(fld==0) = NaN;
[C1,h1] = contourf(YC_lat,RC/1000,fld);
set(h1,'LineColor','none')
% cmin = min(Z1(:)); cmax = max(Z1(:)); cint = (max(Z1(:)) - min(Z1(:)))/10;
cmin = 2.05; cmax = 2.35; cint = 0.005;
colormap(ax(1),bmap1(cmin,cmax,cint))
hcb = colorbar;
% set(hcb1,'YTick',cmin:cint:cmax)
caxis([cmin cmax])
set(gca,'units','points','position',[ax_xLow,ax_yLow,ax_xLength,ax_yLength])
xlabel('Latitude (°)','FontSize',14)
xlabel(hcb,'DIC [mol C/m]','FontSize',12)
ylabel('Depth (km)','FontSize',14)
% xlim([0 max(YC_lat )])

% THETA  
ax(2)=axes;
fld=THETA'; fld(fld==0) = NaN;
[C,h] = contour(YC_lat,RC,fld,'LineWidth',1.5);
axis(ax(2),'off')
set(gca,'units','points','position',[ax_xLow,ax_yLow,ax_xLength,ax_yLength])
colormap(ax(2),bmap3(cmin,cmax,cint))
v = [-2 0 2 4 10 16];
clabel(C,h,v,'Color',[237,107,14]/256)

% psiR 
ax(4)=axes;
cmin = -80; cint = 0.5; cmax = 80;
fld=1e-6*psiR'; fld(fld==0) = NaN;
[C,h] = contourf(YC_lat,RC/1000,fld,'LineWidth',1.5);
set(h,'LineColor','none')
set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength+45,ax_yLength])
colormap(ax(4),bluepalered(cmin,cmax,cint))
% set(gca,'xticklabel',{[]})
hcb = colorbar;
ylabel('Depth (km)','FontSize',14)
xlabel(hcb,'Residual streamfunction [Sv]','FontSize',12)
caxis([cmin cmax])

ax(5)=axes;
fld=sigma'; fld(fld==0) = NaN;
vect = [26 26.5 27 27.5:0.05:28];
[C,h] = contour(YC_lat,RC,fld,vect,'LineWidth',1.5);
axis(ax(5),'off')
set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength,ax_yLength])
colormap(ax(5),bmap2(cmin,cmax,cint))
v = [26 26.5 27 27.5 28];
clabel(C,h,v)

% Isopycnals (sigma)
ax(5)=axes;
fld=sigma'; fld(fld==0) = NaN;
vect = [26 26.5 27 27.5:0.05:28];
[C,h] = contour(YC_lat,RC,fld,vect,'LineWidth',1.5);
axis(ax(5),'off')
set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength,ax_yLength])
colormap(ax(5),bmap2(cmin,cmax,cint))
v = [26 26.5 27 27.5 28];
clabel(C,h,v)

% Maximum ice edge
ax2=axes;
fld = repmat(SI_max,1,nz_ice); 
h1 = pcolor(yc_lat,z_ice,fld');
colormap(ax2,0.7*ones(1,3));
set(h1, 'EdgeColor', 'none');
set(gca, 'visible', 'off')
set(gca,'units','points','position',[ax_xLow,ax_yLow+ax_yLength,ax_xLength,10])

% Minimum ice edge
ax3=axes;
fld = repmat(SI_min,1,nz_ice); 
h1 = pcolor(yc_lat,z_ice,fld');
colormap(ax3,0.3*ones(1,3));
set(h1, 'EdgeColor', 'none');
set(gca, 'visible', 'off')
set(gca,'units','points','position',[ax_xLow,ax_yLow+ax_yLength,ax_xLength,10])


% %% Zoomed in figure
% 
% zmin = -300;
% [ ~, zidx ] = min( abs(zmin-RC) );
% 
% figure
% fig.Units = 'points';
% x0=50;
% y0=50;
% width=560;
% height=300;
% fnt_size=11;
% fnt_name = 'Times New Roman';
% ln_width = 0.5;
% ax_xLength = 400;
% ax_yLength = 200;
% ax_xLow = 55;
% ax_yLow = 45;
% ax_yHigh = ax_yLow + ax_yLength + 20;
% ax_yHigh2 = ax_yHigh + ax_yLength + 20;
% ax_xHigh = ax_xLow + ax_xLength + 45;
% 
% set(gcf,'units','points','position',[x0,y0,width,height])
% set(gcf,'color','white')
% 
% % psiR 
% ax(4)=axes;
% cmin = -80; cint = 1; cmax = 80;
% fld=1e-6*psiR; fld(fld==0) = NaN;
% [C,h] = contourf(YC_lat,RC(1:zidx),fld(:,1:zidx)','LineWidth',1.5);
% set(h,'LineColor','none')
% % set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength+45,ax_yLength])
% colormap(ax(4),bluepalered(cmin,cmax,cint))
% % set(gca,'xticklabel',{[]})
% hcb = colorbar;
% ylabel('Depth (m)')
% xlabel('Latitude (°)')
% xlabel(hcb,'Residual streamfunction [Sv]')
% caxis([cmin cmax])

% % Isopycnals (sigma)
% ax(5)=axes;
% fld=sigma; fld(fld==0) = NaN;
% vect = [26 26.5 27 27.5:0.05:28];
% [C,h] = contour(YC_lat,RC(1:zidx),fld(:,1:zidx)',vect,'LineWidth',1.5);
% axis(ax(5),'off')
% set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength,ax_yLength])
% colormap(ax(5),bmap2(cmin,cmax,cint))
% v = [26 26.5 27 27.5 28];
% clabel(C,h,v)

% % Mixed layer depth
% ax(6)=axes;
% plot(yc_lat,-MXLDEPTH/1000,'r','linewidth',2)
% axis(ax(6),'off')
% set(gca,'units','points','position',[ax_xLow,ax_yHigh,ax_xLength,ax_yLength])
% ylim([0 4])

%%
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
    r = linspace(0,0,n);
    g = linspace(0,0,n);
    b = linspace(0,0,n);
    map = [r(:), g(:), b(:)];
end

function map = bmap3(cmin,cmax,cint)
    n=floor(abs( (cmax-cmin)  /cint));
    
    r1= 237/256;
    g1= 107/256;
    b1 = 14/256;
    
    r = linspace(r1,r1,n);
    g = linspace(g1,g1,n);
    b = linspace(b1,b1,n);
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

function [ map ] = temp_map(cmin,cmax,cint)
%     n=floor(abs(cmin/cint));
%     m=floor(abs(cmax/cint));

    n = floor(abs(cmax-cmin)/cint/2);
    
    st1 = [24, 80, 186]/255;
    en1 = [200, 216, 247]/255;
    st2 = [242, 220, 208]/255;
    en2 = [193, 11, 14]/255;
    
    r = [linspace(st1(1),en1(1),n) linspace(st2(1),en2(1),n)];
    g = [linspace(st1(2),en1(2),n) linspace(st2(2),en2(2),n)];
    b = [linspace(st1(3),en1(3),n) linspace(st2(3),en2(3),n)];
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


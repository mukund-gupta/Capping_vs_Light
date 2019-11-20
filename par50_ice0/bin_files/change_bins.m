clear all; clc; close all;

my_pwd = pwd();
 
% Adding library
addpath('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/library')

% Grid stuff
load('grid2d.mat')
DRC = grid2d.DRC;
DRF = grid2d.DRF;
DXG = grid2d.DXG;
DXY = grid2d.DYG;
XC = grid2d.XC;
YC = grid2d.YC;
RC = grid2d.RC;
RF = grid2d.RF;
yc = mean(YC,1);
y2 = 5:10:(3200-5);
y3 = 2:4:(3200-2);
nx = 300;
ny2 = length(y2);
ny3 = length(y3);
nz=50;

% Changing the boundary conditions

% Reading par file
fname = 'par_monthly.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
par = reshape(var,ny2,12); 

% Modifying par file
par_idx = 160;
par_ice_conc = 0.50;
par_new = par;
par_new(1:par_idx,:) = (1 - par_ice_conc)*par_new(1:par_idx,:);
fname = 'par_50.bin';
fileID = fopen(fname,'w');
fwrite(fileID,par_new,'real*8','ieee-be');


% Creating sea ice file
sea_ice = zeros(ny2,12);
fname = 'ice_0.bin';
fileID = fopen(fname,'w');
fwrite(fileID,sea_ice,'real*8','ieee-be');

% Saving sea ice file
par_ice = zeros(ny2,12);
par_ice(1:par_idx,:) = par_ice_conc*ones(par_idx,12);
save('BCs.mat','sea_ice','par_ice')

% Plotting
figure
set(gcf,'color','w')

subplot(2,1,1)
fld = par;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 
colorbar

subplot(2,1,2)
fld = par_new;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 
colorbar

figure

subplot(2,1,1)
plot(yc/1000,mean(par,2),'r','linewidth',1.5)
hold on
plot(yc/1000,mean(par_new,2),'b','linewidth',1.5)

subplot(2,1,2)
plot(yc/1000,mean(sea_ice,2),'b','linewidth',1.5)
hold on






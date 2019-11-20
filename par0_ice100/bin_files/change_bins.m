clear all; clc; close all;
 
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

% Creating sea ice file
sea_ice = zeros(ny2,12);
ice_idx = 160;
sea_ice(1:ice_idx,:) = ones(ice_idx,12);

fname = 'sea_ice_2G8.bin';
fileID = fopen(fname,'w');
fwrite(fileID,sea_ice,'real*8','ieee-be');

% Reading ice file
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
sea_ice_test = reshape(var,ny2,12); % Specific humidity [kg/kg I believe]

% Saving sea ice file
save('sea_ice.mat','sea_ice')

% Plotting
figure
set(gcf,'color','w')
fld = squeeze(mean(sea_ice_test,2));
plot(yc/1000,fld,'b','linewidth',1.5)
title('Glacial sea ice (for DIC only)')
xlabel('y (km)') 
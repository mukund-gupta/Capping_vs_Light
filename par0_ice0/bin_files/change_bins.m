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

% Checking the theta profile file
fname = 'rbcs_temp_25.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
THETA = reshape(var,ny2,nz);
THETA_prof = mean(THETA,1)';

    % ctrl
cd('/Users/guptam/Dropbox (MIT)/MIT/research/work/glacial_interglacial/2d_channel/ctrl/bin_files')
fname = 'rbcs_temp_25.bin';
my_h2=fopen(fname,'r','b');
var=fread(my_h2,'real*8');
fclose(my_h2);
THETAc = reshape(var,ny2,nz);
THETAc_prof = mean(THETAc,1)';
cd(my_pwd)

figure
plot(THETAc_prof,RC,'color','r','linewidth',1.5)
hold on
plot(THETAc_prof,RC,'color','b','linewidth',1.5)
xlabel('Temperature (C)')
ylabel('Depth (m)')

% Changing the boundary conditions

% Creating sea ice file
sea_ice = zeros(ny2,12);
ice_idx = 160;
sea_ice(1:ice_idx,:) = 0*ones(ice_idx,12);

fname = 'sea_ice_2G12.bin';
fileID = fopen(fname,'w');
fwrite(fileID,sea_ice,'real*8','ieee-be');

% Reading ice file
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
sea_ice_test = reshape(var,ny2,12); % Specific humidity [kg/kg I believe]

% Saving sea ice file
par_ice = zeros(ny2,12);
save('BCs.mat','sea_ice','par_ice')

% Plotting
figure
set(gcf,'color','w')
fld = squeeze(mean(sea_ice,3));
plot(yc/1000,fld,'b','linewidth',1.5)
title('Glacial sea ice (for DIC only)')
xlabel('y (km)') 
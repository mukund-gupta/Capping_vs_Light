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

% Reading PAR file
fname = 'par_monthly.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
par = reshape(var,ny2,12); 

% Reading dwnSw file
fname = 'core_dwnSw_sec30E_320.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
dwnSw = reshape(var,ny2,12); 



% Plotting
figure
set(gcf,'color','w')

subplot(2,1,1)
fld = par;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 

subplot(2,1,2)
fld = dwnSw;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 

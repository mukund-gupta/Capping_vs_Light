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

% Loading sea ice file
load('../../ctrl/sim_1/SI_Fract.mat')
alpha = MonthlyMean(SI_Fract);

% Changing the boundary conditions

% Reading par file
fname = 'par_monthly.bin';
my_h=fopen(fname,'r','b');
var=fread(my_h,'real*8');
fclose(my_h);
par = reshape(var,ny2,12); 

% Modifying par file
par_new = par.*(1-alpha);
fname = 'par_real.bin';
fileID = fopen(fname,'w');
fwrite(fileID,par_new,'real*8','ieee-be');


% Plotting
figure
set(gcf,'color','w')

subplot(2,1,1)
fld = par;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 

subplot(2,1,2)
fld = par_new;
contourf(1:12,yc/1000,fld,'b','linewidth',1.5)
title('par')
xlabel('y (km)') 

function v_out = MonthlyMean(v)
    % Calculates yearly means from monthly data
    % Assumes time is the last dimension of the input matrix
    % Assumes there is an exact number of years.

    sl = 12; % length of slice (num of months in a year)
    dim = length(size(v))+1; % dimension of input array
    
    m = size(v);
    m2 = [m(1:end-1) sl m(end)/sl];

    v2 = reshape(v,m2);
    v_out = squeeze(mean(v2,dim));

end



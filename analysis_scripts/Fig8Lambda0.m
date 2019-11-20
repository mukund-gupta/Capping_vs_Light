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

%% Model data

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


%% Ventilation and light (with L)

av_idx = 160;
idx_st = 1;

my_runs = {'par0_ice0','par25_ice25','par50_ice50','par75_ice75','par90_ice90','par100_ice100'};
alpha_L = [0 .25 .5 .75 0.9 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    Fint_L(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
end


%% Ventilation and light (with L/2)

av_idx = 80;
idx_st = 1;

my_runs = {'par0_ice0','length_sens/both_0.25','length_sens/both_0.5','length_sens/both_0.75','par100_ice100'};
alpha_L2 = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    Fint_L2(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
end


%% Analytical model

% Initial parameters
load('fit_params.mat')

a_plt = 0:0.01:1;

% FL1 = getFnorm(a_plt,ke,kb,V,H,2*L,Cin); FL1(end) = 0;
% lambda_1 = (ke+kb)*2*L/(V*H);
% 
% FL2 = getFnorm(a_plt,ke,kb,V,H,L/2,Cin); FL2(end) = 0;
% lambda_2 = (ke+kb)*L/(V*H)/2;
% 
% FL3 = getFnorm(a_plt,ke,kb,V,H,L/4,Cin); FL3(end) = 0;
% lambda_3 = (ke+kb)*L/(V*H)/4;

FL1 = getFtot(a_plt,ke,kb,V,H,2*L,Cin); FL1(end) = 0;
lambda_1 = (ke+kb)*L/(V*H);

FL2 = getFtot(a_plt,ke,kb,V,H,L/2,Cin); FL2(end) = 0;
lambda_2 = (ke+kb)*L/(V*H)/2;

FL3 = getFtot(a_plt,ke,kb,V,H,L/4,Cin); FL3(end) = 0;
lambda_3 = (ke+kb)*L/(V*H)/4;


%%
% Plotting

figure
width = 270;
height = 220;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);
mksz = 130;
clrs = [0 0 0; 19, 98, 235]/255;

leg = ['Channel: \lambda_0 = ',num2str(round(lambda_1))];
sc1 = scatter(alpha_L,Fint_L/Fint_L(1),mksz,clrs(1,:),'filled','Displayname',leg);
hold on

leg = ['Channel: \lambda_0 = ',num2str(round(lambda_2))];
sc2 = scatter(alpha_L2,Fint_L2/Fint_L2(1),mksz,clrs(2,:),'filled','Displayname',leg);

leg = ['Analytical: \lambda_0 = ',num2str(round(lambda_1))];
plt1 = plot(a_plt,FL1/FL1(1),'DisplayName',leg,'linewidth',1.5,'color',clrs(1,:));

leg = ['Analytical: \lambda_0 = ',num2str(round(lambda_2))];
plt2 = plot(a_plt,FL2/FL2(1),'DisplayName',leg,'linewidth',1.5,'color',clrs(2,:));

% leg = ['\lambda_0 = ',num2str(round(lambda_3))];
% plt3 = plot(a_plt,FL3/FL3(1),'DisplayName',leg,'linewidth',1.5,'color',clrs(3,:));

xlabel('Sea ice fraction')
ylabel('Normalized flux')
legend show


% figure
% sc1 = scatter(alpha_L,Fint_L,mksz,clrs(1,:),'filled','Displayname','L = 1650 km');
% hold on
% sc2 = scatter(alpha_L2,Fint_L2,mksz,clrs(2,:),'filled','Displayname','L = 850 km');
% 
% leg = ['\lambda_0 = ',num2str(round(lambda_1))];
% plt1 = plot(a_plt,FL1,'DisplayName',leg,'linewidth',1.5,'color',clrs(1,:));
% 
% leg = ['\lambda_0 = ',num2str(round(lambda_2))];
% plt2 = plot(a_plt,FL2,'DisplayName',leg,'linewidth',1.5,'color',clrs(2,:));
% 
% xlabel('Sea ice fraction')
% ylabel('Normalized flux')
% legend show

%%
% Functions

function F = getFbio(alpha,ke,kb,V,H,L,Cin)
    F = -ke*Cin*V*H ./ (ke + kb*(1-alpha)) .* ( exp(-L/V/H*(ke+kb*(1-alpha))) -1 );
end

function F = getFblock(alpha,ke,kb,V,H,L,Cin)
    F = -ke*(1-alpha)*Cin*V*H ./ (ke*(1-alpha) + kb) .* ( exp(-L/V/H*(ke*(1-alpha) + kb)) -1 );
end

function F = getFtot(alpha,ke,kb,V,H,L,Cin)
    F = -ke*(1-alpha)*Cin*V*H ./ (ke*(1-alpha) + kb*(1-alpha)) .* ( exp(-L/V/H*(ke*(1-alpha) + kb*(1-alpha))) -1 );
end

function F = getFnorm(alpha,ke,kb,V,H,L,Cin)
    lambda = (ke*(1-alpha) + kb*(1-alpha))*L/V/H;
    lambda_0 = (ke + kb)*L/V/H;
    F = (exp(-lambda) - 1)/(exp(-lambda_0) - 1);
end

function [dFdke,dFdkb] = getGradfe(alpha,ke,kb,V,H,L,Cin)
    f = @(ke_v,kb_v) -ke_v*(1-alpha)*Cin*V*H ./ (ke_v*(1-alpha) + kb) .* ( exp(-L/V/H*(ke_v*(1-alpha) + kb)) -1 );
    [ke_v,kb_v] = ainit(ke,kb,3);   
    z = f(ke_v,kb_v);                               
    dFdke = z{1,0};
    dFdkb = z{0,1};
end

function [dFdke,dFdkb] = getGradfb(alpha,ke,kb,V,H,L,Cin)
    f = @(ke_v,kb_v) -ke_v*Cin*V*H ./ (ke_v + kb_v*(1-alpha)) .* ( exp(-L/V/H*(ke_v+kb_v*(1-alpha))) -1 );
    [ke_v,kb_v] = ainit(ke,kb,3);   
    z = f(ke_v,kb_v);                               
    dFdke = z{1,0};
    dFdkb = z{0,1};
end
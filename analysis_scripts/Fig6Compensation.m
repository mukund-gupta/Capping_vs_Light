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

av_idx = length(yc)/2;
idx_st = 1;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation
alpha_vent = [0 .25 .5 .75 1];

my_runs = {'par0_ice0','par0_ice25','par0_ice50','par0_ice75fr','par0_ice100'};

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
   
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_vent(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_vent_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;
    
    
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light

my_runs = {'par0_ice0','par25_ice0','par50_ice0','par75fr_ice0','par100_ice0'};
alpha_light = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_light(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_light_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation and light

my_runs = {'par0_ice0','par25_ice25','par50_ice50','par75_ice75','par90_ice90','par100_ice100'};
alpha_both = [0 .25 .5 .75 .9 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    tmp = mean(DICCFLX(idx_st:av_idx,tt),2);
    Fint_both(i) = -trapz(yc(idx_st:av_idx),tmp)*yr*Mc;
    
    tmp = mean(DICCFLX(1:ny,tt),2);
    Fint_both_tot(i) = -trapz(yc(1:ny),tmp)*yr*Mc;
    
end



%% Fitting
xi = [0 .25 .5 .75 1];
alpha = xi;
yb_i = Fint_light;
ye_i = Fint_vent;
N = length(xi);

% Initial parameters
H = 100; % Depth in [m]
V = 0.02; % Meridional velocity [m/s]
L = 1650*10^3; % length [m]
ke_ref = 0.95*H/yr; % Exchange coefficient (ensures 1yr e-folding timescale)
kb_ref = 2.2*ke_ref;
ke = ke_ref;
kb = kb_ref;

Cin_ref = yb_i(1)*(ke+kb) / (V*H*ke);
Cin = Cin_ref;
yb = getFbio(xi,ke,kb,V,H,L,Cin);
ye = getFblock(xi,ke,kb,V,H,L,Cin);

% Gradient descent on linear problem
G = sum((yb - yb_i).^2)/N + sum((ye - ye_i).^2)/N;
dG = 1e-25; % Learning rate

% figure
% width = 315;
% height = 250;
% set(gcf,'color','w')
% set(gcf,'units','points','position',[50,50,width,height]);
% mksz = 130;
% clrs = [48, 124, 244; 12, 160, 22; 255,0,0]/255;
% 
% a_plt = 0:0.01:1;
% Fe = getFblock(a_plt,ke,kb,V,H,L,Cin);
% Fb = getFbio(a_plt,ke,kb,V,H,L,Cin);
% Ftot = getFtot(a_plt,ke,kb,V,H,L,Cin); Ftot(end) = 0;
% 
% sc1 = scatter(xi,ye_i,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
% hold on
% plt1 = plot(a_plt,Fe,'color',clrs(1,:));
% sc2 = scatter(xi,yb_i,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
% pl2t2 = plot(a_plt,Fb,'color',clrs(2,:));
% sc3 = scatter(alpha_both,Fint_both,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
% plt3 = plot(a_plt,Ftot,'color',clrs(3,:));
% xlabel('Sea ice fraction')
% ylabel('Outgassing carbon flux (gCm^{-1}yr^{-1})')
% legend([sc1 sc2 sc3],{'Ventilation','Light','Ventilation + Light'},'Location','NorthWest')

% G
% for i=1:10
%     
%     [dfedke,dfedkb] = getGradfe(alpha,ke,kb,V,H,L,Cin); 
%     [dfbdke,dfbdkb] = getGradfb(alpha,ke,kb,V,H,L,Cin);
%     
%     dGdke = 2/N * sum(dfedke.*(ye - ye_i)) + 2/N * sum(dfbdke.*(yb - yb_i));
%     dGdkb = 2/N * sum(dfedkb.*(ye - ye_i)) + 2/N * sum(dfbdkb.*(yb - yb_i));
%     
%     ke = ke - dGdke*dG;
%     kb = kb - dGdkb*dG;
%     
%     ye = getFblock(xi,ke,kb,V,H,L,Cin);
%     yb = getFbio(xi,ke,kb,V,H,L,Cin);
% %     Cin = yb_i(1)*(ke + kb) / (V*H*ke);
%  
%     
%     G = sum((yb - yb_i).^2)/N + sum((ye - ye_i).^2)/N
%     
% end

tau_e = H/ke/yr
tau_b = H/kb/yr
tau_a = L/V/yr

%%
% Plotting

a_plt = 0:0.01:1;
Fe = getFblock(a_plt,ke,kb,V,H,L,Cin);
Fb = getFbio(a_plt,ke,kb,V,H,L,Cin);
Ftot = getFtot(a_plt,ke,kb,V,H,L,Cin); Ftot(end) = 0;

figure
width = 615;
height = 250;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);
mksz = 130;
clrs = [237, 141, 14; 12, 160, 22; 0,0,0]/255;

subplot(1,2,1)
sc1 = scatter(alpha_vent,ye_i,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
hold on
plt1 = plot(a_plt,Fe,'color',clrs(1,:));
sc2 = scatter(alpha_light,yb_i,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
pl2t2 = plot(a_plt,Fb,'color',clrs(2,:));
sc3 = scatter(alpha_both,Fint_both,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
plt3 = plot(a_plt,Ftot,'color',clrs(3,:));
xlabel('Sea ice fraction')
ylabel('Integrated flux (gC.m^{-1}.yr^{-1})')
ylim([0 6e7])
legend([sc1 sc2 sc3],{'Capping','Light Attenuation','Capping + Light Atten.'},'Location','NorthWest')
title('(a) Ice zone')


subplot(1,2,2)
sc1 = scatter(alpha_vent,Fint_vent_tot,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
hold on
sc2 = scatter(alpha_light,Fint_light_tot,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
sc3 = scatter(alpha_both,Fint_both_tot,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
xlabel('Sea ice fraction')
ylabel('Integrated flux (gC.m^{-1}.yr^{-1})')
% legend([sc1 sc2 sc3],{'Capping','Light Attenuation','Capping + Light Atten.'},'Location','NorthWest')
title('(b) Whole domain')



%% Sensitivity to V,L,ke and kb

a_plt = 0:0.01:1;
Fc = getFtot(a_plt,ke,kb,V,H,L,Cin); Fc(end) = 0;
Fv1 = getFtot(a_plt,ke,kb,V/2,H,L,Cin); Fv1(end) = 0;
Fv2 = getFtot(a_plt,ke,kb,2*V,H,L,Cin); Fv2(end) = 0;
FL1 = getFtot(a_plt,ke,kb,V,H,L/2,Cin); FL1(end) = 0;
FL2 = getFtot(a_plt,ke,kb,V,H,2*L,Cin); FL2(end) = 0;
Fke1 = getFtot(a_plt,ke/2,kb,V,H,L,Cin); Fke1(end) = 0;
Fke2 = getFtot(a_plt,2*ke,kb,V,H,L,Cin); Fke2(end) = 0;
Fkb1 = getFtot(a_plt,ke,kb/2,V,H,L,Cin); Fkb1(end) = 0;
Fkb2 = getFtot(a_plt,ke,2*kb,V,H,L,Cin); Fkb2(end) = 0;

figure
width = 600;
height = 500;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

ymin = 0;
ymax = 25;

subplot(2,2,1)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fv1,'DisplayName','V/2','linewidth',1.5)
plot(a_plt,Fv2,'DisplayName','2V','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to V')

subplot(2,2,2)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,FL1,'DisplayName','L/2','linewidth',1.5)
plot(a_plt,FL2,'DisplayName','2L','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to L')

subplot(2,2,3)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fke1,'DisplayName','k_e/2','linewidth',1.5)
plot(a_plt,Fke2,'DisplayName','2k_e','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_e')

subplot(2,2,4)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fkb1,'DisplayName','k_b/2','linewidth',1.5)
plot(a_plt,Fkb2,'DisplayName','2k_b','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_b')

supertitle(['Ventilation and Light',newline],'Fontsize',14)

% Normalized

Fc = getFnorm(a_plt,ke,kb,V,H,L,Cin); Fc(end) = 0;

Fv1 = getFnorm(a_plt,ke,kb,V/2,H,L,Cin); Fv1(end) = 0;
Fv2 = getFnorm(a_plt,ke,kb,2*V,H,L,Cin); Fv2(end) = 0;

FL1 = getFnorm(a_plt,ke,kb,V,H,L/2,Cin); FL1(end) = 0;
FL2 = getFnorm(a_plt,ke,kb,V,H,2*L,Cin); FL2(end) = 0;

Fke1 = getFnorm(a_plt,ke/2,kb,V,H,L,Cin); Fke1(end) = 0;
Fke2 = getFnorm(a_plt,2*ke,kb,V,H,L,Cin); Fke2(end) = 0;

Fkb1 = getFnorm(a_plt,ke,kb/2,V,H,L,Cin); Fkb1(end) = 0;
Fkb2 = getFnorm(a_plt,ke,2*kb,V,H,L,Cin); Fkb2(end) = 0;



figure
width = 600;
height = 500;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

ymin = 0;
ymax = 1;

subplot(2,2,1)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fv1,'DisplayName','V/2','linewidth',1.5)
plot(a_plt,Fv2,'DisplayName','2V','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (normalized)')
ylim([ymin ymax])
title('Sensitivity to V')

subplot(2,2,2)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,FL1,'DisplayName','L/2','linewidth',1.5)
plot(a_plt,FL2,'DisplayName','2L','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (normalized)')
ylim([ymin ymax])
title('Sensitivity to L')

subplot(2,2,3)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fke1,'DisplayName','k_e/2','linewidth',1.5)
plot(a_plt,Fke2,'DisplayName','2k_e','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (normalized)')
ylim([ymin ymax])
title('Sensitivity to k_e')

subplot(2,2,4)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fkb1,'DisplayName','k_b/2','linewidth',1.5)
plot(a_plt,Fkb2,'DisplayName','2k_b','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (normalized)')
ylim([ymin ymax])
title('Sensitivity to k_b')

supertitle(['Ventilation and Light (normalized)',newline],'Fontsize',14)


%% Sensitivity to V,L,ke and kb (Ventilation only)

a_plt = 0:0.01:1;
Fc = getFblock(a_plt,ke,kb,V,H,L,Cin); Fc(end) = 0;
Fv1 = getFblock(a_plt,ke,kb,V/2,H,L,Cin); Fv1(end) = 0;
Fv2 = getFblock(a_plt,ke,kb,2*V,H,L,Cin); Fv2(end) = 0;
FL1 = getFblock(a_plt,ke,kb,V,H,L/2,Cin); FL1(end) = 0;
FL2 = getFblock(a_plt,ke,kb,V,H,2*L,Cin); FL2(end) = 0;
Fke1 = getFblock(a_plt,ke/2,kb,V,H,L,Cin); Fke1(end) = 0;
Fke2 = getFblock(a_plt,2*ke,kb,V,H,L,Cin); Fke2(end) = 0;
Fkb1 = getFblock(a_plt,ke,kb/2,V,H,L,Cin); Fkb1(end) = 0;
Fkb2 = getFblock(a_plt,ke,2*kb,V,H,L,Cin); Fkb2(end) = 0;

figure
width = 600;
height = 500;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

ymin = 0;
ymax = 25;

subplot(2,2,1)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fv1,'DisplayName','V/2','linewidth',1.5)
plot(a_plt,Fv2,'DisplayName','2V','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to V')

subplot(2,2,2)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,FL1,'DisplayName','L/2','linewidth',1.5)
plot(a_plt,FL2,'DisplayName','2L','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to L')

subplot(2,2,3)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fke1,'DisplayName','k_e/2','linewidth',1.5)
plot(a_plt,Fke2,'DisplayName','2k_e','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_e')

subplot(2,2,4)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fkb1,'DisplayName','k_b/2','linewidth',1.5)
plot(a_plt,Fkb2,'DisplayName','2k_b','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_b')

supertitle(['Ventilation only',newline],'Fontsize',14)



%% Sensitivity to V,L,ke and kb (Light only)

a_plt = 0:0.01:1;
Fc = getFbio(a_plt,ke,kb,V,H,L,Cin); 
Fv1 = getFbio(a_plt,ke,kb,V/2,H,L,Cin); 
Fv2 = getFbio(a_plt,ke,kb,2*V,H,L,Cin);
FL1 = getFbio(a_plt,ke,kb,V,H,L/2,Cin); 
FL2 = getFbio(a_plt,ke,kb,V,H,2*L,Cin); 
Fke1 = getFbio(a_plt,ke/2,kb,V,H,L,Cin); 
Fke2 = getFbio(a_plt,2*ke,kb,V,H,L,Cin); 
Fkb1 = getFbio(a_plt,ke,kb/2,V,H,L,Cin);
Fkb2 = getFbio(a_plt,ke,2*kb,V,H,L,Cin); 

figure
width = 600;
height = 500;
set(gcf,'color','w')
set(gcf,'units','points','position',[50,50,width,height]);

ymin = 0;
ymax = 100;

subplot(2,2,1)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fv1,'DisplayName','V/2','linewidth',1.5)
plot(a_plt,Fv2,'DisplayName','2V','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to V')

subplot(2,2,2)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,FL1,'DisplayName','L/2','linewidth',1.5)
plot(a_plt,FL2,'DisplayName','2L','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to L')

subplot(2,2,3)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fke1,'DisplayName','k_e/2','linewidth',1.5)
plot(a_plt,Fke2,'DisplayName','2k_e','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_e')

subplot(2,2,4)
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Fkb1,'DisplayName','k_b/2','linewidth',1.5)
plot(a_plt,Fkb2,'DisplayName','2k_b','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (gCm^{-2}yr^{-1})')
ylim([ymin ymax])
title('Sensitivity to k_b')

supertitle(['Light only',newline],'Fontsize',14)




%% Sensitivity to lambda_0

Fc = getFnorm(a_plt,ke,kb,V,H,L,Cin); Fc(end) = 0;
Flambda1 = getFnorm(a_plt,ke,kb,V,H,L/2,Cin); Flambda1(end) = 0;
Flambda2 = getFnorm(a_plt,ke,kb,V,H,2*L,Cin); Flambda2(end) = 0;


ymin = 0;
ymax = 1;

figure
set(gcf,'color','w')
plot(a_plt,Fc,'DisplayName','ctrl','linewidth',2,'color',[0 0 0]/256)
hold on
plot(a_plt,Flambda1,'DisplayName','\lambda_0/2','linewidth',1.5)
plot(a_plt,Flambda2,'DisplayName','2\lambda_0','linewidth',1.5)
legend show
xlabel('Ice fraction')
ylabel('Flux (normalized)')
ylim([ymin ymax])
title('Ventilation & Light - Sensitivity to \lambda_0')



%% Integrating over the whole domain of the GCM (sensitivity) 

av_idx = length(yc);
idx_st = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation

my_runs = {'par0_ice0','par0_ice25','par0_ice50','par0_ice75fr','par0_ice100'};
alpha_vent = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    Fint_vent(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light

my_runs = {'par0_ice0','par25_ice0','par50_ice0','par75fr_ice0','par100_ice0'};
alpha_light = [0 .25 .5 .75 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    Fint_light(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ventilation and light

my_runs = {'par0_ice0','par25_ice25','par50_ice50','par75_ice75','par90_ice90','par100_ice100'};
alpha_both = [0 .25 .5 .75 .9 1];

for i=1:length(my_runs)
    % Loading data
    cd(['../',my_runs{i},'/sim_1'])
    %load('DICPCO2.mat')
    load('DICCFLX.mat') % [mol/m2/sec]
    cd(my_pwd)
    
    Fint_both(i) = -mean(mean(DICCFLX(idx_st:av_idx,tt),2),1)*yr*Mc;
    
end

figure
sc1 = scatter(alpha_vent,Fint_vent,mksz,clrs(1,:),'filled','Displayname','Blocking (No PAR)');
hold on
sc2 = scatter(alpha_light,Fint_light,mksz,clrs(2,:),'filled','Displayname','PAR (No Blocking)');
sc3 = scatter(alpha_both,Fint_both,mksz,clrs(3,:),'filled','Displayname','PAR and Blocking');
xlabel('Sea ice fraction')
ylabel('Outgassing carbon flux (gCm^{-2}yr^{-1})')
legend([sc1 sc2 sc3],{'Ventilation','Light','Ventilation + Light'},'Location','NorthWest')
title('Integrating over whole domain')
% The light/ventilation compensation still plays a big role. Perhaps not as
% big as over the ice region only because the bio part is partially limited
% over the northern part of the domain. 

%% Ventilation and light - meridional structure

alpha_b = 0.5;
alpha_e = 0.5;
dy = L/100;
y = 0:dy:L;

[Fc,~] = getFlat(alpha_e,alpha_b,ke,kb,V,H,L,Cin,y);

[Fv1,~] = getFlat(alpha_e,alpha_b,ke,kb,V/2,H,L,Cin,y);
[Fv2,~] = getFlat(alpha_e,alpha_b,ke,kb,2*V,H,L,Cin,y);

y1 = 0:(dy/2):L/2;
y2 = 0:(2*dy):2*L;
[FL1,~] = getFlat(alpha_e,alpha_b,ke,kb,V,H,L/2,Cin,y1);
[FL2,~] = getFlat(alpha_e,alpha_b,ke,kb,V,H,2*L,Cin,y2);

[Fke1,~] = getFlat(alpha_e,alpha_b,ke/2,kb,V,H,L,Cin,y);
[Fke2,~] = getFlat(alpha_e,alpha_b,2*ke,kb,V,H,L,Cin,y);

[Fkb1,~] = getFlat(alpha_e,alpha_b,ke,kb/2,V,H,L,Cin,y);
[Fkb2,~] = getFlat(alpha_e,alpha_b,ke,2*kb,V,H,L,Cin,y);

figure

clrs = [19, 140, 37; 0 0 0; 239, 97, 9]/255;

set(gcf,'color','w')

subplot(2,2,1)
plot(y/1000,Fv2*yr*Mc,'DisplayName','2V','linewidth',2.5,'color',clrs(3,:))
hold on
plot(y/1000,Fc*yr*Mc,'DisplayName','ctrl','linewidth',2.5,'color',clrs(2,:))
plot(y/1000,Fv1*yr*Mc,'DisplayName','V/2','linewidth',2.5,'color',clrs(1,:))

legend show
title('Sensitivity to V')
xlabel('Length (km)')

subplot(2,2,2)
plot(y2/1000,FL2,'DisplayName','2L','linewidth',2.5,'color',clrs(3,:))
hold on
plot(y/1000,Fc,'DisplayName','ctrl','linewidth',2.5,'color',clrs(2,:))
plot(y1/1000,FL1,'DisplayName','L/2','linewidth',2.5,'color',clrs(1,:))
legend show
title('Sensitivity to L')
xlabel('Length (km)')

subplot(2,2,3)
plot(y/1000,Fke2,'DisplayName','2k_e','linewidth',2.5,'color',clrs(3,:))
hold on
plot(y/1000,Fc,'DisplayName','ctrl','linewidth',2.5,'color',clrs(2,:))
plot(y/1000,Fke1,'DisplayName','k_e/2','linewidth',2.5,'color',clrs(1,:))
legend show
title('Sensitivity to k_e')
xlabel('Length (km)')

subplot(2,2,4)
plot(y/1000,Fkb2,'DisplayName','2k_b','linewidth',2.5,'color',clrs(3,:))
hold on
plot(y/1000,Fc,'DisplayName','ctrl','linewidth',2.5,'color',clrs(2,:))
plot(y/1000,Fkb1,'DisplayName','k_b/2','linewidth',2.5,'color',clrs(1,:))
legend show
title('Sensitivity to k_b')
xlabel('Length (km)')



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
    lambda_0 = (ke + kb)*L/V/H;
    lambda = lambda_0*(1-alpha);   
    F = (exp(-lambda) - 1)/(exp(-lambda_0) - 1);
end

function [F,C] = getFlat(alpha_e,alpha_b,ke,kb,V,H,L,Ci,y)
    lambda = (ke*(1-alpha_e) + kb*(1-alpha_b))*L/V/H;
    C = Ci*exp(-lambda*y/L);
    F = ke*(1-alpha_e)*C;
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

function supertitle(mytitle,varargin)

axes('Units','Normal');
h = title(mytitle,varargin{1:length(varargin)});
set(gca,'visible','off')
set(h,'visible','on')

end
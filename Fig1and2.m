%% Figure 1+2
% Generates figure 1 showing input curves
% Generates figure 2 showing validation against experimental data
% (both figures are generated with the peak height = 0.6 and with an
% ensemble model)
% Generates supplemental Figures: 
%       - S2 (tissue-level model parameters)
%       - S3 (network dynamics)
%       - S4 (alterations to peak height and correlations)
% Generates output figure for additional validation (not shown in paper)
% Generates plots of all input curves from the ensemble mdoel (not shown in
% paper)

% Last updated by ACZ 2.18.2020

clear
clc

peak = 0.6;
stdev = 0.05;


%% validation data

% collagen expression data (Deten 2001, Zimmerman 2001)
cData1_ave = [0.969, 1.607; 4.778, 3.75; 11.713, 16.964; 11.46, 12.857; 3.915, 6.964; 3.35, 6.25; 2.263, 3.571; 0.715, 2.679];
d1 = max(cData1_ave,[],1);
cData1_ave = cData1_ave./d1;
cData1_std = [0.834, 0.714; 1.354, 1.25; 1.041, 4.643; 1.667, 2.5; 1.52, 2.5; 0.937, 0.893; 0.625, 1.25; 0.521, 0.892];
cData1_std = cData1_std./d1;
cData1_time = [0.5, 1, 3, 6, 9, 12, 21, 55]./7;

cData2_ave = [0.22, 0.32; 0.756, 0.408; 0.76, 0.756; 0.368, 0.268; 0.228, 0.224];
d2 = max(cData2_ave,[],1);
cData2_ave = cData2_ave./d2;
cData2_std = [0.04, 0.08; 0.08, 0.032; 0.052, 0.152; 0.048, 0.048; 0.036, 0.032];
cData2_std = cData2_std./d2;
cData2_time = [1, 4, 7, 14, 28]./7;


% Collagen content data (Fomovsky AJP 2010)
time_rat = [7;14;21;42];
AF_rat = [9.0; 16.5; 21.0; 26.5];
AFsd_rat = [4.0; 8.0; 9.0; 7.0];


%% Simulations with peak height = 0.6

% generate input curves and plot curves for Figure 1
[InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve_12_19(peak, peak);

% extract the parameters from logic based ODE model
[params,y0] = fib617_params(peak);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};

% simulation using logic based ODE model (signal network)
outputIns = [29, 31, 83, 86, 88, 51, 96, 97, 98];
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
yI = real(interp1(t1,y1,tInSim));

% extract parameters from logic based model to feed into tissue-level model
Cmrna = sum(yI(:,[101,102]),2);
peakCol = max(Cmrna);

% simulation using tissue-level model
[c1_nom,days] = MISimODEplot(Cmrna,tInSim,peakCol); % plots figure used in S2

% prepare results of simulation for plotting
yOuts = yI(:,outputIns);
yIns = yI(:,[2, 27, 50, 55, 58, 17, 10, 37, 34]);
yCols = yI(168:end,[101,102])./max(yI(168:end,101));
tWeek = (tInSim(168:end) - 168)./168;

col1_nom = yI(168:end,101);
col3_nom = yI(168:end,102);

% plots similar to Figure 2 for peak height set to 0.6
figure
subplot(2,1,1)
yyaxis left
plot(tWeek,col1_nom);
ylim([0 0.3])
yyaxis right
hold on; errorbar(cData1_time,cData1_ave(:,1),cData1_std(:,1),'bo');
hold on; errorbar(cData2_time,cData2_ave(:,1),cData2_std(:,1),'b*');
ylim([0 1.4]);
xlim([0 9]);
legend('model','PMID: 11444923','PMID: 11557576');
subplot(2,1,2)
yyaxis left
plot(tWeek,col3_nom);
ylim([0 0.3])
yyaxis right
hold on; errorbar(cData1_time,cData1_ave(:,2),cData1_std(:,2),'bo');
hold on; errorbar(cData2_time,cData2_ave(:,2),cData2_std(:,2),'b*');
ylim([0 1.4]);
xlim([0 9]);
legend('model','PMID: 11444923','PMID: 11557576');


figure
subplot(1,2,1)
plot(tInSim,yOuts(:,1:5))
legend(speciesNames(outputIns(1:5)));
set(gca,'XTick',0:168:2329);
set(gca,'XTickLabel',-1:1:13,'fontsize',7);
xlabel('Time (weeks)');ylabel('Activity');
axis([0 2328 0 0.7]);
subplot(1,2,2)
plot(tInSim,yOuts(:,6:end))
legend(speciesNames(outputIns(6:end)));
set(gca,'XTick',0:168:2329);
set(gca,'XTickLabel',-1:1:13,'fontsize',7);
xlabel('Time (weeks)');ylabel('Activity');
axis([0 2328 0 0.7]);


% Outputs for additional validation (not shown in paper)
figure
plot(tInSim,yIns)
legend(speciesNames([2, 27, 50, 55, 58, 17, 10, 37, 34]));
set(gca,'XTick',0:168:2329);
set(gca,'XTickLabel',-1:1:13,'fontsize',7);
xlabel('Time (weeks)');ylabel('Activity');
axis([0 2328 0 1]);


%% Ensemble modeling

for i = 1:500
    disp(i)
    [InputCsim,tInSim,inputNode,resNorm,resNormConvert,weightR(i,:)] = InputCurve_12_19rand(peak, stdev);
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
    
    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
    yI2 = real(interp1(t1,y1,tInSim));
    Cmrna = sum(yI2(:,[101,102]),2);
    
    [c1,days] = MISimODE(Cmrna,tInSim,peakCol); %peakCol defined from standard simulation
    
    Net(:,:,i) = yI2;
    Carea(i,:) = c1;
    c1mrna(:,i) = yI2(168:end,101);
    c3mrna(:,i) = yI2(168:end,102);
    ins_ens(:,:,i) = InputCsim(:,:);
  
    
end


%% plot inputs from ensemble data
% experimental data
y1_tgfb = [-10, 10, 174, 149, 72];
t1_tgfb = [1, 24, 96, 672, 1344];
y1_il6 = [16.1, 38.4, 45.8, 10.8, 2.8, 2.5, 1.3, 1.3, 1.3, 2.5, 1.6];
t1_il6 = [3, 6, 12, 24, 72, 144, 216, 504, 672, 1344, 2016];
y1_il1 = [5.7, 25.6, 32.7, 15.7, 2.6, 1.9, 2.6, 3.1, 2.1, 2.1];
t1_il1 = [3, 6, 12, 24, 72, 144, 216, 672, 1344, 2016];
y1_tnfa = [0.16, 0.34, 0.31, 0.22, 0.19]./0.01; %normalized to the control not to sham
t1_tnfa = [24, 96, 264, 672, 960];
y1_ne = [1.24];
t1_ne = [1334];
y2_ne = [1.80];
t2_ne = [1176];
y1_et = [20, 36, 22, 5];
t1_et = [72, 120, 168, 1008];
y1_bnp = [14.6, 36.2, 22, 4.6];
t1_bnp = [72, 120, 168, 1008];
y1_ang = [1.31]./0.20;
t1_ang = [1334];
y2_ang = [21.8]./5.2;
t2_ang = [504];
y1_pdgf = [0.25, 2.76, 3.64, 4.52, 2.80, 2.11];
t1_pdgf = [24, 72, 168, 336, 672, 1008];


% ensemble inputs
mean_ins = mean(ins_ens,3);
ins_ens = permute(ins_ens,[3,2,1]);

set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
cmap = [0.9,0.9,0.9];


% AngII
subplot(3,3,1)
set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
yyaxis left
plot(tInSim, mean_ins(8,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,8),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_ang+168,y1_ang,'bo');hold on
plot(t2_ang+168,y2_ang,'b*');hold on
lgd=legend({'Idealized Curve','Hu 2014','Yamagishi 1993'},'Location','southeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('Angiotensin II');

%TGFb
subplot(3,3,2)
yyaxis left
plot(tInSim, mean_ins(1,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,1),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level')
yyaxis right
plot(t1_tgfb+168,y1_tgfb,'bo');hold on
lgd=legend({'Idealized Curve','Stavropoulou 2010'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('TGFB');

%IL6
subplot(3,3,3)
yyaxis left
plot(tInSim, mean_ins(2,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,2),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_il6+168,y1_il6,'bo');hold on
lgd=legend({'Idealized Curve','Deten 2002'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('IL6');

%IL1
subplot(3,3,4)
yyaxis left
plot(tInSim, mean_ins(3,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,3),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_il1+168,y1_il1,'bo');hold on
lgd=legend({'Idealized Curve','Deten 2002'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('IL1');

% TNFa
subplot(3,3,5)
yyaxis left
plot(tInSim, mean_ins(4,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,4),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_tnfa+168,y1_tnfa,'bo');hold on
lgd=legend({'Idealized Curve','Heba 2001'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('TNFa');

% NE
subplot(3,3,6)
yyaxis left
plot(tInSim, mean_ins(5,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,5),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right 
plot(t1_ne+168,y1_ne,'bo');hold on
plot(t2_ne+168,y2_ne,'b*');hold on
ylabel('Fold Change \it in vivo');
lgd=legend({'Idealized Curve','Veldhuisen 1995', 'Feng 2016'},'Location','southeast');
lgd.FontSize=14;
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('NE');

% ET1
subplot(3,3,7)
yyaxis left
plot(tInSim, mean_ins(6,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,6),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_et+168,y1_et,'bo');hold on
lgd=legend({'Idealized Curve','Loennechen 2001'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('ET1');

% BNP
subplot(3,3,8)
yyaxis left
plot(tInSim, mean_ins(7,:),'Color', 'k'); hold on 
b=plot(tInSim, ins_ens(:,:,7),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_bnp+168,y1_bnp,'bo');hold on
lgd=legend({'Idealized Curve','Loennechen 2001'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('NP');

% PDGF
subplot(3,3,9)
yyaxis left
plot(tInSim, mean_ins(9,:),'Color', 'k'); hold on 
b = plot(tInSim, ins_ens(:,:,9),'Color',cmap); hold on
uistack(b,'bottom')
ylim([0 peak*1.05]);
ylabel('Normalized Input Level');
yyaxis right
plot(t1_pdgf+168,y1_pdgf,'bo');hold on
lgd=legend({'Idealized Curve','Zhao 2011'},'Location','northeast');
lgd.FontSize=14;
ylabel('Fold Change \it in vivo');
set(gca,'XTick',0:336:2025);
set(gca,'XTickLabel',0:2:12,'fontsize',14);
xlabel('Time (Weeks)');
title('PDGF');

%% plots in Figure 2


% set up plot dependencies
tWeek = (tInSim(168:end) - 168)./168;
cmap = [0.2 0.2 0.2 0.05];

% find means of simulations
avec1 = mean(c1mrna,2);
avec3 = mean(c3mrna,2);
aveNet = mean(Net,3);
save('aveNet.mat','aveNet')



figure
subplot(2,1,1)
yyaxis left
a = plot(tWeek,c1mrna,'-','Color',cmap);hold on
b = plot(tWeek,avec1,'-','Color','k','LineWidth',2);hold on
c = plot(tWeek,col1_nom,'-','Color','m','LineWidth',2);
uistack(a,'bottom');
ylim([0 0.3])
yyaxis right
errorbar(cData1_time,cData1_ave(:,1),cData1_std(:,1),'bo');hold on; 
errorbar(cData2_time,cData2_ave(:,1),cData2_std(:,1),'b*');hold on; 
ylim([0 1.4]);
xlim([0 9]);
legend('ensemble','mean of ensemble','set peak','Deten 2001','Zimmerman 2001');
subplot(2,1,2)
yyaxis left
a = plot(tWeek,c3mrna,'-','Color',cmap);hold on
b = plot(tWeek,avec3,'-','Color','k','LineWidth',2);hold on
c = plot(tWeek,col3_nom,'-','Color','m','LineWidth',2);
uistack(a,'bottom');
ylim([0 0.3])
yyaxis right
hold on; errorbar(cData1_time,cData1_ave(:,2),cData1_std(:,2),'bo');
hold on; errorbar(cData2_time,cData2_ave(:,2),cData2_std(:,2),'b*');
ylim([0 1.4]);
xlim([0 9]);
legend('ensemble','mean of ensemble','set peak','Deten 2001','Zimmerman 2001');
title(strcat('Randomly Sampled Peak Height with normal range mean = ',num2str(peak),'st dev = ',num2str(stdev)));


tweek = days./7;

figure; 
errorbar(time_rat./7,AF_rat,AFsd_rat,'ko');
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
hold on; a = plot(tweek,mean(Carea,1),'-','Color','k','LineWidth',2);
hold on; b = plot(tweek,Carea,'-','Color',cmap);
hold on; c = plot(tweek,c1_nom,'--','Color','m','LineWidth',2);
uistack(b,'bottom');
title('nonlinear collagen model compared to data');
legend({'Fomovsky 2010','ensemble','mean of ensemble','set peak'},'Location','southeast');
title(strcat('Randomly Sampled Peak Height with normal range mean = ',num2str(peak),'st dev = ',num2str(stdev)));



inputLabel = {'TGFB'; 'IL6'; 'IL1'; 'TNFa'; 'NE'; 'ET1'; 'BNP'; 'Ang'; 'PDGF'};
figure
for i = 1:9
    subplot(3,3,i)
    plot(weightR(:,i),Carea(:,end),'o');
    xlabel('input weight'); ylabel('final area fraction');
    title(inputLabel(i))
end



% Option to generate supplemental figure with all activation for ensemble
% figure;
% imagesc(aveNet')
% colormap(flipud(bone))
% caxis([0 1])
% xlim([0,1682])
% ylabel('Network Node')
% xlabel('Time Post-MI (Weeks)')
% colorbar
% set(gca,'YTick',1:length(speciesNames),'FontSize',10);
% set(gca,'YTickLabel',[speciesNames],'FontSize',7);
% set(gca,'XTick',0:168:1512);
% set(gca,'XTickLabel',-1:1:8);
% disp('done')

% supplemental figure S3 for just 0.6
figure;
imagesc(yI')
colormap(flipud(bone))
caxis([0 1])
xlim([0,1682])
ylabel('Network Node')
xlabel('Time Post-MI (Weeks)')
colorbar
set(gca,'YTick',1:length(speciesNames),'FontSize',10);
set(gca,'YTickLabel',[speciesNames],'FontSize',7);
set(gca,'XTick',0:168:1512);
set(gca,'XTickLabel',-1:1:8);
disp('done')

% not shown figure of outputs from ensemble model
aveOuts = aveNet(:,outputIns);
figure
plot(tInSim,yOuts)
legend(speciesNames(outputIns));
set(gca,'XTick',0:168:2329);
set(gca,'XTickLabel',-1:1:13,'fontsize',7);
xlabel('Time (weeks)');ylabel('Activity');
axis([0 2328 0 0.11]);


%% Simulations to demonstrate how peak height affects simulations
% used in Figure S4

pk = [0.35:0.05:0.65];
numC = length(pk);



for i = 1:numC
    % extract the parameters
    [params,y0] = fib617_params(pk(i));
    [rpar,tau,ymax,speciesNames,KI]=params{:};

    disp(strcat('Simulation Number',num2str(i)))
    [InputCsim,tInSim,inputNode] = InputCurve_12_19NP(pk(i), pk(i));
    
    %running standard simulations
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
    
    % standard simulation
    outputIns = [101, 102, 103, 104, 95, 94, 96];
    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
    yI = real(interp1(t1,y1,tInSim));
    Cmrna = sum(yI(:,[101,102]),2);
    peakCol = max(Cmrna);
    [c1,days] = MISimODE(Cmrna,tInSim,peakCol);

    Carea(i,:) = c1;
    c1mrna(:,i) = yI(168:end,101);
    c3mrna(:,i) = yI(168:end,102);
    tWeek = (tInSim(168:end) - 168)./168;
end




cmap = colormap(bone(numC+1));

figure
title('Alter peak height')
subplot(2,1,1)
yyaxis left
for i = 1:numC
    plot(tWeek,c1mrna(:,i),'-','Color',cmap(i,:));hold on
end
ylim([0 0.3])
legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65');
subplot(2,1,2)
yyaxis left
for i = 1:numC
    plot(tWeek,c3mrna(:,i),'-','Color',cmap(i,:));hold on
end
ylim([0 0.3])
legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65');

figure; 
title('Alter peak height');
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
plot(tweek,Carea,'k');
 for i = 1:numC
    plot(tweek,Carea(i,:),'-','Color',cmap(i,:));hold on
 end
legend({'0.35','0.4','0.45','0.5','0.55','0.6','0.65'},'Location','southeast');


%% Simulations to demostrate how baseline levels vary
% used in figure S4

pk = [0:0.05:0.25];
numC = length(pk);



for i = 1:numC
    % extract the parameters
    [params,y0] = fib617_paramsMod(pk(i));
    [rpar,tau,ymax,speciesNames,KI]=params{:};

    disp(strcat('Simulation Number',num2str(i)))
    [InputCsim,tInSim,inputNode] = InputCurve_12_19base(pk(i), 0.6);
    
    %running standard simulations
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
    
    % standard simulation
    outputIns = [101, 102, 103, 104, 95, 94, 96];
    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
    yI = real(interp1(t1,y1,tInSim));
    Cmrna = sum(yI(:,[101,102]),2);
    peakCol = max(Cmrna);
    [c1,days] = MISimODE(Cmrna,tInSim,peakCol);

    Carea2(i,:) = c1;
    c1mrna2(:,i) = yI(:,101);
    c3mrna2(:,i) = yI(:,102);
    tWeek2 = (tInSim - 168)./168;
end




cmap = colormap(bone(numC+1));

figure
title('Alter baseline')
subplot(2,1,1)
yyaxis left
for i = 1:numC
    plot(tWeek2,c1mrna2(:,i),'-','Color',cmap(i,:));hold on
end
ylim([0 0.3])
xlim([-1 12])
legend('0','0.05','0.1','0.15','0.2','0.25','0.2','0.25');

subplot(2,1,2)
yyaxis left
for i = 1:numC
    plot(tWeek2,c3mrna2(:,i),'-','Color',cmap(i,:));hold on
end
ylim([0 0.3])
xlim([-1 12])
legend('0','0.05','0.1','0.15','0.2','0.25','0.2','0.25');



figure; 
title('Alter baseline')
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
plot(tweek,Carea2,'k');
 for i = 1:numC
    plot(tweek,Carea2(i,:),'-','Color',cmap(i,:));hold on
 end
 xlim([-1 12])
legend({'0','0.05','0.1','0.15','0.2','0.25','0.2','0.25','0.3','0.35'},'Location','southeast');


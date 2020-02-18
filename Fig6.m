%% %% figure 6 for dynamic MI paper, perturbation dosing effects on collagen area accumulation
clear
clc

max_values=[0.04,.1,1,10,25];
legendName = {'ymax=0.04','ymax=0.1','ymax=1 (Control)','ymax=10','ymax=25'};
pHeight = 0.5;

%% simulation with no alteration


% get default input curves and establish parameters
[InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve_12_19NP(pHeight, pHeight);
week2=tInSim./7./24;
week2=week2-1;
% extract the parameters
paramName = 'fib617_params';
eval(strcat('[params,y0] = ',paramName,';'));
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
numSpec = length(speciesNames);
% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));

yIcol=sum(yI(:,[101,102]),2);


%% collagen accumulation calculations control

days =  linspace(-7,89,900)';
fibro = [ones(1,71)*0.1,linspace(0.1,1,70),ones(1,759)]';
kgen = 0.142*24;         % area fraction / day
kdeg = 0.002*24;        % / day
rawCol1 = yIcol';
%rawCol1 = control_collagen_mRNA';
 tday = (tInSim -168)./24;
 week = days./7;

rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
figure;
plot(week,c1,'k')
xlim([-1,9])
ylim([0,40])
ylabel('Area Fraction %')
xlabel('Time (weeks)')


%% Run simulations under various node upregulations
% %% B1integrin simulations
% 
% %assign ymax values
% max_values=[.1,.5,1,5,10];
% yIcol=[]; %Summed CImRNA+CIIImRNA values
% B1_course=[]; %ACE timecourse
% for i = 1:length(max_values)
% [rpar,tau,ymax,speciesNames,KI]=params{:};
% ymax_new=ymax;
% %ymax B1int
% ymax_new(43)=max_values(i);
% 
% paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
% 
% options = [];
% [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
% yI = real(interp1(t1,y1,tInSim));
% 
% yIcol(:,i)=sum(yI(:,[101,102]),2);
% B1_course(:,i)=yI(:,43);
% end
% 
% %collagen calculations
% col_acc=[]; %accumulation
% 
% for i =1:length(yIcol(1,:))
% rawCol1 = yIcol(:,i)';
%  tday = (tInSim -168)./24;
%  week = days./7;
% rawCol2 = interp1(tday,rawCol1,days);
% % raw collagen is scaled to the peak of the standard simulation and the
% % minimum of the current simulation
% kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 
% 
% % MMP production
% % idealized curve based on post MI data in PMID 8531210
% kd0 = 0.2;            % baseline (minimum) value
% kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
% kd3 = 0.00432*24;       % time constant for fall (/day)
% tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
% kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
% kdt = [ones(70,1).*kd0;kdt];
% % run ODE
% [t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
% c1(c1 > 100) = 100;
% col_acc(:,i)=c1;
% end
% 
% 
% 
% 
% %% make plots for B1 ymax manipulation
% 
% figure;
% plot(week2,B1_course(:,1),':','color',[33,113,181]./256)
%  hold on
%  plot(week2,B1_course(:,2),'--','color',[107,174,214]./256)
%  hold on
%  plot(week2,B1_course(:,3),'k-')
%  hold on
% plot(week2,B1_course(:,4),':','color',[251,106,74]./256)
%  hold on
% plot(week2,B1_course(:,5),'--','color',[165,15,21]./256)
% legend({'ymax=0.1','ymax=0.5','ymax=1 (Control)','ymax=5','ymax=10'});
% title('B1int Expression Over Time')
% xlim([-1,9])
% ylabel('Node Expression')
% xlabel('Time (Weeks)')
% 
% 
% figure;
%  plot(week,col_acc(:,1),':','color',[33,113,181]./256)
%  hold on
%  plot(week,col_acc(:,2),'--','color',[107,174,214]./256)
%  hold on
%  plot(week,col_acc(:,3),'k-')
%  hold on
% plot(week,col_acc(:,4),':','color',[251,106,74]./256)
%  hold on
% plot(week,col_acc(:,5),'--','color',[165,15,21]./256)
% legend({'ymax=0.1','ymax=0.5','ymax=1 (Control)','ymax=5','ymax=10'});
% title('B1int effect on collagen accumulation')
% xlim([-1,9])
% ylabel('Area Fraction (%)')
% xlabel('Time (Weeks)')
% ylim([0,110])
% 
% 
% figure;
% % Collagen mRNA over time
% plot(week2,yIcol(:,1),':','color',[33,113,181]./256)
%  hold on
%  plot(week2,yIcol(:,2),'--','color',[107,174,214]./256)
%  hold on
%  plot(week2,yIcol(:,3),'k-')
%  hold on
% plot(week2,yIcol(:,4),':','color',[251,106,74]./256)
%  hold on
% plot(week2,yIcol(:,5),'--','color',[165,15,21]./256)
% legend({'ymax=0.1','ymax=0.5','ymax=1 (Control)','ymax=5','ymax=10'});
% title('B1int effect on collagen mRNA')
% xlim([-1,9])
% ylabel('Sum(CImRNA & CIIImRNA)')
% xlabel('Time (Weeks)')
% ylim([0,2.1])

%% smad7 simulations

%assign ymax values
yIcol=[]; %Summed CImRNA+CIIImRNA values
smad7_course=[]; %smad7 timecourse
for i = 1:length(max_values)
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_new=ymax;
%ymax smad7
ymax_new(30)=max_values(i);

paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));

yIcol(:,i)=sum(yI(:,[101,102]),2);
smad7_course(:,i)=yI(:,30);
end

%collagen calculations
col_acc=[]; %accumulation

for i =1:length(yIcol(1,:))
rawCol1 = yIcol(:,i)';
 tday = (tInSim -168)./24;
 week = days./7;
rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
col_acc(:,i)=c1;
end

% smad7 sensitivity analysis
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_up = ymax;
ymax_up(30) = 10;

for i = 1:length(speciesNames)
    ymax_new = ymax_up;
    ymax_new(i) = 0;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
    yI = real(interp1(t1,y1,tInSim));

    yIcol_ko(:,i)=sum(yI(:,[101,102]),2);
    
end

for j =1:length(yIcol_ko(1,:))
rawCol1 = yIcol_ko(:,j)';
 tday = (tInSim -168)./24;
 week = days./7;
rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
col_acc1(:,j)=c1(end);
end

fig=figure;
imagesc(real(col_acc1),[0 100]);

colormap(flipud(bone));
caxis([0 100]);
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',13);
xlabel('Species');
xtickangle(270);
colorbar('Location','eastoutside');



%% make plots for smad7 ymax manipulation

figure;
plot(week2,smad7_course(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,smad7_course(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,smad7_course(:,3),'k-')
 hold on
plot(week2,smad7_course(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,smad7_course(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad7 Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
 plot(week,col_acc(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week,col_acc(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week,col_acc(:,3),'k-')
 hold on
plot(week,col_acc(:,4),':','color',[251,106,74]./256)
 hold on
plot(week,col_acc(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad7 effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
% Collagen mRNA over time
plot(week2,yIcol(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,yIcol(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,yIcol(:,3),'k-')
 hold on
plot(week2,yIcol(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,yIcol(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad7 effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')

%% Il1 simulations

%assign ymax values
yIcol=[]; %Summed CImRNA+CIIImRNA values
Il1_course=[]; %Il1 timecourse
for i = 1:length(max_values)
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_new=ymax;
%ymax Il1
ymax_new(54)=max_values(i);

paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));

yIcol(:,i)=sum(yI(:,[101,102]),2);
Il1_course(:,i)=yI(:,54);
end

%collagen calculations
col_acc=[]; %accumulation

for i =1:length(yIcol(1,:))
rawCol1 = yIcol(:,i)';
 tday = (tInSim -168)./24;
 week = days./7;
rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
col_acc(:,i)=c1;
end




%% make plots for Il1 ymax manipulation

figure;
plot(week2,Il1_course(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,Il1_course(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,Il1_course(:,3),'k-')
 hold on
plot(week2,Il1_course(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,Il1_course(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
 plot(week,col_acc(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week,col_acc(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week,col_acc(:,3),'k-')
 hold on
plot(week,col_acc(:,4),':','color',[251,106,74]./256)
 hold on
plot(week,col_acc(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 effect on collagen accumulation')
xlim([-1,9])
ylim([0,110])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
% Collagen mRNA over time
plot(week2,yIcol(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,yIcol(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,yIcol(:,3),'k-')
 hold on
plot(week2,yIcol(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,yIcol(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')

%% PKG simulations

%assign ymax values
yIcol=[]; %Summed CImRNA+CIIImRNA values
PKG_course=[]; %PKG timecourse
for i = 1:length(max_values)
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_new=ymax;
%ymax PKG
ymax_new(40)=max_values(i);

paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));

yIcol(:,i)=sum(yI(:,[101,102]),2);
PKG_course(:,i)=yI(:,40);
end

%collagen calculations
col_acc=[]; %accumulation

for i =1:length(yIcol(1,:))
rawCol1 = yIcol(:,i)';
 tday = (tInSim -168)./24;
 week = days./7;
rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
col_acc(:,i)=c1;
end




%% make plots for PKG ymax manipulation

figure;
plot(week2,PKG_course(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,PKG_course(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,PKG_course(:,3),'k-')
 hold on
plot(week2,PKG_course(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,PKG_course(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('PKG Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
 plot(week,col_acc(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week,col_acc(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week,col_acc(:,3),'k-')
 hold on
plot(week,col_acc(:,4),':','color',[251,106,74]./256)
 hold on
plot(week,col_acc(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('PKG effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
% Collagen mRNA over time
plot(week2,yIcol(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,yIcol(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,yIcol(:,3),'k-')
 hold on
plot(week2,yIcol(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,yIcol(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('PKG effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')

%% NOX simulations

%assign ymax values
yIcol=[]; %Summed CImRNA+CIIImRNA values
ACE_course=[]; %NOX timecourse
for i = 1:length(max_values)
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_new=ymax;
%ymax NOX
ymax_new(7)=max_values(i);

paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));

yIcol(:,i)=sum(yI(:,[101,102]),2);
NOX_course(:,i)=yI(:,7);
end

%collagen calculations
col_acc=[]; %accumulation

for i =1:length(yIcol(1,:))
rawCol1 = yIcol(:,i)';
 tday = (tInSim -168)./24;
 week = days./7;
rawCol2 = interp1(tday,rawCol1,days);
% raw collagen is scaled to the peak of the standard simulation and the
% minimum of the current simulation
kgt = (1-0.0555).*(rawCol2 - min(rawCol2)) / (0.3173 - min(rawCol2)) + 0.0555; 

% MMP production
% idealized curve based on post MI data in PMID 8531210
kd0 = 0.2;            % baseline (minimum) value
kd2 = 0.0466*24/2;        % time constant for rise (/day) %JWH modified this 2x!
kd3 = 0.00432*24;       % time constant for fall (/day)
tpd = log(kd2/kd3)/(kd2-kd3); kd1 = (1-kd0)/(exp(-kd3*tpd)-exp(-kd2*tpd));
kdt = kd0 + kd1*(exp(-kd3*days(71:end))-exp(-kd2*days(71:end)));
kdt = [ones(70,1).*kd0;kdt];
% run ODE
[t,c1] = ode45(@nlincoll,days,3.0,[],kgen,kdeg,days,fibro,kgt,kdt);
c1(c1 > 100) = 100;
col_acc(:,i)=c1;
end




%% make plots for NOX ymax manipulation

figure;
plot(week2,NOX_course(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,NOX_course(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,NOX_course(:,3),'k-')
 hold on
plot(week2,NOX_course(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,NOX_course(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')



figure;
 plot(week,col_acc(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week,col_acc(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week,col_acc(:,3),'k-')
 hold on
plot(week,col_acc(:,4),':','color',[251,106,74]./256)
 hold on
plot(week,col_acc(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX effect on collagen accumulation')
xlim([-1,9])
ylim([0,110])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')

figure;
% Collagen mRNA over time
plot(week2,yIcol(:,1),':','color',[33,113,181]./256)
 hold on
 plot(week2,yIcol(:,2),'--','color',[107,174,214]./256)
 hold on
 plot(week2,yIcol(:,3),'k-')
 hold on
plot(week2,yIcol(:,4),':','color',[251,106,74]./256)
 hold on
plot(week2,yIcol(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')


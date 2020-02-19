%% Generates figure 6: highlighted regulators of collagen
% Generates figure S8: other highlighted regulators 
% also runs sensitivity analysis for schematic development
% last updated ACZ 2.18.2020

clear
clc


peak = 0.6;

nodeIndx = [30, 40, 7, 55, 43, 12, 29];
max_values=[0.1,0.5,1,5,10]; %
legendName = {'ymax=0.1','ymax=0.5','ymax=1 (Control)','ymax=5','ymax=10'};

%% control simulation


[InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve_12_19NP(peak, peak);

% extract the parameters
[params,y0] = fib617_params(peak);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};

options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
yI = real(interp1(t1,y1,tInSim));
Cmrna = sum(yI(:,[101,102]),2);
peakCol = max(Cmrna);
[c1_nom,days] = MISimODE(Cmrna,tInSim,peakCol);


week2=(tInSim-168)./168;
week = days./7;

%% smad7 simulations Figure 6a

smadMax = max_values./2; % must divide by 2 because default ymax of smad7 = 0.5
    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(1))=smadMax(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(1));
        Carea(:,j) = c1;
        
    end
    
    

% plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend({'ymax=0.05','ymax=0.1','ymax=0.5 (control)','ymax=2.5','ymax=5'});
title('smad7 Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend({'ymax=0.05','ymax=0.1','ymax=0.5 (control)','ymax=2.5','ymax=5'});
title('smad7 effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend({'ymax=0.05','ymax=0.1','ymax=0.5 (control)','ymax=2.5','ymax=5'});
title('smad7 effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')




% smad7 sensitivity analysis using all peak = set peak height
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_up = ymax;
ymax_up(nodeIndx(1)) = smadMax(end-1);
disp('smad7 sens')

for i = 1:length(speciesNames)
    ymax_new = ymax_up;
    ymax_new(i) = 0;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
    yI = real(interp1(t1,y1,tInSim));

    yIcol_ko(:,i)=sum(yI(:,[101,102]),2);
    
    [c_ko,days] = MISimODE(yIcol_ko(:,i),tInSim,peakCol);
    colAcc_ko(i) = c_ko(end);
    
end

figure
imagesc(real(colAcc_ko),[0 100]);
colormap(flipud(bone));
caxis([0 50]);
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',13);
xlabel('Species');
xtickangle(270);
colorbar('Location','eastoutside');
title('smad7 sensitivity');


%% PKG simulations Figure 6b

    
    for j = 1:length(max_values)
        disp(j)

        ymax_new=ymax;
        ymax_new(nodeIndx(2))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(2));
        Carea(:,j) = c1;



    end
    

 % plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('pkg Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('pkg effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('pkg effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')




% pkg sensitivity analysis using all peak = set peak height
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_up = ymax;
ymax_up(nodeIndx(2)) = max_values(end-1);
disp('pkg sens');
for i = 1:length(speciesNames)
    ymax_new = ymax_up;
    ymax_new(i) = 0;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
    yI = real(interp1(t1,y1,tInSim));

    yIcol_ko(:,i)=sum(yI(:,[101,102]),2);
    
    [c_ko,days] = MISimODE(yIcol_ko(:,i),tInSim,peakCol);
    colAcc_ko(i) = c_ko(end);
    
end

figure
imagesc(real(colAcc_ko),[0 100]);
colormap(flipud(bone));
caxis([0 50]);
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',13);
xlabel('Species');
xtickangle(270);
colorbar('Location','eastoutside');
title('pkg sensitivity');


%% NOX simulations Fig 6c

    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(3))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(3));
        Carea(:,j) = c1;



    end
    


% plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('NOX effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')


% nox sensitivity analysis using all peak = set peak height
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_up = ymax;
ymax_up(nodeIndx(3)) = max_values(end-1);


for i = 1:length(speciesNames)
    ymax_new = ymax_up;
    ymax_new(i) = 0;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
    yI = real(interp1(t1,y1,tInSim));

    yIcol_ko(:,i)=sum(yI(:,[101,102]),2);
    
    [c_ko,days] = MISimODE(yIcol_ko(:,i),tInSim,peakCol);
    colAcc_ko(i) = c_ko(end);
    
end

figure
imagesc(real(colAcc_ko),[0 100]);
colormap(flipud(bone));
caxis([0 100]);
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',13);
xlabel('Species');
xtickangle(270);
colorbar('Location','eastoutside');
title('NOX sensitivity');

%% inputIL1 simulations Fig6d
    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(4))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(4));
        Carea(:,j) = c1;



    end


% plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('IL1 effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')




% inputIL1 sensitivity analysis using all peak = set peak height
[rpar,tau,ymax,speciesNames,KI]=params{:};
ymax_up = ymax;
ymax_up(nodeIndx(4)) = max_values(end-1);
disp('IL1 sens')

for i = 1:length(speciesNames)
    ymax_new = ymax_up;
    ymax_new(i) = 0;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};

    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); %y0 loaded previously
    yI = real(interp1(t1,y1,tInSim));

    yIcol_ko(:,i)=sum(yI(:,[101,102]),2);
    
    [c_ko,days] = MISimODE(yIcol_ko(:,i),tInSim,peakCol);
    colAcc_ko(i) = c_ko(end);
    
end

figure
imagesc(real(colAcc_ko),[0 100]);
colormap(flipud(bone));
caxis([0 100]);
set(gca,'XTick',1:length(speciesNames));
set(gca,'XTickLabel',speciesNames,'fontsize',13);
xlabel('Species');
xtickangle(270);
colorbar('Location','eastoutside');
title('inputIL1 sensitivity');



%% sup: B1int simulations Figure S8a

    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(5))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(5));
        Carea(:,j) = c1;



    end
    



% plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('B1int Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('B1int effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('B1int effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')

%% sup: ET1R simulations Fig S8b
    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(6))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(6));
        Carea(:,j) = c1;



    end
    


% plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('ETAR Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('ETAR effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('ETAR effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')


%% sup: smad3 simulations FigS8c

    
    for j = 1:length(max_values)
        disp(j)
        % re-define smad7
        ymax_new=ymax;
        ymax_new(nodeIndx(7))=max_values(j);
        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};


        options = [];
        [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,paramsNew);
        yI = real(interp1(t1,y1,tInSim));
        Cmrna(:,j) = sum(yI(:,[101,102]),2);

        [c1,days] = MISimODE(Cmrna(:,j),tInSim,peakCol); %peakCol defined from standard simulation

        nodeAct(:,j) = yI(:,nodeIndx(7));
        Carea(:,j) = c1;



    end
    

   % plot results
figure;
plot(week2,nodeAct(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,nodeAct(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,nodeAct(:,3),'k-'); hold on;
plot(week2,nodeAct(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,nodeAct(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad3 Expression Over Time')
xlim([-1,9])
ylabel('Node Expression')
xlabel('Time (Weeks)')


figure;
plot(week,Carea(:,1),':','color',[33,113,181]./256);hold on;
plot(week,Carea(:,2),'--','color',[107,174,214]./256);hold on;
plot(week,Carea(:,3),'k-');hold on;
plot(week,Carea(:,4),':','color',[251,106,74]./256);hold on;
plot(week,Carea(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad3 effect on collagen accumulation')
xlim([-1,9])
ylabel('Area Fraction (%)')
xlabel('Time (Weeks)')


figure;
plot(week2,Cmrna(:,1),':','color',[33,113,181]./256);hold on;
plot(week2,Cmrna(:,2),'--','color',[107,174,214]./256);hold on;
plot(week2,Cmrna(:,3),'k-');hold on;
plot(week2,Cmrna(:,4),':','color',[251,106,74]./256);hold on;
plot(week2,Cmrna(:,5),'--','color',[165,15,21]./256)
legend(legendName);
title('smad3 effect on collagen mRNA')
xlim([-1,9])
ylabel('Sum(CImRNA & CIIImRNA)')
xlabel('Time (Weeks)')
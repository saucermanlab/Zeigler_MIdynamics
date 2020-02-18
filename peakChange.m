%% Two different ways to test peak height

% set up validation data
cData1_ave = [0.969, 1.607; 4.778, 3.75; 11.713, 16.964; 11.46, 12.857; 3.915, 6.964; 3.35, 6.25; 2.263, 3.571; 0.715, 2.679];
d1 = max(cData1_ave(:));
cData1_ave = cData1_ave./d1;
cData1_std = [0.834, 0.714; 1.354, 1.25; 1.041, 4.643; 1.667, 2.5; 1.52, 2.5; 0.937, 0.893; 0.625, 1.25; 0.521, 0.892];
cData1_std = cData1_std./d1;
cData1_time = [0.5, 1, 3, 6, 9, 12, 21, 55]./7;

cData2_ave = [0.22, 0.32; 0.756, 0.408; 0.76, 0.756; 0.368, 0.268; 0.228, 0.224];
d2 = max(cData2_ave(:));
cData2_ave = cData2_ave./d2;
cData2_std = [0.04, 0.08; 0.08, 0.032; 0.052, 0.152; 0.048, 0.048; 0.036, 0.032];
cData2_std = cData2_std./d2;
cData2_time = [1, 4, 7, 14, 28]./7;

%% Simulations to demonstrate how peak height affects simulations

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
% colororder(bone);
ylim([0 0.3])
legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65');
subplot(2,1,2)
yyaxis left
for i = 1:numC
    plot(tWeek,c3mrna(:,i),'-','Color',cmap(i,:));hold on
end
ylim([0 0.3])
legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65');


% Collagen content data to match (Fomovsky AJP 2010)
time_rat = [7;14;21;42];
AF_rat = [9.0; 16.5; 21.0; 26.5];
AFsd_rat = [4.0; 8.0; 9.0; 7.0];

tweek = days./7;

figure; 
title('Alter peak height');
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
plot(tweek,Carea,'k');
 for i = 1:numC
    plot(tweek,Carea(i,:),'-','Color',cmap(i,:));hold on
 end
legend({'0.35','0.4','0.45','0.5','0.55','0.6','0.65'},'Location','southeast');


%% Simulations to demostrate how baseline levels vary

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


% Collagen content data to match (Fomovsky AJP 2010)
time_rat = [7;14;21;42];
AF_rat = [9.0; 16.5; 21.0; 26.5];
AFsd_rat = [4.0; 8.0; 9.0; 7.0];

tweek = days./7;

figure; 
title('Alter baseline')
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
plot(tweek,Carea2,'k');
 for i = 1:numC
    plot(tweek,Carea2(i,:),'-','Color',cmap(i,:));hold on
 end
 xlim([-1 12])
legend({'0','0.05','0.1','0.15','0.2','0.25','0.2','0.25','0.3','0.35'},'Location','southeast');


%% Simulations that draw from a random distribution
sRand = [0.35,0.65];

for i = 1:100
    disp(strcat('Simulation Number',num2str(i)))
    [InputCsim, tInSim] = InputCurve_12_19rand(sRand(1),sRand(2));
    %running standard simulations
    % extract the parameters
    paramName = 'fib617_params';
    eval(strcat('[params,y0] = ',paramName,';'));
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
    
    % standard simulation
    outputIns = [101, 102, 103, 104, 95, 94, 96];
    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
    yI = real(interp1(t1,y1,tInSim));
    [c1,days] = MISimODEplotNP(sum(yI(:,[101,102]),2),tInSim);

    Carea(i,:) = c1;
    c1mrna(:,i) = yI(168:end,101);
    c3mrna(:,i) = yI(168:end,102);
    tWeek = (tInSim(168:end) - 168)./168;
    
    tgfb_Input(i,:) = InputCsim(1,:);
end



cmap = [0.250000000000000,0.250000000000000,0.375000000000000];


figure
subplot(2,1,1)
yyaxis left
plot(tWeek,c1mrna,'-','Color',cmap);hold on
ylim([0 0.11])
yyaxis right
errorbar(cData1_time,cData1_ave(:,1),cData1_std(:,1),'bo');hold on; 
errorbar(cData2_time,cData2_ave(:,1),cData2_std(:,1),'b*');hold on; 
ylim([0 1.2]);
xlim([0 9]);
%legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65','PMID: 11444923','PMID: 11557576');
subplot(2,1,2)
yyaxis left
plot(tWeek,c3mrna,'-','Color',cmap);hold on
ylim([0 0.11])
yyaxis right
hold on; errorbar(cData1_time,cData1_ave(:,2),cData1_std(:,2),'bo');
hold on; errorbar(cData2_time,cData2_ave(:,2),cData2_std(:,2),'b*');
ylim([0 1.2]);
xlim([0 9]);
%legend('0.35','0.4','0.45','0.5','0.55','0.6','0.65','PMID: 11444923','PMID: 11557576');
title(strcat('Randomly Sampled Peak Height with uniform range',num2str(sRand)));


% Collagen content data to match (Fomovsky AJP 2010)
time_rat = [7;14;21;42];
AF_rat = [9.0; 16.5; 21.0; 26.5];
AFsd_rat = [4.0; 8.0; 9.0; 7.0];

tweek = days./7;

figure; 
errorbar(time_rat./7,AF_rat,AFsd_rat,'ko');
axis([-1 9 0 40]); xlabel('Time (weeks)'); ylabel('Area Fraction (%)');
 hold on; plot(tweek,Carea,'k');
    plot(tweek,Carea,'-','Color',cmap);hold on
%title('nonlinear collagen model compared to data');
%legend({'Fomovsky 2010','0.35','0.4','0.45','0.5','0.55','0.6','0.65'},'Location','southeast');
title(strcat('Randomly Sampled Peak Height with uniform range',num2str(sRand)));



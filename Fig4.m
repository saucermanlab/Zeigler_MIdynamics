%% FIG4 driver effects under sttaic stimuli that most align with days 0,1,7,42
% RUN THE SCRIPT FOR FIGURE 5 FIRST TO GENERATE A DEPENDENT VARIABLE FOR
% THIS SCRIPT
%'ignoreDex' variable generated in script for Fig5
% updated ACZ 1.20.2020

clear 
clc
%% set up parameters
perturb = 10;
pHeight = 0.6;

% define the static inputs that mimic each day (from fig3 code)
% order of Inputs:
% 1:TGFB, 2:IL6, 3:IL1, 4:TNF, 5:NE, 6:ET1, 7:NP, 8:AngII, 9:PDGF
inputNames = {'TGFB','IL6','IL1','TNF','NE','ET1','NP','AngII','PDGF'};
numStatic = 2; % number of static inputs used to mimic dynamic
day1ins = [3,7];
day7ins = [1,6];
day42ins = [5,8];

[InputCsim,tInSim,inputNode,resNorm,resNormConvert]= InputCurve_12_19NP(pHeight,pHeight);


InputCsimcontrol = ones(9,length(tInSim))*0.1;

InputCsimDay1=InputCsimcontrol;
InputCsimDay7=InputCsimcontrol;
InputCsimDay42=InputCsimcontrol;


highInput=ones(numStatic,length(tInSim))*pHeight;

%day1
InputCsimDay1(day1ins,:)=highInput;
%day7
InputCsimDay7(day7ins,:)=highInput;
%day42
InputCsimDay42(day42ins,:)=highInput;





%% control (no upregulation) static conditions


%control day0 simulation
options=[];
modelName = 'dynamicODE';
[params,y0] = fib617_params(pHeight);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};

[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));
col_control=(yI(end,101)+yI(end,102)); %get control collagen level

figure;
imagesc(yI)

% control day 0.5
params = {rpar,tau,ymax,speciesNames,KI,InputCsimDay1,inputNode,tInSim};


[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));
col_control_day1=(yI(end,101)+yI(end,102)); %get control day 1 collagen level

figure;
imagesc(yI)

% control day 7
params = {rpar,tau,ymax,speciesNames,KI,InputCsimDay7,inputNode,tInSim};


[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));
col_control_day7=(yI(end,101)+yI(end,102)); 

figure;
imagesc(yI)


% control day 42
params = {rpar,tau,ymax,speciesNames,KI,InputCsimDay42,inputNode,tInSim};


[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));
col_control_day42=(yI(end,101)+yI(end,102)); 

figure;
imagesc(yI)

%% static screens
numSpec = length(speciesNames);

% day 0 screen
col_score_control = []; %initiate col_score_control

for i = 1:numSpec
    disp(['day0 # ',num2str(i),' of ', num2str(length(ymax))]) 
    ymax_new = ymax;
    ymax_new(i)=perturb;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsimcontrol,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    yI2 = real(interp1(t2,y2,tInSim));
    col_score_control(i) = (yI2(end,101)+yI2(end,102));%get CImRNA values overtime for each perturbation
end

col_score_control_fold = col_score_control-col_control; %get delta of perturbations on CImRNA


%day 1 screen
col_score_day1 = []; %initiate col_score_control

for i = 1:numSpec
    disp(['day1 # ',num2str(i),' of ', num2str(length(ymax))]) 
    ymax_new = ymax;
    ymax_new(i)=perturb;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsimDay1,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    yI2 = real(interp1(t2,y2,tInSim));
    col_score_day1(i) = (yI2(end,101)+yI2(end,102));%get CImRNA values overtime for each perturbation
end

col_score_day1_fold = col_score_day1-col_control_day1; %get delta of perturbations on CImRNA



%day 7 screen
col_score_day7 = []; %initiate col_score_control

for i = 1:numSpec
    disp(['day7 # ',num2str(i),' of ', num2str(length(ymax))]) 
    ymax_new = ymax;
    ymax_new(i)=perturb;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsimDay7,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    yI2 = real(interp1(t2,y2,tInSim));
    col_score_day7(i) = (yI2(end,101)+yI2(end,102));%get CImRNA values overtime for each perturbation
end

col_score_day7_fold = col_score_day7-col_control_day7; %get delta of perturbations on CImRNA

%day 42 screen
col_score_day42 = []; %initiate col_score_control

for i = 1:numSpec
    disp(['day42 # ',num2str(i),' of ', num2str(length(ymax))]) 
    ymax_new = ymax;
    ymax_new(i)=perturb;
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsimDay42,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    yI2 = real(interp1(t2,y2,tInSim));
    col_score_day42(i) = (yI2(end,101)+yI2(end,102));%get CImRNA values overtime for each perturbation
end

col_score_day42_fold = col_score_day42-col_control_day42; %get delta of perturbations on CImRNA


%% plots with ymax perturbations

%%remove the input nodex
%input_index=[2,10,17,20,27,34,37,42,50,55,58];
index=[1:107];
load('ignoreDex')
%index([input_index])=[];%remove input nodes
index([ignoreDex])=[];%remove input nodes
col_score_control_fold=col_score_control_fold(index);
col_score_day1_fold=col_score_day1_fold(index);
col_score_day7_fold=col_score_day7_fold(index);
col_score_day42_fold=col_score_day42_fold(index);

speciesNames=speciesNames(index);

load('ignoreDex') % index of nodes to ignore based on heatmap values in fig4


%% Plot screens

numSpec = length(speciesNames);
l = ones(1,numSpec);
[sort_col_score_control_fold,I] = sort(col_score_control_fold);

%Day 0
figure
bar(real(col_score_control_fold(I)))
hold on
plot(0,'--r');
title('Drivers of collagen I+III in Control context');
ylabel('Delta Collagen mRNA sum(CImRNA&CIIImRNA)');
xlabel('overexpressed node');
set(gca,'XTick',1:numSpec);
set(gca,'XTickLabel',speciesNames(I));
ylim([-0.5 2]);
xtickangle(270);
l = ones(1,numSpec);

%Day 1
figure
bar(real(col_score_day1_fold(I)))
hold on
plot(0,'--r');
title(['Drivers of collagen I+III', inputNames(day1ins)]);
ylabel('Delta Collagen mRNA sum(CImRNA&CIIImRNA)');
xlabel('overexpressed node');
set(gca,'XTickLabel',[]);
set(gca,'XTick',1:numSpec);
set(gca,'XTickLabel',speciesNames(I));
ylim([-0.5 2]);
xtickangle(270);

%Day 7
figure
bar(real(col_score_day7_fold(I)))
hold on
plot(0,'--r');
title(['Drivers of collagen I+III', inputNames(day7ins)]);
ylabel('Delta Collagen mRNA sum(CImRNA&CIIImRNA)');
xlabel('overexpressed node');
set(gca,'XTick',1:numSpec);
set(gca,'XTickLabel',speciesNames(I));
ylim([-0.5 2]);
xtickangle(270);

%Day 42
figure
bar(real(col_score_day42_fold(I)))
hold on
plot(0,'--r');
title(['Drivers of collagen I+III', inputNames(day42ins)]);
ylabel('Delta Collagen mRNA sum(CImRNA&CIIImRNA)');
xlabel('overexpressed node');
set(gca,'XTick',1:numSpec);
set(gca,'XTickLabel',speciesNames(I));
ylim([-0.5 2]);
xtickangle(270);

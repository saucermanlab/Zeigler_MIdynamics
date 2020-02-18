%% Generate figure 3: PCA analysis and comparison to static stimuli
% PCA plots are generated in R from exported files
% also generates figure S5: network activation in static conditions
% last updated ACZ 1.19.2020

clear
clc

dayTick = [168,180,216,240,264,288,312,336,504,672,840,1008,1176]; %defining the 0day, 1 day, 7day, and 42 day time points, day 1 tick is 
outputIns = [92:98, 84, 86, 87, 91, 100:107]; % all phenotypic outputs

%deltaIn = peak input weight or static input weight
deltaIn = 0.6;


%% standard dynamic simulation
[InputCsim,tInSim,inputNode] = InputCurve_12_19NP(deltaIn,deltaIn);

% % extract the parameters
[params,y0] = fib617_params(deltaIn);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
% 
% FOR 0.6 MODEL
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
yI = real(interp1(t1,y1,tInSim));

% % FOR ENSEMBLE MODEL
% load aveNet_st05
% yI = aveNet;

yP = yI(dayTick(1),outputIns);
yE = yI(dayTick(2),outputIns);
yM = yI(dayTick(8),outputIns);
yL = yI(dayTick(end),outputIns);

dayLabels={'Day0','Day0.5','Day7','Day42'};


%% Set up labels and pairwise combinations
inputNodeW = {2,4,5,6,7,9,10,1,8};
inputNames={'TGFB','IL6','IL1','TNFa','NE','ET1','NP','AngII','PDGF'};
inputNode=[2 4 5 6 7 9 10 1 8];

C2=nchoosek(1:9,2);
comboNames = inputNames;


%get label names and values for stimulus pairs
for i = 1:36
    index=C2(i,:);
    inputNodeW{end+1}=inputNode(index);
    comboNames{end+1}=[inputNames{index(1)},'+',inputNames{index(2)}];
end
comboNames{end+1} = 'Control';


%% run static simulations 
numIn = length(inputNodeW);
yOut=[];
for i = 1:numIn
    disp(i)
    tInSim = [0:1:500];
    InputCsim = ones(1,length(tInSim)).*deltaIn;
    Node = inputNodeW{i};
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,Node,tInSim};
    options = [];
    [t2,y2] = ode15s(@dynamicODE,[0 500],y0,options,params);
    y2 = real(y2);
    yOut(:,i) = y2(end,outputIns);
    yAll(:,i) = y2(end,:);
end 
% No stimulation
    tInSim = [0:1:500];
    InputCsim = ones(length(inputNode),length(tInSim)).*0.1;
    Node = [inputNode];    
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,Node,tInSim};
    options = [];
    [t3,y3] = ode15s(@dynamicODE,[0 500],y0,options,params);
    y3 = real(y3);
    yOut(:,numIn+1) = y3(end,outputIns);
    yAll(:,numIn+1) = y3(end,:);
    %combine all outputs
    allOuts = [yP;yE;yM;yL;yOut'];
    
allOuts2 = [yI(dayTick,outputIns);yOut']; %only used for PCA because has more time points


 %% euclidean distance from each day 

D=squareform(pdist(real(allOuts)));

day0dist=D(5:end,1);
day1dist=D(5:end,2);
day7dist=D(5:end,3);
day42dist=D(5:end,4);

%% display the lowest eucD values for each dynamic point
doubleIndex = [find(day0dist==min(day0dist)), find(day1dist==min(day1dist)), ...
    find(day7dist==min(day7dist)), find(day42dist==min(day42dist))];

doubleMatch = comboNames(doubleIndex);

sL = length(inputNames);
singleIndex = [find(day0dist==min(day0dist([1:sL,end]))), find(day1dist==min(day1dist([1:sL,end]))),...
    find(day7dist==min(day7dist([1:sL,end]))), find(day42dist==min(day42dist([1:sL,end])))];

singleMatch =comboNames(singleIndex);

%doubles
disp(strcat('day0 double match = ',doubleMatch{1}))
disp(strcat('day0.5 double match = ',doubleMatch{2}))
disp(strcat('day7 double match = ',doubleMatch{3}))
disp(strcat('day42 double match = ',doubleMatch{4}))
%singles
disp(strcat('day0 single match = ',singleMatch{1}))
disp(strcat('day0.5 single match = ',singleMatch{2}))
disp(strcat('day7 single match = ',singleMatch{3}))
disp(strcat('day42 single match = ',singleMatch{4}))

%% barplots of euclidean distance vs all 4 timepoints
fig=figure;
bar(real(day0dist))
set(gca,'XTick',1:length(comboNames));
set(gca,'XTickLabel',comboNames,'fontsize',10);
xtickangle(270)
title('Euclidean Distance from day 0')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 16 4];

fig=figure;
bar(real(day1dist))
set(gca,'XTick',1:length(comboNames));
set(gca,'XTickLabel',comboNames,'fontsize',10);
xtickangle(270)
title('Euclidean Distance from day 0.5')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 16 4];

fig=figure;
bar(real(day7dist))
set(gca,'XTick',1:length(comboNames));
set(gca,'XTickLabel',comboNames,'fontsize',10);
xtickangle(270)
title('Euclidean Distance from day 7')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 16 4];

fig=figure;
bar(real(day42dist))
set(gca,'XTick',1:length(comboNames));
set(gca,'XTickLabel',comboNames,'fontsize',10);
xtickangle(270)
title('Euclidean Distance from day 42')
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 16 4];
    
   



%% reduced heatmap, only nodes of fibrotic interest
fig=figure;
subplot(1,3,1);
imagesc(real(allOuts(1:4,:)'),[0 1]);

colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(outputIns));
set(gca,'YTickLabel',speciesNames(outputIns),'fontsize',13);
xlabel('Timepoint or Stimulus');
set(gca,'XTick',1:length(dayLabels));
set(gca,'XTickLabel',dayLabels,'fontsize',10);
ylabel('Output');
title(['Input Simulations and Timepoints']);
xtickangle(70)
colorbar('Location','eastoutside');

subplot(1,3,2);
imagesc(real(allOuts(singleIndex+4,:)'),[0 1]); % add 4 to index to account for dynamic time points

colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(outputIns));
set(gca,'YTickLabel',speciesNames(outputIns),'fontsize',13);
xlabel('Timepoint or Stimulus');
set(gca,'XTick',1:length(singleMatch));
set(gca,'XTickLabel',singleMatch,'fontsize',10);
ylabel('Output');
title(['Input Simulations and Timepoints']);
xtickangle(70)
colorbar('Location','eastoutside');

subplot(1,3,3);
imagesc(real(allOuts(doubleIndex+4,:)'),[0 1]); % add 4 to index to account for dynamic time points

colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(outputIns));
set(gca,'YTickLabel',speciesNames(outputIns),'fontsize',13);
xlabel('Timepoint or Stimulus');
set(gca,'XTick',1:length(doubleMatch));
set(gca,'XTickLabel', doubleMatch,'fontsize',10);
ylabel('Output');
title(['Input Simulations and Timepoints']);
xtickangle(70)
colorbar('Location','eastoutside');


%% establish labels for all timepoints/stimuli

dayLabels2={'Day0','Day0.5','Day2','Day3','Day4','Day5','Day6','Day7','Day14','Day21','Day28','Day35','Day42'};

labels2 = [dayLabels2,comboNames];

% create variable that defines different color groups
colorCode = [0,2,1];
highlight = repelem(colorCode,[length(dayLabels2), length(inputNames), (length(comboNames) - length(inputNames))]);
highlight(doubleIndex+length(dayLabels2)) = 3; % paired input matches
highlight(singleIndex+length(dayLabels2)) = 4; % single input matches
highlight(end) = 5; %control

nodeMatch = zeros(length(labels2),1);
nodeMatch([1:length(dayLabels2), doubleIndex+length(dayLabels2), singleIndex+length(dayLabels2)]) = 1;




%% create and write csv for R plotting
[coeff,score,latent,~,explained] = pca(real(allOuts2),'Centered','on','NumComponents',4); %Run PCA on 4 components
%create Cell array to export to R
exportScore=table(score(:,1),score(:,2),score(:,3),labels2',highlight',nodeMatch);
exportLoad = table(coeff(:,1),coeff(:,2),coeff(:,3),speciesNames(outputIns)');
% export table as CSV 
write(exportScore,'PCAscores.csv') %this is fininshed in R
writetable(exportLoad,'PCAload.csv')

% 
%     %from here, open the R script 'Fig 3A' and run the script to generate
%     %the plots used in the paper

%% Generate figure S5 showing network activity for all static conditions

AllNet = [yI(dayTick,:)',yAll];
figure
imagesc(real(AllNet),[0 1]);
colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',speciesNames,'fontsize',13);
xlabel('Stimulus');
set(gca,'XTick',1:length(labels2));
set(gca,'XTickLabel',labels2,'fontsize',7);
ylabel('Phenotypic Output');
title(['Input Simulations and Timepoints']);
xtickangle(270);
colorbar('Location','eastoutside');


%% Generate figure showing similarity of close time points

dayTick2 = [50, 504] +168;

closeNet = yI(dayTick2,outputIns);
figure
imagesc(real(closeNet'),[0 1]);
colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',speciesNames(outputIns),'fontsize',13);
xlabel('Stimulus');
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'day 2.1','day21'},'fontsize',7);
ylabel('Phenotypic Output');
title(['Input Simulations and Timepoints']);
xtickangle(270);
colorbar('Location','eastoutside');


%% supplemental figure 3
%get curves
[InputCsim,tInSim,inputNode] = InputCurve_1_19(0.6,0.6);

% extract the parameters
modelName = 'dynamicODE';
paramName = 'fib617_params';
eval(strcat('[params,y0] = ',paramName,';'));
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};

%% standard dynamic simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
yI = real(interp1(t1,y1,tInSim));
disp('done')

inputNodeW = {2,4,5,6,7,9,10,1,8};
inputNames={'TGFB','IL6','IL1','TNFa','NE','ET1','NP','AngII','PDGF'};
inputNode=[2 4 5 6 7 9 10 1 8];

%% get all combinations without redundant combos of stimulus
%pairs
C2=nchoosek(1:9,2);

%get label names and values for stimulus pairs
for i = 1:36
    index=C2(i,:);
    inputNodeW{end+1}=inputNode(index);
    inputNames{end+1}=[inputNames{index(1)},'+',inputNames{index(2)}];
end

dayTick = [150,180,336,1176]; %defining the 0day, 1 day, 7day, and 42 day time poin

yP = yI(dayTick(1),:);
yE = yI(dayTick(2),:);
yM = yI(dayTick(3),:);
yL = yI(dayTick(4),:);

labels ={};
lables2={};
dayLabels={'Day0','Day1','Day7','Day42'};
for i =1:4
labels{i}=dayLabels{i};
end
for i = 1:length(inputNames)
labels{i+4}=inputNames{i};
end
labels{4+length(inputNames)+1}='Control';

%% run static simulations iterively
deltaIn = 0.6;
numIn = length(inputNodeW);

for i = 1:numIn
    disp(i)
    tInSim = [0:1:500];
    InputCsim = ones(1,length(tInSim)).*deltaIn;
    Node = inputNodeW{i};
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,Node,tInSim};
    options = [];
    [t2,y2] = ode15s(@dynamicODE,[0 500],y0,options,params);
    y2 = real(y2);
    yOut(:,i) = y2(end,:);
end 
% No stimulation
    tInSim = [0:1:500];
    InputCsim = ones(length(inputNode),length(tInSim)).*0.1;
    Node = [inputNode];    
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,Node,tInSim};
    options = [];
    [t3,y3] = ode15s(@dynamicODE,[0 500],y0,options,params);
    y2 = real(y3);
    yOut(:,numIn+1) = y3(end,:);
    %combine all outputs
    allOuts = [yP;yE;yM;yL;yOut'];

    
    %% establish labels for all timepoints/stimuli, this part just indexes through names to apply name labels for the allOuts data    
labels ={};
dayLabels={'Day0','Day1','Day7','Day42'};
for i =1:4
labels{i}=dayLabels{i};
end
for i = 1:length(inputNames)
labels{i+4}=inputNames{i};
end
labels{4+length(inputNames)+1}='Control';


    
      %% make heatmap
     
    
fig=figure;

imagesc(real(allOuts)',[0 1]);
colormap(flipud(bone));
caxis([0 0.5]);
set(gca,'YTick',1:length(speciesNames));
set(gca,'YTickLabel',speciesNames,'fontsize',13);
xlabel('Stimulus');
set(gca,'XTick',1:length(labels));
set(gca,'XTickLabel',labels,'fontsize',7);
ylabel('Phenotypic Output');
title(['Input Simulations and Timepoints']);
xtickangle(270);
colorbar('Location','eastoutside');
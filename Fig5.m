%% figure 5: dynamic driver screen
% updated ACZ 1.20.2020


perturb = 10;
pHeight = 0.6;

%% standard simulation

[InputCsim,tInSim,inputNode] = InputCurve_12_19NP(pHeight,pHeight);


[params,y0] = fib617_params(pHeight);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
numSpec = length(speciesNames);


options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); 
yI = real(interp1(t1,y1,tInSim));

%% Iterative simulations with ymax upregulation
numSpec = length(speciesNames);
sens_up = zeros(2329,107); %initiate sens, keep line for increased speed in loop
pert_node_up=zeros(2329,107); %get expression overtime over each perturbed node
pert_node_up_norm=zeros(2329,107);
%for upregulgation, raise ymax of perturbed node to specified value
for i = 1:numSpec
    disp(['Upregulation # ',num2str(i),' of ', num2str(length(ymax))]) 
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    ymax_new = ymax;
    ymax_new(i)=perturb;
   
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    yI2 = real(interp1(t2,y2,tInSim));
     sens_up(:,i) = (yI2(:,101)+yI2(:,102)); %get CImRNA and CIIImRNA values overtime averaged for each perturbation
    pert_node_up(:,i)=(yI2(:,i)); %get level of perturbed node overtime
    pert_node_up_norm(:,i)=(yI2(:,i))-(yI(:,i)); %get level of perturbed node overtime minus the control
end
%% Normalize to delta vs normal simulation values
sens_up_norm = sens_up-(yI(:,101)+yI(:,102)); %subtract control CImRNA values from perturbed simulation to get difference at each timepoint

%% Get indicies of nodes that have a minimal effect below threshold, these will be removed from
%screening figure

ignoreDex=find(max(abs(sens_up_norm))<0.1); %nodes to remove due to minimal effect
save('ignoreDex','ignoreDex')
sens2=sens_up_norm;
sens2(:,ignoreDex)=[];
names2=speciesNames;
names2(ignoreDex)=[];

%% Custom Colormap (Blue => White => Red)
redd=[ones(1,100),1:-0.01:0];
bluu=[0:0.01:1,ones(1,100)];
whit=[0:.01:1,1,1:-.01:0];
cmap=[bluu',whit(1:201)',redd'];


%%
% plot the clustered sensitivity
tree = linkage(sens_up_norm(168:1008,:)'); %tree
dist = pdist(sens_up_norm(168:1008,:)'); %distances
leafOrder = optimalleaforder(tree,dist);
%background subtracted
figure  
colormap(cmap)
imagesc([sens_up_norm(1:1177,leafOrder)']) %1176=24(hr)*7(days)*7(weeks)
c=colorbar('southoutside', 'FontSize', 15);
c.Label.String = 'Delta Collagen I/III mRNA Sum';
%set(c,'position',[0.15 0.025 0.4 0.025])
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',[speciesNames(leafOrder)],'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Altered Node','FontSize',20);
caxis([-1,1]); 


% plot the clustered sensitivity
tree = linkage(sens2(200:1008,:)'); %tree
dist = pdist(sens2(200:1008,:)'); %distances
leafOrder = optimalleaforder(tree,dist);
%background subtracted
figure  
colormap(cmap)
imagesc([sens2(1:1177,leafOrder)']) %1176=24(hr)*7(days)*7(weeks)
c=colorbar('southoutside', 'FontSize', 15);
c.Label.String = 'Delta Collagen I/III mRNA Sum';
%set(c,'position',[0.15 0.025 0.4 0.025])
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',[names2(leafOrder)],'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Overexpressed Node','FontSize',20);
caxis([-1,1]); 


% figure S6 network activity in each overexpression

figure
colormap(flipud(bone))
imagesc(pert_node_up_norm(1:1177,:)')
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',names2(leafOrder),'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Altered Node','FontSize',20);


%% Ensemble modeling
peak = 0.6;
stdev = 0.05;

for i = 1:100
    [InputCsim,tInSim,inputNode,resNorm,resNormConvert,weightR(i,:)] = InputCurve_12_19rand(peak, stdev);
    params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};
    
    options = [];
    [t1,y1] = ode15s(@dynamicODE,[0 2329],y0,options,params);
    yI = real(interp1(t1,y1,tInSim));
    
    for j = 1:numSpec
        disp(['Ensemble',num2str(i),'Upregulation # ',num2str(j),' of 107']) 
        [rpar,tau,ymax,speciesNames,KI]=params{:};
        ymax_new = ymax;
        ymax_new(j)=perturb;

        paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
        [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
        yI2(:,:,i) = real(interp1(t2,y2,tInSim))';
        sens_ens(j,:,i) = (yI2(:,101)+yI2(:,102)) - (yI(:,101)+yI(:,102)); %get CImRNA and CIIImRNA values overtime averaged for each perturbation
     end
  
end

sens_mean = mean(sens_ens(:,1:end-1,:),3);


tree = linkage(sens_mean); %tree
dist = pdist(sens_mean); %distances
leafOrder = optimalleaforder(tree,dist);
%background subtracted
figure;  
colormap(cmap)
imagesc(sens_mean(leafOrder,1:1177)) %1176=24(hr)*7(days)*7(weeks)
c=colorbar('southoutside', 'FontSize', 15);
c.Label.String = 'Delta Collagen I/III mRNA Sum';
%set(c,'position',[0.15 0.025 0.4 0.025])
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',speciesNames(leafOrder),'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Altered Node','FontSize',20);
caxis([-1,1]); 



sens_meanDim = sens_mean;
sens_meanDim(ignoreDex,:) = [];
tree = linkage(sens_meanDim); %tree
dist = pdist(sens_meanDim); %distances
leafOrder = optimalleaforder(tree,dist);
%background subtracted
figure;  
colormap(cmap)
imagesc(sens_meanDim(leafOrder,1:1177)) %1176=24(hr)*7(days)*7(weeks)
c=colorbar('southoutside', 'FontSize', 15);
c.Label.String = 'Delta Collagen I/III mRNA Sum';
%set(c,'position',[0.15 0.025 0.4 0.025])
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',names2(leafOrder),'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Altered Node','FontSize',20);
caxis([-1,1]); 


% figure S6 network activity in each overexpression
yImean = mean(yI2,3);

figure
colormap(flipud(bone))
imagesc(yImean(:,1:1177))
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',names2(leafOrder),'FontSize',7);
set(gca,'XTick',0:168:1176);
set(gca,'XTickLabel',-1:1:6);
xlabel('Time (Weeks)','FontSize',20);
ylabel('Altered Node','FontSize',20);




% Sensitivity analysis

perturb = 0.5;

% extract the parameters

paramName = 'fib617_params';
eval(strcat('[params,y0] = ',paramName,';'));
[rpar,tau,ymax,speciesNames,KI]=params{:};

w = rpar(1,:);
n = rpar(2,:);
EC50 = rpar(3,:);


%% Control


tInSim = [0:1:500];
InputCsim = ones(1,length(tInSim)).*0.1;
inputNode = [10];    
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};  

% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yEnd0 = sum(y1(end,[101,102]));






%full sensitivity analysis of full KIn
for i = 1:length(ymax)
    disp(['Upregulation # ',num2str(i),' of ', num2str(length(ymax))]) 
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    ymax_new = ymax;
    ymax_new(i)=10;
   
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 

    colCon(i,:) = sum(y2(end,[101,102]))/yEnd0;
end


%% IL6+AngII
%inputs 1 and 4
%perturbations above the control conditions
tInSim = [0:1:500];
InputCsim = ones(1,length(tInSim)).*perturb;
inputNode = [1,4];    
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};  

% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yEnd0 = sum(y1(end,[101,102]));



%full sensitivity analysis of full KIn
for i = 1:length(ymax)
    disp(['Upregulation # ',num2str(i),' of ', num2str(length(ymax))]) 
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    ymax_new = ymax;
    ymax_new(i)=10;
   
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    col05(i,:) = sum(y2(end,[101,102]))/yEnd0;
end


%% TGFB+IL1
% 2 and 5 

%perturbations above the control conditions
tInSim = [0:1:500];
InputCsim = ones(1,length(tInSim)).*perturb;
inputNode = [2,5];    
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim}; 

% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yEnd0 = sum(y1(end,[101,102]));



%full sensitivity analysis of full KIn
for i = 1:length(ymax)
    disp(['Upregulation # ',num2str(i),' of ', num2str(length(ymax))]) 
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    ymax_new = ymax;
    ymax_new(i)=10;
   
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    col7(i,:) = sum(y2(end,[101,102]))/yEnd0;
end


%% AngII + PDGF
% 1 and 8

%perturbations above the control conditions
tInSim = [0:1:500];
InputCsim = ones(1,length(tInSim)).*perturb;
inputNode = [1,8];    
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim}; 

% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yEnd0 = sum(y1(end,[101,102]));



%full sensitivity analysis of full KIn
for i = 1:length(ymax)
    disp(['Upregulation # ',num2str(i),' of ', num2str(length(ymax))]) 
    [rpar,tau,ymax,speciesNames,KI]=params{:};
    ymax_new = ymax;
    ymax_new(i)=10;
   
    paramsNew = {rpar,tau,ymax_new,speciesNames,KI,InputCsim,inputNode,tInSim};
    [t2,y2] = ode15s(@dynamicODE,[0 2328],y0,options,paramsNew); 
    col42(i,:) = sum(y2(end,[101,102]))/yEnd0;
end


%% plots


% plot the single stimulus
numSpec = length(speciesNames);
l = ones(1,numSpec);
[scolCon,I] = sort(colCon);
figure
subplot(4,1,1)
bar(real(colCon(I)))
hold on
plot(l,'--r');
title('drivers of collagen I in Control context');
ylabel('fold change of Collagen mRNA');
xlabel('Overexpressed node');
set(gca,'XTick',1:numSpec);
ylim([0 2]);  

l = ones(1,numSpec);
subplot(4,1,2)
bar(real(col05(I)))
hold on
plot(l,'--r');
title('drivers of collagen I in AngII+IL6 context');
ylabel('fold change of fibrotic outputs');
xlabel('knocked down node');
set(gca,'XTick',1:numSpec);
ylim([0 2]);

subplot(4,1,3)
bar(real(col7(I)))
hold on
plot(l,'--r');
title('drivers of collagen I in TGFB+IL1 context');
ylabel('fold change of fibrotic outputs');
xlabel('knocked down node');
set(gca,'XTick',1:numSpec);
ylim([0 2]);

subplot(4,1,4)
bar(real(col42(I)))
hold on
plot(l,'--r');
title('drivers of collagen I in AngII+PDGF context');
ylabel('fold change of fibrotic outputs');
xlabel('knocked down node');
set(gca,'XTick',1:numSpec);
ylim([0 2]);
set(gca,'XTickLabel',speciesNames(I));
xticklabel_rotate;

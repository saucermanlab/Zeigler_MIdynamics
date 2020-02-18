% Supplemental figure 4 for dynamic MI
% updated ACZ 1.20.2020

peak = 0.5;

%generate input curves
[InputCsim,tInSim,inputNode] = InputCurve_12_19NP(peak,peak);


[params,y0] = fib617_params(peak);
[rpar,tau,ymax,speciesNames,KI]=params{:};
params = {rpar,tau,ymax,speciesNames,KI,InputCsim,inputNode,tInSim};

numSpec = length(speciesNames);

% standard simulation
options = [];
[t1,y1] = ode15s(@dynamicODE,[0 2328],y0,options,params); %y0 loaded previously
yI = real(interp1(t1,y1,tInSim));
%confirmatory plot (not saved) to show that inputs are changing overtime
figure;
imagesc(real(yI)')
caxis([0 1])
xlim([0,1682])
ylabel('Network Node')
xlabel('Time Post-MI (Weeks)')
colorbar
set(gca,'YTick',1:numSpec,'FontSize',10);
set(gca,'YTickLabel',[speciesNames],'FontSize',7);
set(gca,'XTick',0:168:1512);
set(gca,'XTickLabel',-1:1:8);
disp('done')
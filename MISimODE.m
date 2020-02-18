function [c1,days,time_rat,AF_rat,AFsd_rat] = MISimODE(yI,tInsim,peakCol)
% Use the model simulations to inform a simple model of collagen accumulation

% Will sent fibroblast density data in WJR_InfarctTurnoverDynamics.xlsx,
%   I approximated curve as linear rise in 1st 7 days, constant afterwards
%   and normalized so that max fibroblast concentration is 1.0
days =  linspace(-7,89,900)';
fibro = [ones(1,71)*0.1,linspace(0.1,1,70),ones(1,759)]';
% 
% days =  linspace(-7,97,4001)';
% fibro = [ones(1,71)*0.1,linspace(0.1,1,70),ones(1,3860)]'; 

%  note: these values are copied from the optimization performed in
%  MISimODEplot.m, see notes there for more info
kgen = 1.92;
kdeg = 0.024;


rawCol1 = yI';
tday = (tInsim -168)./24;
rawCol2 = interp1(tday,rawCol1,days);

% raw collagen is scaled to the peak of the standard simulation
kgt = (1-0.05).*(rawCol2 - min(rawCol2)) / (peakCol - min(rawCol2)) + 0.05; 

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


dcf = kgen*kgt.*fibro;
week = days./7;


end

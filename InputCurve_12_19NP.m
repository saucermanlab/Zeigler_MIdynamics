
function [InputCsim,tInSim,inputNode,resNorm,resNormConvert] = InputCurve_12_19NP(maxWeightLow, maxWeightHigh)
% InputCsim is a 9x2329 vector with the weights of the 9 input cytokines at all 2329 time points. 
% inputNode is a vector listing the indicies of the input reactions. 

tIn = [0:1:2160];
resNormConvert = [];

weightAxisDifferenceHigh = maxWeightHigh - 0.1;
weightAxisDifferenceLow = maxWeightLow - 0.1;
calc_FCmin = @(fcMax, maxWeight, weightAxisDifference) (-fcMax + 1)*(maxWeight/weightAxisDifference)+fcMax;

%% TGFB
y1_tgfb = [-10, 10, 174, 149, 72];
t1_tgfb = [1, 24, 96, 672, 1344];


y = [y1_tgfb];
t = [t1_tgfb];

[tSorted, tSortOrder] = sort(t);
ySorted = y(tSortOrder);

tSorted = tSorted(2:end);
ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fTGFB = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));
% fTGFB = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)))+...
%     p(4)*exp(-t/p(5)).*(1-exp(-t/p(6)));

paramsHandTuned_TGFB = [209.126999999910,1000,50];
% paramsHandTuned_TGFB = [440.673000000000, 170, 1500, 25, 2000, 5900];
%plot(tSorted,ySorted,'o',tIn,fTGFB(paramsHandTuned_TGFB,tIn))
%yTGFB = fTGFB(paramsHandTuned_TGFB,tIn);

FCmax_TGFB = max(y);
FCmin_TGFB = calc_FCmin(FCmax_TGFB, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_TGFB, 1, FCmax_TGFB*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 

% x = [-2.6298, 1, 20.1063];
% y = [0, 0.1, 0.63];

fTGFB_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fTGFB_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yTGFB = fTGFB_Convert(p, fTGFB(paramsHandTuned_TGFB,tIn)); 

while abs(max(yTGFB) - maxWeightHigh) > 0.00001
    %disp('TGFB');
    %disp(max(yTGFB));
    if (max(yTGFB) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_TGFB(1) - 0.001;
    else
        magTerm = paramsHandTuned_TGFB(1) + 0.001;
    end
    paramsHandTuned_TGFB(1) = magTerm;
    yTGFB = fTGFB_Convert(p, fTGFB(paramsHandTuned_TGFB,tIn));
end
    
% while  abs(yTGFB(t1_tgfb(4)) - 0.3) > 0.001
%     if yTGFB(t1_tgfb(4)) - 0.3 < 0
%         magTerm = paramsHandTuned_TGFB(4) + 0.1;
%     else
%         magTerm = paramsHandTuned_TGFB(4) - 0.1;
%     end
%     paramsHandTuned_TGFB(4) = magTerm;
%     yTGFB = fTGFB_Convert(p, fTGFB(paramsHandTuned_TGFB,tIn));
% end
% 
% while abs(max(yTGFB) - maxWeightHigh) > 0.00001
%     %disp('TGFB');
%     %disp(max(yTGFB));
%     if (max(yTGFB) - maxWeightHigh) > 0
%         magTerm = paramsHandTuned_TGFB(1) - 0.001;
%     else
%         magTerm = paramsHandTuned_TGFB(1) + 0.001;
%     end
%     paramsHandTuned_TGFB(1) = magTerm;
%     yTGFB = fTGFB_Convert(p, fTGFB(paramsHandTuned_TGFB,tIn));
% end

resNormTGFB = sum((fTGFB(paramsHandTuned_TGFB,tSorted) - ySorted).^2);
rmse_TGFB = sqrt(mean((fTGFB(paramsHandTuned_TGFB,tSorted) - ySorted).^2));

%% IL6
y1_il6 = [16.1, 38.4, 45.8, 10.8, 2.8, 2.5, 1.3, 1.3, 1.3, 2.5, 1.6];
t1_il6 = [3, 6, 12, 24, 72, 144, 216, 504, 672, 1344, 2016];

y = [y1_il6];
t = [t1_il6];

[tSorted, tSortOrder] = sort(t);
ySorted = y(tSortOrder);

tSorted = tSorted(2:end);
ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fIL6 = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));%+...
    %p(4)*exp(-t/p(5)).*(1-exp(-t/p(6)));

paramsHandTuned_IL6 = [6687.72799901963, 10, 550];
%plot(tSorted,ySorted,'o',tIn,fIL6(paramsHandTuned_IL6,tIn))
%yIL6 = fIL6(paramsHandTuned_IL6,tIn);

FCmax_IL6 = max(y);
FCmin_IL6 = calc_FCmin(FCmax_IL6, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_IL6, 1, FCmax_IL6*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 

% x = [-14, 1, 79.8];
% y = [0, 0.1, 0.63];

fIL6_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fIL6_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yIL6 = fIL6_Convert(p, fIL6(paramsHandTuned_IL6,tIn)); 

while abs(max(yIL6) - maxWeightHigh) > 0.00001
    %disp('IL6');
    %disp(max(yIL6));
    if (max(yIL6) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_IL6(1) - 0.001;
    else
        magTerm = paramsHandTuned_IL6(1) + 0.001;
    end
    paramsHandTuned_IL6(1) = magTerm;
    yIL6 = fIL6_Convert(p, fIL6(paramsHandTuned_IL6,tIn));
end

resNormIl6 = sum((fIL6(paramsHandTuned_IL6,tSorted) - ySorted).^2);
rmse_IL6 = sqrt(mean((fIL6(paramsHandTuned_IL6,tSorted) - ySorted).^2));

%% IL1 
y1_il1 = [5.7, 25.6, 32.7, 15.7, 2.6, 1.9, 2.6, 3.1, 2.1, 2.1];
t1_il1 = [3, 6, 12, 24, 72, 144, 216, 672, 1344, 2016];

y = [y1_il1];
t = [t1_il1];

[tSorted, tSortOrder] = sort(t);
ySorted = y(tSortOrder);

tSorted = tSorted(2:end);
ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fIL1 = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

paramsHandTuned_IL1 = [66.1659999999568, 18, 6];
%plot(tSorted,ySorted,'o',tIn,fIL1(paramsHandTuned_IL1,tIn))
%yIL1 = fIL1(paramsHandTuned_IL1,tIn);

FCmax_IL1 = max(y);
FCmin_IL1 = calc_FCmin(FCmax_IL1, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_IL1, 1, FCmax_IL1*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 

% x = [-2.88, 1, 21.42];
% y = [0, 0.1, 0.63];

fIL1_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fIL1_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yIL1 = fIL1_Convert(p, fIL1(paramsHandTuned_IL1,tIn)); 

while abs(max(yIL1) - maxWeightHigh) > 0.00001
    %disp('IL1');
    %disp(max(yIL1));
    if (max(yIL1) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_IL1(1) - 0.001;
    else
        magTerm = paramsHandTuned_IL1(1) + 0.001;
    end
    paramsHandTuned_IL1(1) = magTerm;
    yIL1 = fIL1_Convert(p, fIL1(paramsHandTuned_IL1,tIn));
end

resNormIL1 = sum((fIL1(paramsHandTuned_IL1,tSorted) - ySorted).^2);
rmse_IL1 = sqrt(mean((fIL1(paramsHandTuned_IL1,tSorted) - ySorted).^2));

%% TNF 
y1_tnfa = [0.16, 0.34, 0.31, 0.22, 0.19]./0.01; %normalized to the control not to sham
t1_tnfa = [24, 96, 264, 672, 960];

y = [y1_tnfa];
t = [t1_tnfa];

[tSorted, tSortOrder] = sort(t);
ySorted = y(tSortOrder);

tSorted = tSorted(2:end);
ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fTNFa = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

% +...
%     p(4)*exp(-t/p(5)).*(1-exp(-t/p(6)));

%paramsHandTuned_TNF = [158, 170, 900, 150, 10, 100];
paramsHandTuned_TNFa = [39.9369999999981,1000,50];
% params_handTuned = [10, 49e2, 25];
%plot(tSorted,ySorted,'o',tIn,fTNF(paramsHandTuned_TNF,tIn))
%yTNFa = fTNFa(paramsHandTuned_TNFa,tIn);

FCmax_TNFa = max(y);
FCmin_TNFa = calc_FCmin(FCmax_TNFa, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_TNFa, 1, FCmax_TNFa*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 

% x = [-1, 1, 11.55];
% y = [0, 0.1, 0.63];

fTNFa_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fTNFa_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yTNFa = fTNFa_Convert(p, fTNFa(paramsHandTuned_TNFa,tIn)); 

while abs(max(yTNFa) - maxWeightHigh) > 0.00001
    %disp('TNFa');
    %disp(max(yTNFa));
    if (max(yTNFa) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_TNFa(1) - 0.001;
    else
        magTerm = paramsHandTuned_TNFa(1) + 0.001;
    end
    paramsHandTuned_TNFa(1) = magTerm;
    yTNFa = fTNFa_Convert(p, fTNFa(paramsHandTuned_TNFa,tIn));
end

resNormTNFa = sum((fTNFa(paramsHandTuned_TNFa,tSorted) - ySorted).^2);
rmse_TNFa = sqrt(mean((fTNFa(paramsHandTuned_TNFa,tSorted) - ySorted).^2));

%% NE 
% y1_ne = [1.04, 1.00, 1.02, 1.61];
% t1_ne = [336, 672, 1008, 1344];

y1_ne = [1.24];
t1_ne = [1334];
y2_ne = [1.80];
t2_ne = [1176];

y = [y1_ne, y2_ne];
t = [t1_ne, t2_ne];

[tSorted, tSortOrder] = sort(t1_ne);
ySorted = y1_ne(tSortOrder);

fNE = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

% fNE = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)))+...
%     p(4)*exp(-t/p(5)).*(1-exp(-t/p(6)));

paramsHandTuned_NE = [0.516969999992107, 3200, 750];
% paramsHandTuned_NE = [1.81459999999999, 10, 1];
% params_handTuned = [10, 49e2, 25];
%plot(tSorted,ySorted,'o',tIn,fNE(paramsHandTuned_NE,tIn))
%yNE = fNE(paramsHandTuned_NE,tIn);

FCmax_NE = max(y1_ne);
FCmin_NE = calc_FCmin(FCmax_NE, maxWeightLow, weightAxisDifferenceLow);

x = [FCmin_NE, 1, FCmax_NE*1.05];
y = [0, 0.1, maxWeightLow*1.05]; 

% x = [0.7441, 1, 2.3934];
% y = [0, 0.1, 0.63];

fNE_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fNE_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yNE = fNE_Convert(p, fNE(paramsHandTuned_NE,tIn));

while abs(max(yNE) - maxWeightLow) > 0.00001
    %disp('NE');
    %disp(max(yNE));
    if (max(yNE) - maxWeightLow) > 0
        magTerm = paramsHandTuned_NE(1) - 0.00001;
    else
        magTerm = paramsHandTuned_NE(1) + 0.00001;
    end
    paramsHandTuned_NE(1) = magTerm;
    yNE = fNE_Convert(p, fNE(paramsHandTuned_NE,tIn));
end

resNormNE = sum((fNE(paramsHandTuned_NE,tSorted) - ySorted).^2);
rmse_NE = sqrt(mean((fNE(paramsHandTuned_NE,tSorted) - ySorted).^2));

%% ET1  
y1_et = [20, 36, 22, 5];
t1_et = [72, 120, 168, 1008];

[tSorted, tSortOrder] = sort(t1_et);
ySorted = y1_et(tSortOrder);

% tSorted = tSorted(2:end);
% ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fET1 = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));


paramsHandTuned_ET1 = [90.0012800177651, 120, 60]; %changed p3 to 60 which is half the hours in 5 days (peak signal)
% params_handTuned = [10, 49e2, 25];
%plot(tSorted,ySorted,'o',tIn,fET1(paramsHandTuned_ET1,tIn))
%yET1 = fET1(paramsHandTuned_ET1,tIn);

FCmax_ET1 = max(y1_et);
FCmin_ET1 = calc_FCmin(FCmax_ET1, maxWeightLow, weightAxisDifferenceLow);

x = [FCmin_ET1, 1, FCmax_ET1*1.05];
y = [0, 0.1, maxWeightLow*1.05]; 


% x = [0.7629, 1, 2.2950];
% y = [0, 0.1, 0.63];

fET1_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fET1_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yET1 = fET1_Convert(p, fET1(paramsHandTuned_ET1,tIn));

while abs(max(yET1) - maxWeightLow) > 0.00001
    %disp('ET1');
    %disp(max(yET1));
    if (max(yET1) - maxWeightLow) > 0
        magTerm = paramsHandTuned_ET1(1) - 0.00001;
    else
        magTerm = paramsHandTuned_ET1(1) + 0.00001;
    end
    paramsHandTuned_ET1(1) = magTerm;
    yET1 = fET1_Convert(p, fET1(paramsHandTuned_ET1,tIn));
end

resNormET1 = sum((fET1(paramsHandTuned_ET1,tSorted) - ySorted).^2);
rmse_ET1 = sqrt(mean((fET1(paramsHandTuned_ET1,tSorted) - ySorted).^2));

%% BNP
y1_bnp = [14.6, 36.2, 22, 4.6];
t1_bnp = [72, 120, 168, 1008];

[tSorted, tSortOrder] = sort(t1_bnp);
ySorted = y1_bnp(tSortOrder);

% tSorted = tSorted(2:end);
% ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fBNP = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

paramsHandTuned_BNP = [90.5148900179280, 120, 60];
% params_handTuned = [10, 49e2, 25];
%plot(tSorted,ySorted,'o',tIn,fBNP(paramsHandTuned_BNP,tIn))
%yBNP = fBNP(paramsHandTuned_BNP,tIn);

FCmax_BNP = max(y1_bnp);
FCmin_BNP = calc_FCmin(FCmax_BNP, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_BNP, 1, FCmax_BNP*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 


% x = [-0.993, 1, 11.5133];
% y = [0, 0.1, 0.63];

fBNP_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fBNP_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yBNP = fBNP_Convert(p, fBNP(paramsHandTuned_BNP,tIn));

while abs(max(yBNP) - maxWeightHigh) > 0.00001
    %disp('BNP');
    %disp(max(yBNP));
    if (max(yBNP) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_BNP(1) - 0.00001;
    else
        magTerm = paramsHandTuned_BNP(1) + 0.00001;
    end
    paramsHandTuned_BNP(1) = magTerm;
    yBNP = fBNP_Convert(p, fBNP(paramsHandTuned_BNP,tIn));
end

resNormBNP = sum((fBNP(paramsHandTuned_BNP,tSorted) - ySorted).^2);
rmse_BNP = sqrt(mean((fBNP(paramsHandTuned_BNP,tSorted) - ySorted).^2));

%% AngII
y1_ang = [1.31]./0.20;
t1_ang = [1334];
y2_ang = [21.8]./5.2;
t2_ang = [504];

y = [y1_ang, y2_ang];
t = [t1_ang, t2_ang];

[tSorted, tSortOrder] = sort(t);
ySorted = y(tSortOrder);

tSorted = tSorted(2:end);
ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fAng = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

paramsHandTuned_Ang = [10.0797999997794, 3200, 750];
% params_handTuned = [10, 49e2, 25];
%plot(tSorted,ySorted,'o',tIn,fAng(paramsHandTuned_Ang,tIn))
%yAng = fAng(paramsHandTuned_Ang,tIn);

FCmax_Ang = max(y);
FCmin_Ang = calc_FCmin(FCmax_Ang, maxWeightLow, weightAxisDifferenceLow);

x = [FCmin_Ang, 1, FCmax_Ang*1.05];
y = [0, 0.1, maxWeightLow*1.05]; 


% x = [0.7259, 1, 2.4889];
% y = [0,0.1,0.63];

fAng_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fAng_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yAng = fAng_Convert(p, fAng(paramsHandTuned_Ang,tIn)); 

while abs(max(yAng) - maxWeightLow) > 0.00001
    %disp('Ang');
    %disp(max(yAng));
    if (max(yAng) - maxWeightLow) > 0
        magTerm = paramsHandTuned_Ang(1) - 0.00001;
    else
        magTerm = paramsHandTuned_Ang(1) + 0.00001;
    end
    paramsHandTuned_Ang(1) = magTerm;
    yAng = fAng_Convert(p, fAng(paramsHandTuned_Ang,tIn));
end

resNormAng = sum((fAng(paramsHandTuned_Ang,tSorted) - ySorted).^2);
rmse_Ang = sqrt(mean((fAng(paramsHandTuned_Ang,tSorted) - ySorted).^2));

%% PDGF
y1_pdgf = [0.25, 2.76, 3.64, 4.52, 2.80, 2.11];
t1_pdgf = [24, 72, 168, 336, 672, 1008];

[tSorted, tSortOrder] = sort(t1_pdgf);
ySorted = y1_pdgf(tSortOrder);

% tSorted = tSorted(2:end);
% ySorted = ySorted(2:end); % Since there are two, t = 0 & y = 1 points. 

fPDGF = @(p,t) 1+p(1)*exp(-t/p(2)).*(1-exp(-t/p(3)));

paramsHandTuned_PDGF = [17.9955999999834,400,560];
%plot(tSorted,ySorted,'o',tIn,fPDGF(paramsHandTuned_PDGF,tIn))
%yPDGF = fPDGF(paramsHandTuned_PDGF,tIn);

FCmax_PDGF = max(y1_pdgf);
FCmin_PDGF = calc_FCmin(FCmax_PDGF, maxWeightHigh, weightAxisDifferenceHigh);

x = [FCmin_PDGF, 1, FCmax_PDGF*1.05];
y = [0, 0.1, maxWeightHigh*1.05]; 


% x = [0.5782, 1, 3.2643];
% y = [0, 0.1, 0.63];

fPDGF_Convert = @(p,xData) p(1).*xData + p(2);

p0 = rand(1,2);
[p, resnorm] = lsqcurvefit(fPDGF_Convert,p0,x,y);
resNormConvert = [resNormConvert;resnorm];
yPDGF = fPDGF_Convert(p, fPDGF(paramsHandTuned_PDGF,tIn));

while abs(max(yPDGF) - maxWeightHigh) > 0.00001
    %disp('PDGF');
    %disp(max(yPDGF));
    if (max(yPDGF) - maxWeightHigh) > 0
        magTerm = paramsHandTuned_PDGF(1) - 0.0001;
    else
        magTerm = paramsHandTuned_PDGF(1) + 0.0001;
    end
    paramsHandTuned_PDGF(1) = magTerm;
    yPDGF = fPDGF_Convert(p, fPDGF(paramsHandTuned_PDGF,tIn));
end

resNormPDGF = sum((fPDGF(paramsHandTuned_PDGF,tSorted) - ySorted).^2);
rmse_PDGF = sqrt(mean((fPDGF(paramsHandTuned_PDGF,tSorted) - ySorted).^2));

% combine all the inputs
InputC = [yTGFB; yIL6; yIL1; yTNFa; yNE; yET1; yBNP; yAng; yPDGF];
InputCsim = [ones(9,168).*0.1, InputC];  % add a week's worth of normal stimulation
tInSim = [0:1:2328];
inputNode = [2, 4, 5, 6, 7, 9, 10, 1, 8]; % based on w index of input reaction
resNorm = [resNormTGFB; resNormIl6; resNormIL1; resNormTNFa;...
    resNormNE; resNormET1; resNormBNP; resNormAng; resNormPDGF]; 
rmse = [rmse_TGFB; rmse_IL6; rmse_IL1; rmse_TNFa;...
    rmse_NE; rmse_ET1; rmse_BNP; rmse_Ang; rmse_PDGF];
    
end


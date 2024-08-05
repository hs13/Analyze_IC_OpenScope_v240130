load('/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/ICexpts_revision23/openscope_psthavgall.mat')
%% log mean vs log var plot slopes and intercepts across trial types
% intercept (log fano factor) is highest for blank
% not significantly different between IC and LC
whichblock = 'ICwcfg1_presentations';
Nsessions = numel(nwbsessions);
slopes = zeros(Nsessions, numel(ICtrialtypes));
yintercepts = zeros(Nsessions, numel(ICtrialtypes));
for ises = 1:Nsessions
tempneu = sesneuall==ises;
tempspkavg = Ronavgall.(whichblock)(tempneu,:)*0.4;
tempspkvar = (Ronstdall.(whichblock)(tempneu,:)*0.4).^2;
for itt = 1:numel(ICtrialtypes)
    tempx = log10(tempspkavg(:,itt));
tempy = log10(tempspkvar(:,itt));
valneu = tempspkavg(:,itt)>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
slopes(ises,itt) = a;
yintercepts(ises,itt) = b;
end
end

ICtt2p = [0 106 107 110 111 1105 1109];
tt2p = ICtrialtypes(ismember(ICtrialtypes, ICtt2p));
figure; 
plot(1:numel(tt2p), 10.^(yintercepts(:,ismember(ICtrialtypes, ICtt2p))), '.-' )
set(gca, 'XTick', 1:numel(tt2p), 'XTickLabel', tt2p)

ICtt2p = [106 107 110 111];
[p,tbl,stats]=friedman( (yintercepts(:,ismember(ICtrialtypes, ICtt2p))) );
disp(p)
multcompare(stats)

ICtt2p = [106 107 110 111];
[p,tbl,stats]=friedman( (slopes(:,ismember(ICtrialtypes, ICtt2p))) );
disp(p)
multcompare(stats)

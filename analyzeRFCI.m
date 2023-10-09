% THIS CODE TREATS THE ENTIRE SPINNING GRATING AS ONE TRIAL, CONSISTENT
% WITH THE ANALYSIS DONE WITH 2P DATASET

% % 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 
% % 10-100s: which RFcenter, 1s: which direction
% RFCI = struct();
% tloi = psthtli>0 & psthtli<=1000;
% tempR = squeeze(1000*mean(psth.RFCI(tloi,:,:), 1))';
% temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
% RFCI.(visblocknames{b}) = analyzeRFCI(tempR, temptrialorder);

% 'RFcentersrel', 'cRFsigneurons', 'cRFind', 'iRFsigneurons', 'iRFind', ...
% 'RFcentersrel9', 'RFcenterinds9', 'cRFind9', 'pRFclassic9', 'cRFsigexcl9', ...
% 'iRFind9', 'pRFinverse9', 'iRFsigexcl9')

function RFCI = analyzeRFCI(R, trialorder, sponFR)

Nneurons = size(R,1);
if numel(sponFR) ~= Nneurons
    error('check sponFR')
end

Nrfs = 9;
Rrfclassic = zeros(Nneurons, Nrfs);
Rrfinverse = zeros(Nneurons, Nrfs);
% 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
for rfi = 1:Nrfs
    crftrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    Rrfclassic(:,rfi) = mean(R(:, crftrials),2);
    
    irftrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    Rrfinverse(:,rfi) = mean(R(:, irftrials),2);
end

[~, RFindclassic] = max(Rrfclassic, [], 2);
[~, RFindinverse] = max(Rrfinverse, [], 2);

Pkw_rfclassic = zeros(Nneurons, 1);
Pkw_rfinverse = zeros(Nneurons, 1);
sigmc_rfclassic = false(Nneurons, 1);
sigmc_rfinverse = false(Nneurons, 1);
for ci = 1:Nneurons
    %     indRFCIsize = floor(mod(rectrialorder, 10000)/1000);
    %     indRFCIloc = floor(mod(rectrialorder, 1000)/10);
    crftrials =  floor(trialorder/10000) == 0;
    [p,tbl,stats] = kruskalwallis(R(ci,crftrials), trialorder(crftrials), 'off');
    Pkw_rfclassic(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindclassic(ci)|c(:,2)==RFindclassic(ci);
    sigmc_rfclassic(ci) = all(c(temprows,6)<0.05);
    
    irftrials =  floor(trialorder/10000) == 1;
    [p,tbl,stats] = kruskalwallis(R(ci,irftrials), trialorder(irftrials), 'off');
    Pkw_rfinverse(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindinverse(ci)|c(:,2)==RFindinverse(ci);
    sigmc_rfinverse(ci) = all(c(temprows,6)<0.05);
end


pRrfclassic = NaN(Nneurons, Nrfs);
for rfi = 1:Nrfs
    crftrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;    
    if nnz(crftrials)==0
        continue
    end
    for ci = 1:Nneurons
        pRrfclassic(ci,rfi) = signrank(R(ci,crftrials)-sponFR(ci));
    end
end
pRrfinverse = NaN(Nneurons, Nrfs);
for rfi = 1:Nrfs
    irftrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    if nnz(irftrials)==0
        continue
    end
    for ci = 1:Nneurons
        pRrfinverse(ci,rfi) = signrank(R(ci,irftrials)-sponFR(ci));
    end
end

pRFclassic = pRrfclassic(sub2ind(size(pRrfclassic), (1:Nneurons)', RFindclassic));
pRFinverse = pRrfinverse(sub2ind(size(pRrfinverse), (1:Nneurons)', RFindinverse));

% find cells where RF is the only significant location. bonferroni-holm corrected
tempsorted = sort(pRrfclassic, 2);
RFsigexclclassic = pRFclassic*Nrfs<0.05 & pRFclassic==tempsorted(:,1) & sum(tempsorted .* [Nrfs:-1:1]<0.05, 2)==1;

tempsorted = sort(pRrfinverse, 2);
RFsigexclinverse = pRFinverse*Nrfs<0.05 & pRFinverse==tempsorted(:,1) & sum(tempsorted .* [Nrfs:-1:1]<0.05, 2)==1;

RFexclsigclassic = pRFclassic<0.05 & sum(pRrfclassic<0.05, 2)==1;
RFexclsiginverse = pRFinverse<0.05 & sum(pRrfinverse<0.05, 2)==1;

RFCIfields = {'Rrfclassic', 'Rrfinverse', 'RFindclassic', 'RFindinverse', ...
    'Pkw_rfclassic', 'Pkw_rfinverse', 'sigmc_rfclassic', 'sigmc_rfinverse', ...
    'pRrfclassic', 'pRrfinverse', 'pRFclassic', 'pRFinverse', ...
    'RFsigexclclassic', 'RFsigexclinverse', 'RFexclsigclassic', 'RFexclsiginverse'};
RFCI = struct();
for f = 1:numel(RFCIfields)
    RFCI.(RFCIfields{f}) = eval(RFCIfields{f});
end

end
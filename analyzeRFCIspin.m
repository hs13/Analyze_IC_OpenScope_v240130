% THIS CODE EACH ORIENTATION OF THE SPINNING GRATING AS ONE TRIAL
% note, orientation order is not radomized, so we expect
% adaptation/prediction

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

function RFCIspin = analyzeRFCIspin(R, trialorder, sponFR)

Nneurons = size(R,1);
if numel(sponFR) ~= Nneurons
    error('check sponFR')
end

Nrfs = 9;
Noris = 4;
trialrfpos = floor(mod(trialorder, 1000) / 10);
% trialori = mod(trialorder, 10);

RFCIspintrialtypes = reshape( 10000*(0:1)+reshape(1000 + 10*(0:Nrfs-1)+(1:4)',[],1) ,[],1);
[v,c]=uniquecnt(trialorder);
if ~isequal(v, RFCIspintrialtypes)
    error('not all trialtypes are represented');
end
% if ~isequal(mod(reshape(v,4,[]),10), repmat((1:4)',1,2*Nrfs))
%     error('revise the code since it assumes that v and RFCIspintrialtypes are equal');
% end
% minrepperori = min(reshape(c,4,[]), [], 2);
% [~, iorimaxrep] = max(minrepperori,[],1);
minrepclassic = min(c(floor(v/10000)==0));
minrepinverse = min(c(floor(v/10000)==1));
Rclassic = NaN(Nneurons, max(c), Noris, Nrfs);
Rinverse = NaN(Nneurons, max(c), Noris, Nrfs);
for rfi = 1:Nrfs
    for iori = 1:Noris
        ctype = 0*10000 + 1000 + (rfi-1)*10 + iori;
        crftrials = trialorder==ctype;
        Rclassic(:, 1:nnz(crftrials), iori, rfi) = R(:, crftrials);
        
        itype = 1*10000 + 1000 + (rfi-1)*10 + iori;
        irftrials = trialorder==itype;
        Rinverse(:, 1:nnz(irftrials), iori, rfi) = R(:, irftrials);
    end
end

Ravgclassic = NaN(Nneurons, Nrfs, Noris);
Ravginverse = NaN(Nneurons, Nrfs, Noris);
% 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
for rfi = 1:Nrfs
    for iori = 1:Noris
        ctype = 0*10000 + 1000 + (rfi-1)*10 + iori;
        crftrials = trialorder==ctype;
        Ravgclassic(:,rfi,iori) = mean(R(:, crftrials),2);
        
        itype = 1*10000 + 1000 + (rfi-1)*10 + iori;
        irftrials = trialorder==itype;
        Ravginverse(:,rfi,iori) = mean(R(:, irftrials),2);
    end
end

if ~( isequaln(Ravgclassic, squeeze(nanmean(permute(Rclassic,[1 2 4 3]),2))) && ...
        isequaln(Ravginverse, squeeze(nanmean(permute(Rinverse,[1 2 4 3]),2))) )
    error('Rclassic and Ravgclassic disagree: check code')
end

Rrfclassic = squeeze(mean(Ravgclassic,3));
Rrfinverse = squeeze(mean(Ravginverse,3));
[~, RFindclassic] = max(Rrfclassic, [], 2);
[~, RFindinverse] = max(Rrfinverse, [], 2);

Pkw_rfclassic = zeros(Nneurons, 1);
Pkw_rfinverse = zeros(Nneurons, 1);
sigmc_rfclassic = false(Nneurons, 1);
sigmc_rfinverse = false(Nneurons, 1);
for ci = 1:Nneurons
    crftrials =  floor(trialorder/10000)==0;
    [p,tbl,stats] = kruskalwallis(R(ci,crftrials), trialrfpos(crftrials), 'off');
    Pkw_rfclassic(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindclassic(ci)|c(:,2)==RFindclassic(ci);
    sigmc_rfclassic(ci) = all(c(temprows,6)<0.05);
    
    irftrials =  floor(trialorder/10000)==1;
    [p,tbl,stats] = kruskalwallis(R(ci,irftrials), trialrfpos(irftrials), 'off');
    Pkw_rfinverse(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindinverse(ci)|c(:,2)==RFindinverse(ci);
    sigmc_rfinverse(ci) = all(c(temprows,6)<0.05);
end

% friedman is non-parametric version of two-way anova 
% Friedman's test is similar to classical balanced two-way ANOVA, 
% but it tests only for column effects after adjusting for possible row effects. 
% It does not test for row effects or interaction effects. 
% Friedman's test is appropriate when columns represent treatments that are under study, 
% and rows represent nuisance effects (blocks) that need to be taken into account 
% but are not of any interest.
% friedman on position+orientation 
if isequal( Rclassic(:,all(~isnan(Rclassic), [1 3 4]),:,:), Rclassic(:,1:minrepclassic,:,:) ) && ...
        isequal( Rinverse(:,all(~isnan(Rinverse), [1 3 4]),:,:), Rinverse(:,1:minrepinverse,:,:) )
else
    error('check code')
end

Pfried_rfclassic = NaN(Nneurons, 1);
Pfried_rfinverse = NaN(Nneurons, 1);
sigfmc_rfclassic = false(Nneurons, 1);
sigfmc_rfinverse = false(Nneurons, 1);
for ci = 1:Nneurons
    tempmat = squeeze(reshape(Rclassic(ci,1:minrepclassic,:,:),[],Nrfs));
    [p,tbl,stats] = friedman(tempmat, minrepclassic, 'off');
    Pfried_rfclassic(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindclassic(ci)|c(:,2)==RFindclassic(ci);
    sigfmc_rfclassic(ci) = all(c(temprows,6)<0.05);

    tempmat = squeeze(reshape(Rinverse(ci,1:minrepinverse,:,:),[],Nrfs));
    [p,tbl,stats] = friedman(tempmat, minrepinverse, 'off');
    Pfried_rfinverse(ci) = p;
    
    [c,m,h]=multcompare(stats,'display','off');
    temprows = c(:,1)==RFindinverse(ci)|c(:,2)==RFindinverse(ci);
    sigfmc_rfinverse(ci) = all(c(temprows,6)<0.05);
end

pRrfclassic = NaN(Nneurons, Nrfs);
for rfi = 1:Nrfs
    crftrials = floor(trialorder/10000)==0 & floor(mod(trialorder, 1000) / 10)==rfi-1 ;    
    if nnz(crftrials)==0
        continue
    end
    for ci = 1:Nneurons
        pRrfclassic(ci,rfi) = signrank(R(ci,crftrials)-sponFR(ci));
    end
end
pRrfinverse = NaN(Nneurons, Nrfs);
for rfi = 1:Nrfs
    irftrials = floor(trialorder/10000)==1 & floor(mod(trialorder, 1000) / 10)==rfi-1 ;
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
    'Pfried_rfclassic', 'Pfried_rfinverse', 'sigfmc_rfclassic', 'sigfmc_rfinverse', ...
    'pRrfclassic', 'pRrfinverse', 'pRFclassic', 'pRFinverse', ...
    'RFsigexclclassic', 'RFsigexclinverse', 'RFexclsigclassic', 'RFexclsiginverse'};
RFCIspin = struct();
for f = 1:numel(RFCIfields)
    RFCIspin.(RFCIfields{f}) = eval(RFCIfields{f});
end

end
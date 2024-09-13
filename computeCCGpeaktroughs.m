function [hpk, tpk, htr, ttr] = computeCCGpeaktroughs(CCG, CCGtli)
% CCG format: Nneurons * Nneurons * Ntimepoints
% CCGplus: take >=0 timepoints, and rearrange so that the first axis is timepont

flanktli = abs(CCGtli)>50 & abs(CCGtli)<=100;
CCGflankstd = squeeze( std(CCG(:,:,flanktli),0,3) );
CCGflankmean = squeeze( mean(CCG(:,:,flanktli),3) );

Nneu = size(CCG,2);
CCGtlirs = CCGtli(CCGtli>=0);
CCGrs = permute(CCG(:,:,CCGtli>=0), [3 1 2]);
dCCGrs = diff(CCGrs);

% t0ind = find(CCGtli>=0, 1,'first');
% if ~isequaln( squeeze(dCCGrs(1,:,:)), squeeze(CCG(:,:,t0ind+1)-CCG(:,:,t0ind)) )
%     error('sanity check failed')
% end

dCCGrssign = dCCGrs;
dCCGrssign(dCCGrs>0) = 1;
dCCGrssign(dCCGrs<0) = -1;
ddCCGrssign = diff(dCCGrssign);

% find first peak: zero crossing of dCCGrs, from plus to minus. has to have
% plus right before this index
pkthresh = CCGflankmean+CCGflankstd;
hpk = zeros(Nneu, Nneu);
tpk = NaN(Nneu, Nneu);
for t = size(ddCCGrssign,1):-1:1
    pairs_cross_at_t = squeeze(CCGrs(t+1,:,:))>pkthresh & squeeze(dCCGrs(t,:,:))>0 & squeeze(ddCCGrssign(t,:,:))<0;
    CCGrst = squeeze(CCGrs(t+1,:,:));
tpk( pairs_cross_at_t ) = CCGtlirs(t+1);
hpk( pairs_cross_at_t ) = CCGrst(pairs_cross_at_t);
end

% sanity check
% figure
% for ix = 1:24
%     pp=randperm(Nneu,2);
% subplot(4,6,ix)
% hold all
% plot(CCGtli, squeeze(CCG(pp(1), pp(2),:)))
% plot([CCGtli(1) CCGtli(end)], CCGflankmean(pp(1), pp(2))*[1 1], 'k--')
% plot([CCGtli(1) CCGtli(end)], (CCGflankmean(pp(1), pp(2))+CCGflankstd(pp(1), pp(2)))*[1 1], 'b:')
% plot([CCGtli(1) CCGtli(end)], (CCGflankmean(pp(1), pp(2))-CCGflankstd(pp(1), pp(2)))*[1 1], 'b:')
% yl = ylim;
% plot([0 0], yl,'k--')
% ylim(yl)
% plot(tpk(pp(1), pp(2)), hpk(pp(1), pp(2)), 'r*')
% xlim([-25 25])
% title(sprintf('presyn %d postsyn %d Tpeak %d',pp(1), pp(2), tpk(pp(1), pp(2))))
% end

% find first trough: zero crossing of dCCGrs, from minus to plus
trthresh = CCGflankmean-CCGflankstd;
htr = zeros(Nneu, Nneu);
ttr = NaN(Nneu, Nneu);
for t = size(ddCCGrssign,1):-1:1
    pairs_cross_at_t = squeeze(CCGrs(t+1,:,:))<trthresh & squeeze(dCCGrs(t,:,:))<0 & squeeze(ddCCGrssign(t,:,:))>0;
    CCGrst = squeeze(CCGrs(t+1,:,:));
ttr( pairs_cross_at_t ) = CCGtlirs(t+1);
htr( pairs_cross_at_t ) = CCGrst(pairs_cross_at_t);
end

% sanity check
% figure
% for ix = 1:24
%     pp=randperm(Nneu,2);
% subplot(4,6,ix)
% hold all
% plot(CCGtli, squeeze(CCG(pp(1), pp(2),:)))
% plot([CCGtli(1) CCGtli(end)], CCGflankmean(pp(1), pp(2))*[1 1], 'k--')
% plot([CCGtli(1) CCGtli(end)], (CCGflankmean(pp(1), pp(2))+CCGflankstd(pp(1), pp(2)))*[1 1], 'b:')
% plot([CCGtli(1) CCGtli(end)], (CCGflankmean(pp(1), pp(2))-CCGflankstd(pp(1), pp(2)))*[1 1], 'b:')
% yl = ylim;
% plot([0 0], yl,'k--')
% ylim(yl)
% plot(ttr(pp(1), pp(2)), htr(pp(1), pp(2)), 'r*')
% xlim([-25 25])
% title(sprintf('presyn %d postsyn %d Ttrough %d',pp(1), pp(2), ttr(pp(1), pp(2))))
% end

end
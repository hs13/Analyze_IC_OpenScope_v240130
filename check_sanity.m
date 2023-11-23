nwbsessions = {'sub-619296', 'sub-620333', 'sub-620334', 'sub-625545', ...
    'sub-625554', 'sub-625555',	'sub-630506', 'sub-630507', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};
pathpp = 'S:\OpenScopeData\00248_v230821\postprocessed\sub-637484\';
load([pathpp 'postprocessed_probeC.mat'])

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};

%% flat psth on blank trials
trialsoi = vis.ICwcfg0_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICwcfg0_presentations(:,trialsoi,:),[2,3]))

trialsoi = vis.ICwcfg1_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICwcfg1_presentations(:,trialsoi,:),[2,3]))




trialsoi = vis.ICwcfg0_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICwcfg0_presentations(:,trialsoi,:),[2,3]))

trialsoi = vis.ICwcfg1_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICwcfg1_presentations(:,trialsoi,:),[2,3]))

trialsoi = vis.ICkcfg0_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICkcfg0_presentations(:,trialsoi,:),[2,3]))

trialsoi = vis.ICkcfg1_presentations.trialorder==0;
figure; plot(psthtli, mean(psth.ICkcfg1_presentations(:,trialsoi,:),[2,3]))

figure; hold all; histogram(ICwcfg1_SP_BK_agg); histogram(ICwcfg0_SP_BK_agg)


%% correlation
figure
for b = 1:4
tempcorr = corr(Ronavgall.(visblocks{b})(neuinarea,:));
subplot(2,2,b)
imagesc(tempcorr)
set(gca, 'XTick', 1:size(tempcorr,2), 'XTickLabel', ICtrialtypes, 'YTick', 1:size(tempcorr,2), 'YTickLabel', ICtrialtypes)
title(visblocks{b})
end

figure; hold all
plot(Ronavgall.ICwcfg1_presentations(neuinarea,ICtrialtypes==106), Ronavgall.ICwcfg1_presentations(neuinarea,ICtrialtypes==111), 'o')
xl = xlim;
plot(xl, xl, 'r-')

figure; hold all
plot(Ronavgall.ICwcfg1_presentations(neuinarea,ICtrialtypes==106), Ronavgall.ICwcfg0_presentations(neuinarea,ICtrialtypes==106), 'o')
xl = xlim;
plot(xl, xl, 'r-')

%% psth over time
xt = -400:200:1500;
figure%('Position',[100 100 2000 800])
for ises = 1:numel(nwbsessions)
pathpp = ['S:\OpenScopeData\00248_v230821\postprocessed\' nwbsessions{ises} '\'];
load([pathpp 'postprocessed_probeC.mat'])
% annotation('textbox', [0.1 0.9 0.8 0.1], 'string', nwbsessions{ises}, 'edgecolor', 'none')
for ii = 1:4
subplot(4,numel(nwbsessions),numel(nwbsessions)*(ii-1)+ises)
imagesc( 1000*squeeze(mean(psth.(visblocks{ii}),3))' )
set(gca, 'XTick', find(ismember(psthtli,xt)), 'XTickLabel', xt)
colormap redblue
xlim(find(ismember(psthtli,[-200 600])))
%colorbar
caxis([0 15])
xlabel('Time (ms)')
ylabel('Trials')
title({nwbsessions{ises}, visblocks{ii}}, 'interpreter', 'none')
end
end

figure
for b = 1:4
Ntrials = size(psth.(visblocks{b}),2);
Nt2p = round(Ntrials/10);
subplot(2,2,b)
hold all
for ii = 1:10
    if ii==10
        trialsoi = Nt2p*(ii-1)+1:Ntrials;
    else
        trialsoi = Nt2p*(ii-1)+1:Nt2p*ii;
    end
    plot(psthtli, 1000*squeeze(mean(psth.ICwcfg1_presentations(:,trialsoi,:),[2 3]))', 'Color', [(ii-1)/9 0 0] )
end
xlim([-200 600])
ylim([4 9])
xlabel('Time (ms)')
ylabel('Rate (Hz)')
title(visblocks{b}, 'interpreter', 'none')
end

Ntrials = size(psth.ICwcfg1_presentations,2);
legs = cell(1,ceil(Ntrials/500));
figure
hold all
for ii = 1:ceil(Ntrials/500)
    if 500*ii<Ntrials
        trialsoi = 500*(ii-1)+1:500*ii;
    else
        trialsoi = 500*(ii-1)+1:Ntrials;
    end
    plot(psthtli, 1000*squeeze(mean(psth.ICwcfg1_presentations(:,trialsoi,:),[2 3]))', 'Color', [(ii-1)/(ceil(Ntrials/500) -1) 0 0] )
    legs{ii} = sprintf('Trials %d-%d', trialsoi(1), trialsoi(end));
end
legend(legs)
xlim([-200 600])
xlabel('Time (ms)')
ylabel('Rate (Hz)')
title([nwbsessions{ises}, ' ICwcfg1'], 'interpreter', 'none')


figure
hold all
Ntrials = size(psth.ICwcfg1_presentations,2);
tloi = psthtli>0 & psthtli<400;
for ii = 1:Ntrials/500
    if 500*ii<Ntrials
        trialsoi = 500*(ii-1)+1:500*ii;
    else
        trialsoi = 500*(ii-1)+1:Ntrials;
    end
    plot(ii, 1000*squeeze(mean(psth.ICwcfg1_presentations(tloi,trialsoi,:),'all')), 'o-', 'Color', [(ii-1)/(Ntrials/500 -1) 0 0] )
end

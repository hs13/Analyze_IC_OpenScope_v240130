%% PCA on LFP features:
% amplitude (range) for each channel
% area under mean-subtracted waveform (i.e., sum of absolute values) for each channel
% peak frequency component for each channel after whitening
% max source channel and timing and max sink channel and timing (or just time-averaged CSD vector)
% try multiple time windows (400ms, 200ms, 100ms, 50ms) 
% categorize sensory evoked activity?
% spike triggered average for each neuron?
% or just trigger based on suprathreshold LFP peaks on each channel (where threshold is set as 6x standard deviation in detrended LFP) -- set threshold so that we get ~2.5 events per second (range between 1 through 5 events per second)
% K-means clustering: decide number of clusters using AIC or BIC or Calinski-Harabasz Index (Jung...Kim Neuron)
% inspect CSD and TFR of clusters

addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('/Users/hyeyoung/Documents/CODE/SpectralEvents')) % for TFR

probes = {'A', 'B', 'C', 'D', 'E'};

pathpp = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/IllusionOpenScope/sub-620333/';
nwbspikefile = '';

load([pathpp 'LFP_CSD_probeC.mat'])
load([pathpp 'LFP_psth_probeC.mat'])
load([pathpp 'LFP_TFR_L23_probeC.mat'])
load([pathpp 'LFP1000Hz_probeC.mat'])

lfpstart = floor((vis.ICwcfg1_presentations.start_time(1)-lfptimeresamp(1))/Tres)+1;
lfpend = floor((vis.ICwcfg1_presentations.stop_time(end)-lfptimeresamp(1))/Tres)+1;

tic
lfpblock = lfpresamp(:,lfpstart:lfpend);

%% CSD: negative values mean sink (influx, depol), positive values mean source (outflux)
% lfpelecspacing = 0.04; % 40micrometers, i.e., 0.04mm
% csdelectinds = 2:Nelec-1;
% csdresamp = -( lfpresamp(csdelectinds+1,:)-2*lfpresamp(csdelectinds,:)+lfpresamp(csdelectinds-1,:) )/(lfpelecspacing.^2);

% smooth with a gaussian kernel first
kerwinhalf = 5; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));
lfpconv = convn(lfpblock, kergauss, 'same');

Nelec = numel(lfpelecid);
lfpelecspacing = 0.04; % 40micrometers, i.e., 0.04mm
csdelectinds = 2:Nelec-1;
csdblock = -( lfpconv(csdelectinds+1,:)-2*lfpconv(csdelectinds,:)+lfpconv(csdelectinds-1,:) )/(lfpelecspacing.^2);

%% TFR : for now, just choose one electrode in layer 2/3
load(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{iprobe}), 'elecL23')
%elecL23 = round(median( find(contains(lfpelecvec.location, '2/3')) ));
electrodesL23 = find(contains(lfpelecvec.location(csdelectinds), '2/3'));

% 30 seconds and 7 GB for each electrode
S = lfpblock(elecL23,:)';
% get rid of NaN values in S
[TFR_L23_block,tVec,fVec] = spectralevents_ts2tfr(S,1:100,1/Tres,5);
toc

%% z-score lfp, csd, tfr
tic
lfpzscore = (lfpblock-mean(lfpblock,'all'))/std(lfpblock,0,'all');
csdzscore = csdblock/std(csdblock,0,'all');
tfrzscore = (TFR_L23_block-mean(TFR_L23_block,'all'))/std(TFR_L23_block,0,'all');
toc

%% find LFP timings that are more than X standard deviations above or below the mean
% 3 std: ~5% peaks and sub threshold each, meaning every 10 ms or so
% 7 std: 25 peaks and 101 troughs in a 72min recording
% 6 std: 100 peaks and 564 troughs in a 72min recording: every 7.6 seconds
maximaopt = 'CSD';
switch maximaopt
    case 'LFP'
        lfpmat = lfpzscore;
    case 'CSD'
        lfpmat = csdzscore;
    case 'TFR'
        lfpmat = tfrzscore;
end

islfppeak = false(size(lfpmat));
islfptrough = false(size(lfpmat));
tic
for e = 1:size(lfpmat,1)
    [pks,locs] = findpeaks(lfpmat(e,:));
    islfppeak(e,locs) = true;
    
    [pks,locs] = findpeaks(-lfpmat(e,:));
    islfptrough(e,locs) = true;
end
toc

Xstd = 5; % factor of standard deviation
tpeaks = any(islfppeak & lfpmat > mean(lfpmat,2)+Xstd*std(lfpmat,0,2), 1);
ttroughs = any(islfptrough & lfpmat < mean(lfpmat,2)-Xstd*std(lfpmat,0,2), 1);
fprintf('%d*std threshold %s: %d troughs interval %.3fs, %d peaks interval %.3fs\n', ...
    Xstd, maximaopt, nnz(ttroughs), Tres*mean(diff(find(ttroughs))), nnz(tpeaks), Tres*mean(diff(find(tpeaks))) )
%4*std threshold LFP: 17519 troughs interval 0.245s, 12586 peaks interval 0.341s
%5*std threshold LFP: 2644 troughs interval 1.621s, 801 peaks interval 5.348s
%6*std threshold LFP: 332 troughs interval 12.7s, 74 peaks interval 54.5s

%% +/-250ms of troughs/peaks
%tmaximas = ttroughs | tpeaks;
tmaximas = ttroughs;

Twin = 250;
trange = -Twin:Twin;
temptinds = find(tmaximas)'+trange;
tic
lfpmaximas = NaN(nnz(tmaximas), length(trange), Nelec);
for e = 1:numel(lfpelecid)
    templfp = lfpblock(e,:);
    lfpmaximas(:,:,e) = templfp(temptinds);
end

csdmaximas = NaN(nnz(tmaximas), length(trange), numel(csdelectinds) );
for e = 1:numel(csdelectinds)
    tempcsd = csdblock(e,:);
    csdmaximas(:,:,e) = tempcsd(temptinds);
end

tfrmaximas = NaN(nnz(tmaximas), length(trange), numel(fVec));
for e = 1:numel(fVec)
    temptfr = TFR_L23_block(e,:);
    tfrmaximas(:,:,e) = temptfr(temptinds);
end
toc


%% PCA on troughs/peaks
lfpzvec = lfpzscore(:,tmaximas)';
csdzvec = csdzscore(:,tmaximas)';
foi = fVec>=7 & fVec <= 80;
tfrzvec = tfrzscore(foi,tmaximas)';

[coeff, score, latent, tsquared, explained] = pca([csdzvec tfrzvec lfpzvec]);
% [coeff, score, latent, tsquared, explained] = pca(csdzvec);

npc = 4;
figure
for ipc = 1:npc
    for jpc = 1:npc
        subplot(npc,npc,npc*(jpc-1)+ipc)
plot(score(:,ipc),score(:,jpc), '.')
xlabel(['PC' num2str(ipc)])
ylabel(['PC' num2str(jpc)])
    end
end


%% how many clusters?
nclusters = 3;
ipc = 1; jpc = 2;

figure
for npc2k=1:12
clustercol = jet(nclusters);
IDX = kmeans(score(:,1:npc2k), nclusters);

subplot(3,4,npc2k)
hold all
%plot(score(:,ipc),score(:,jpc), 'o')
for c = 1:nclusters
plot(score(IDX==c,ipc),score(IDX==c,jpc), 'x', 'Color', clustercol(c,:))
end
xlabel(['PC' num2str(ipc)])
ylabel(['PC' num2str(jpc)])
title(npc2k)
end
%%
npc2k=3;
IDX = kmeans(score(:,1:npc2k), nclusters);

clustercol = jet(nclusters);
figure
hold all
%plot(score(:,ipc),score(:,jpc), 'o')
for c = 1:nclusters
plot(score(IDX==c,ipc),score(IDX==c,jpc), 'x', 'Color', clustercol(c,:))
end
xlabel(['PC' num2str(ipc)])
ylabel(['PC' num2str(jpc)])
legend

%% plot average LFP, CSD and TFR of clusters
yl = [0.5 Nelec+0.5];
figure
for c = 1:nclusters
    trinds = IDX==c;
    subplot(nclusters, 3, (c-1)*3+1)
    imagesc(trange, 1:Nelec, squeeze(mean(lfpmaximas(trinds,:,:),1))' )
    ylim(yl)
    set(gca, 'XGrid', 'on', 'YTick', 1:Nelec, 'YTickLabel', lfpelecvec.location, 'YDir', 'normal')
    colorbar
    caxis(prctile(lfpmaximas(:), [2.5 97.5]));
    title(sprintf('Cluster%d LFP', c))
    
    
    subplot(nclusters, 3, (c-1)*3+2)
    imagesc(trange, csdelectinds, squeeze(mean(csdmaximas(trinds,:,:),1))' )
    ylim(yl)
    set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
    colorbar
    caxis(prctile(csdmaximas(:), [2.5 97.5]));
    title(sprintf('Cluster%d CSD', c))
    
    subplot(nclusters, 3, c*3)
    imagesc(trange, fVec, squeeze(mean(tfrmaximas(trinds,:,:),1))' );
    ylim([7 80]);
    caxis(prctile(tfrmaximas(:), [2.5 97.5]));
    colorbar
    title(sprintf('Cluster%d TFR', c))
end
colormap jet


ctxelec = contains(lfpelecvec.location, 'VIS');
ctxelectop = find(ctxelec, 1, 'last');
ctxelecbottom = find(ctxelec, 1, 'first');
yl = [ctxelecbottom ctxelectop]+.5;
figure
for c = 1:nclusters
    trinds = IDX==c;
    subplot(nclusters, 3, (c-1)*3+1)
    imagesc(trange, 1:Nelec, squeeze(mean(lfpmaximas(trinds,:,:),1))' )
    ylim(yl)
    set(gca, 'XGrid', 'on', 'YTick', 1:Nelec, 'YTickLabel', lfpelecvec.location, 'YDir', 'normal')
    colorbar
    caxis(prctile(lfpmaximas(:), [2.5 97.5]));
    title(sprintf('Cluster%d LFP', c))
    
    
    subplot(nclusters, 3, (c-1)*3+2)
    imagesc(trange, csdelectinds, squeeze(mean(csdmaximas(trinds,:,:),1))' )
    ylim(yl)
    set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
    colorbar
    caxis(prctile(csdmaximas(:), [2.5 97.5]));
    title(sprintf('Cluster%d CSD', c))
    
    subplot(nclusters, 3, c*3)
    imagesc(trange, fVec, squeeze(mean(tfrmaximas(trinds,:,:),1))' );
    ylim([7 80]);
    caxis(prctile(tfrmaximas(:), [2.5 97.5]));
    caxis(prctile(tfrmaximas(:), [0.5 99.5]));
    colorbar
    title(sprintf('Cluster%d TFR', c))
end
colormap jet

figure
imagesc(trange, fVec, squeeze(mean(tfrmaximas(trinds,:,:),1))' );
ylim([7 80]);
%     caxis(prctile(tfrmaximas(:), [2.5 97.5]));
colorbar
title(sprintf('Cluster%d TFR', c))

figure;plot(fVec, squeeze(mean(tfrmaximas(trinds,trange==0,:),1)) )
figure;hold all
plot(fVec, squeeze(tfrmaximas(trinds,trange==0,:)) )
plot(fVec, squeeze(mean(tfrmaximas(trinds,trange==0,:),1)),'k-','linewidth', 3 )
plot(fVec, median(TFR_L23_block,2),'r-','linewidth', 2 )

tfrmedvec = median(TFR_L23_block,2);
normtfrmaximas = tfrmaximas./reshape(tfrmedvec,1,1,length(fVec));
figure;hold all
plot(fVec, squeeze(normtfrmaximas(trinds,trange==0,:)) )
plot(fVec, squeeze(mean(normtfrmaximas(trinds,trange==0,:),1)),'k-','linewidth', 3 )


[mv,mi]= max(squeeze(normtfrmaximas(trinds,trange==0,foi)), [],2);
[sv,si]=sort(mi);

clustinds = find(trinds);
figure
imagesc( squeeze(tfrmaximas(clustinds(si),trange==0,:)) )
xlim([7 80]);

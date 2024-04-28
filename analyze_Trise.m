% find the first spike timing that exceeds 97.5 percentile
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
    clearvars -except ises Nsessions datadir nwbsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    sesclk = tic;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];

    probes = {'A', 'B', 'C', 'D', 'E', 'F'};
    whichblock = 'ICwcfg1_presentations';

    load([pathpp 'postprocessed_probeC.mat'], 'psthtli', 'vis')
    Twin = 5;
    blanktrialinds = find( vis.(whichblock).trialorder==0 );
    Nblank = numel(blanktrialinds);

    for iprobe = 1:numel(probes)
        clearvars neuoind psth
        Trisefn = [pathpp 'Trise' num2str(Twin) 'msbins_probe' probes{iprobe} '.mat'];
        if exist(Trisefn, 'file')
            disp([Trisefn ' exists'])
            continue
        end
        
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        Nneurons = numel(neuoind);
        psthblank = 1000*psth.(whichblock)(:,blanktrialinds,:);
        psthblankconv = convn(psthblank, ones(Twin,1)/Twin, 'same');

        % 5 min
        tic
        Nboot = 1000;
        psthblankconvboot = NaN(length(psthtli), Nboot, Nneurons);
        for b = 1:Nboot
            temptrials = randi(Nblank, Nblank,1);
            psthblankconvboot(:,b,:) = mean( psthblankconv(:,temptrials,:),2);
        end
        toc

        tliim = find(psthtli>=0 & psthtli<400);
        psthtliim = psthtli(tliim);
        [mv,mi] = max(psthblankconvboot(tliim,:,:),[],1);
        psthblankconvpeak = squeeze(mv);
        psthblankconvTpeak = squeeze(psthtliim(mi));

        [mv,mi]=max( psthblankconvboot(tliim,:,:)-psthblankconvboot(tliim-Twin,:,:),[],1 );
        blankmaxdrate = squeeze(mv);
        blankTmaxrise = squeeze(psthtliim(mi));

        vistrialtypes = unique(vis.(whichblock).trialorder);
        Ntt = numel(vistrialtypes);
        vistrialrep = zeros(Ntt,1);
        psthconvavg = NaN(length(psthtli), Ntt, Nneurons);
        psthconvprct = NaN(length(psthtli), Ntt, Nneurons);
        psthconvglobprct = NaN(length(psthtli), Ntt, Nneurons);
        temppeakrep = repmat( reshape(psthblankconvpeak,1,Nboot,Nneurons),length(psthtli),1);
        tic
        for ii = 1:Ntt
            trialsoi = vis.(whichblock).trialorder==vistrialtypes(ii);
            vistrialrep(ii) = nnz(trialsoi);
            temppsthconv = convn(1000*psth.(whichblock)(:,trialsoi,:), ones(Twin,1)/Twin, 'same');
            psthconvavg(:,ii,:) = mean(temppsthconv, 2);

            psthconvprct(:,ii,:) = mean( psthblankconvboot<mean(temppsthconv, 2), 2);
            psthconvglobprct(:,ii,:) = mean( temppeakrep<mean(temppsthconv, 2), 2);
        end
        toc

        [mv,mi] = max(psthconvavg(tliim,:,:),[],1);
        psthconvavgpeak = squeeze(mv);
        psthconvavgTpeak = squeeze(psthtliim(mi));
        psthconvavgpeakprct = NaN(Ntt, Nneurons);
        for ii = 1:Ntt
            psthconvavgpeakprct(ii,:) = mean(psthblankconvpeak<psthconvavgpeak(ii,:), 1);
        end
        % Tpeakind = sub2ind(size(psthconvglobprct), squeeze(tliim(mi)), repmat((1:Ntt)',1,Nneurons), repmat(1:Nneurons,Ntt,1) );
        % isequaln(psthconvglobprct(Tpeakind), psthconvavgpeakprct) % sanity check

        [mv,mi]=max( psthconvavg(tliim,:,:)-psthconvavg(tliim-Twin,:,:),[],1 );
        maxdrate = squeeze(mv);
        Tmaxrise = squeeze(psthtliim(mi));
        maxdrateprct = NaN(Ntt, Nneurons);
        for ii = 1:Ntt
            maxdrateprct(ii,:) = mean(blankmaxdrate<maxdrate(ii,:), 1);
        end

        save([pathpp 'Trise' num2str(Twin) 'msbins_probe' probes{iprobe} '.mat'], ...
            'psthtli', 'tliim', 'psthtliim', ...
            'psthblankconvpeak', 'psthblankconvTpeak', 'blankmaxdrate', 'blankTmaxrise', ...
            'vistrialtypes', 'vistrialrep', 'psthconvavg', 'psthconvprct', 'psthconvglobprct', ...
            'psthconvavgpeak', 'psthconvavgTpeak', 'psthconvavgpeakprct', 'maxdrate', 'Tmaxrise', 'maxdrateprct')
    end
    toc(sesclk)
end

%% test plots
ci = 15; ii = find(ICtrialtypes ==1109);
figure; 
for ii = 1:20
    ci = 7*(ii+10);
    subplot(4,5,ii)
yyaxis left
hold all
plot(psthtli, psthconvavg(:,ii,ci))
plot(psthtli, psthconvavg(:,1,ci), 'k-')
g97pt5 = prctile(psthblankconvpeak,97.5,1);
g2pt5 = prctile(psthblankconvpeak,2.5,1);
plot([psthtli(1) psthtli(end)], g97pt5(ci)*[1 1], 'r-')
% plot([psthtli(1) psthtli(end)], g2pt5(ci)*[1 1], 'r-')
% yyaxis right
% plot(psthtli, psthconvprct(:,ii,ci))
xlim([0 400])
end

ci = 140;
figure
hold all
plot(psthtli, psthconvavg(:,ii,ci))
plot(psthtli, psthconvavg(:,1,ci), 'k-')
plot(psthtli, squeeze(psthblankconvboot(:,:,ci)))
g97pt5 = prctile(psthblankconvpeak,97.5,1);
plot([psthtli(1) psthtli(end)], g97pt5(ci)*[1 1], 'r-')
xlim([0 400])


load([pathpp 'postprocessed.mat'])
load([pathpp 'postprocessed_probeC.mat'], 'psthtli', 'vis')

Nneuall = length(neuallloc);
psthall = false(length(psthtli), length(vis.(whichblock).trialorder), Nneuall);
for iprobe = 1:numel(probes)
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
    psthall(:,:,neuoind) = psth.(whichblock);
end

% neuctx = contains(neuallloc, 'VIS');
% psthctx = psthall(:,:,neuctx);

% if I compute psthblankconv bootstrapped 1000X, the resulting variable for
% one session will be ~25gb
Twin = 5;
blanktrialinds = find( vis.(whichblock).trialorder==0 );
psthblank = psthall(:,blanktrialinds,:);
psthblankconv = convn(psthblank, ones(Twin,1)/Twin, 'same');

% this takes 70 min -- it's probably better to calculate for
% every probe separately...
tic
Nblank = numel(blanktrialinds);
Nboot = 1000;
psthblankconvboot = NaN(length(psthtli), Nboot, Nneuall);
for b = 1:Nboot
    temptrials = randi(Nblank, Nblank,1);
    psthblankconvboot(:,b,:) = mean( psthblankconv(:,temptrials,:),2);
end
toc

% this takes 72 seconds
tliim = find(psthtli>=0 & psthtli<400);
vistrialtypes = unique(vis.(whichblock).trialorder);
Ntt = numel(vistrialtypes);
vistrialrep = zeros(Ntt,1);
psthconvavg = NaN(length(psthtli), Ntt, Nneuall);
psthconvprct = NaN(length(psthtli), Ntt, Nneuall);
tic
for ii = 1:Ntt
    trialsoi = vis.(whichblock).trialorder==vistrialtypes(ii);
    vistrialrep(ii) = nnz(trialsoi);
    temppsthconv = convn(psthall(:,trialsoi,:), ones(Twin,1)/Twin, 'same');
    psthconvavg(:,ii,:) = mean(1000*temppsthconv, 2);
end
toc

tliim = find(psthtli>=0 & psthtli<400);
psthtliim = psthtli(tliim);
[mv,mi]=max( psthconvavg(tliim,:,:)-psthconvavg(tliim-5,:,:),[],1 );
maxdeltarate = squeeze(mv);
Tmaxrise = squeeze(psthtliim(mi));

ci = 488;
neuinarea = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
cgroup = find(neuinarea & ICsigall.ICwcfg1_presentations.ICencoder1);
cgroup = find(neuinarea);
c2p = cgroup(randperm(numel(cgroup),20));
tt2p = 1109;
psthblankconv97pt5 = prctile(psthblankconvboot,97.5,2);
psthblankconv2pt5 = prctile(psthblankconvboot,2.5,2);
figure
for ii = 1:20
    ci = c2p(ii);
    subplot(4,5,ii)
hold all
tempmean = squeeze(mean(psthblankconv(:,:,ci),2));
temp97pt5 = squeeze(psthblankconv97pt5(:,:,ci));
temp2pt5 = squeeze(psthblankconv2pt5(:,:,ci));
shadedErrorBar(psthtli, 1000*tempmean, 1000*[temp97pt5-tempmean tempmean-temp2pt5], 'k-')
plot(psthtli, squeeze( psthconvavg(:,ICtrialtypes==tt2p,ci)) )
% tp = Tmaxrise(ICtrialtypes==tt2p,ci);
% plot( tp, squeeze( psthconvavg(psthtli==tp,ICtrialtypes==tt2p,ci)), 'r*' )
xlim([0 400])
end

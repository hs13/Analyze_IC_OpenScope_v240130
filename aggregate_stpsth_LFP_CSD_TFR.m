datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

%% display V1 probe
for ises = 1:numel(nwbsessions)
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat'])
    V1probes = electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
    disp(nwbsessions{ises})
    disp(unique(V1probes))
end

%% spike triggered CSD during visual presentation period of each trial type
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
for ises = 1:numel(nwbsessions)
    clearvars -except ises datadir nwbsessions probes
    sesclk = tic;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat'])
    V1probes = electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
    if numel(unique(V1probes))>1
        error('this session has more than one V1 probe')
    end
    iprobe = mode(V1probes)+1;

    %for iprobe = 1:numel(probes)
    fprintf('%d/%d %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, probes{iprobe})
    if ~exist(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{iprobe}), 'file')
        fprintf('LFP_CSD_probe%s.mat does not exit!!!\n', probes{iprobe})
        continue
    end
    load(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{iprobe}))
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))

    whichblock = 'ICwcfg1_presentations';
    blocksplit = strsplit(whichblock);
    blockname = blocksplit{1};
    if isfield(vis.(whichblock), 'ICtrialtypes')
        trialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
    else
        trialorder = vis.(whichblock).trialorder;
    end
    vistrialtypes = unique(trialorder);
    Nneu = size(psth.(whichblock), 3);
    tloi = psthtli>=0 & psthtli<400;

    %sttrange = -200:200; % 160s for each trialtype
    sttrange = -100:100; % 80s for each trialtype
    stCSDvisprobe = cell(size(vistrialtypes));
    for typi = 1:numel(vistrialtypes)
        tic
        stCSDvisprobe{typi} = NaN(Nneu, length(sttrange), numel(csdelectinds));
        trialsoi = trialorder==vistrialtypes(typi);
        for t = 1:length(sttrange)
            temptloind = find(tloi)+sttrange(t);
            tempcsdpsth = reshape( csdvispsth.(whichblock)(temptloind, trialsoi, :), nnz(tloi)*nnz(trialsoi), numel(csdelectinds));
            for ci = 1:Nneu
                temppsth = reshape( psth.(whichblock)(tloi, trialsoi, ci), [],1 );
                stCSDvisprobe{typi}(ci,t,:) = mean(tempcsdpsth(temppsth,:),1);
            end
        end
        toc
    end
    %end
    save( sprintf('%sV1_stCSD%s_probe%s.mat', pathpp, blockname, probes{iprobe}), ...
        'lfpelecspacing', 'csdelectinds', 'vistrialtypes', 'sttrange', 'stCSDvisprobe', '-v7.3')
    toc(sesclk)
end

%%
load(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{iprobe}), 'lfpelecvec')
ctxelec = contains(lfpelecvec.location, 'VIS');
ctxelectop = find(ctxelec, 1, 'last');
ctxelecbottom = find(ctxelec, 1, 'first');
yl = [ctxelecbottom ctxelectop]+.5;
lfpeleclocation = lfpelecvec.location;
neurandord = randperm(Nneu);
figure; 
%imagesc(trange, csdelectinds, squeeze(nanmean(stCSDvisprobe{typi},1))')
for ii = 1:32
    ci = neurandord(ii);
    subplot(4,8,ii)
imagesc(sttrange, csdelectinds, squeeze(stCSDvisprobe{typi}(ci,:,:))')
set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpeleclocation(csdelectinds), 'YDir', 'normal')
caxis([-0.015 0.015])
ylim(yl)
xlim([-100 100])
end

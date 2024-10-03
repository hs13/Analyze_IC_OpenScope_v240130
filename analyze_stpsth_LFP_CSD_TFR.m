datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visctxareas = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};

lowpassopt = false;
whichneuarea = 'V1';
lfpareas = {'V1', 'LM'};%
whichblock = 'ICwcfg1_presentations';

%% display V1 probe
for ises = 1:numel(nwbsessions)
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat'])
    disp(nwbsessions{ises})
    for a = 1:numel(lfpareas)
        if strcmp(lfpareas{a}, 'V1')
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
        else
            lfparea = visctxareas{strcmp(visareas, lfpareas{a})};
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, lfparea) );
        end
        disp(lfpareas{a})
        if numel( unique(lfpareaprobeinds) ) ==1
            disp(reshape( unique(lfpareaprobeinds), 1,[]))
        else
            warning('expected one probe for one area')
            [v,c] = uniquecnt(lfpareaprobeinds);
            disp(reshape(v,1,[]))
            disp(reshape(c,1,[]))
        end
    end
end

%% spike triggered CSD during visual presentation period of each trial type
for ises = 1:numel(nwbsessions)
    clearvars -except ises datadir nwbsessions probes visareas visctxareas ...
        lowpassopt whichneuarea lfpareas whichblock



    sesclk = tic;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat'])

    if strcmp(whichneuarea, 'V1')
        neuareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
    else
        neuarea = visctxareas{strcmp(visareas, whichneuarea)};
        neuareaprobeinds = 1+electrode_probeid( contains(electrode_location, neuarea) );
    end
    if numel(unique(neuareaprobeinds))>1
        warning('this session has more than one %s probe', neuarea)
    end
    neuprobe = mode(neuareaprobeinds);
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{neuprobe}))

    blocksplit = strsplit(whichblock);
    blockname = blocksplit{1};

    for a = 1:numel(lfpareas)
        lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
        if strcmp(lfpareas{a}, 'V1')
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
        else
            lfparea = visctxareas{strcmp(visareas, lfpareas{a})};
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, lfparea) );
        end
        if numel(unique(lfpareaprobeinds))>1
            warning('this session has more than one %s probe', lfpareas{a})
        end
        lfpprobe = mode(lfpareaprobeinds);

        %for iprobe = 1:numel(probes)
        fprintf('%d/%d %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, probes{lfpprobe})
        if ~exist(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{lfpprobe}), 'file')
            fprintf('LFP_CSD_probe%s.mat does not exit!!!\n', probes{lfpprobe})
            continue
        end
        if lowpassopt
            stpsthCSDfn = sprintf('%s%s_stCSD%s_lowpassprobe%s.mat', pathpp, whichneuarea, blockname, probes{lfpprobe});
        else
            stpsthCSDfn = sprintf('%s%s_stCSD%s_probe%s.mat', pathpp, whichneuarea, blockname, probes{lfpprobe});
        end
        if exist(stpsthCSDfn, 'file')
            fprintf('%s exists\n', stpsthCSDfn)
            continue
        end

        if lowpassopt
            load(sprintf('%sLFP_psth_probe%s.mat', pathpp, probes{lfpprobe}), 'lfpvispsth')
            Nelec = size(lfpvispsth.(whichblock),3);
            csdelectinds = 2:Nelec-1;
            lfpelecspacing = 0.04; % 40micrometers, i.e., 0.04mm

            lfplpblockpsth = NaN(size(lfpvispsth.(whichblock)));
            for e = 1:Nelec
                lfplpblockpsth(:,:,e) = lowpass( squeeze(lfpvispsth.(whichblock)(:,:,e)), 100,1000);
            end

            csdblockpsth = -( lfplpblockpsth(:,:,csdelectinds+1) - 2*lfplpblockpsth(:,:,csdelectinds) + ...
                lfplpblockpsth(:,:,csdelectinds-1) )/(lfpelecspacing.^2);
        else
            load(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{lfpprobe}), 'csdelectinds', 'csdvispsth', 'lfpelecspacing')
            csdblockpsth = csdvispsth.(whichblock);
        end
        if exist(stpsthCSDfn, 'file')
            fprintf('%s exists\n', stpsthCSDfn)
            continue
        end

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
                tempcsdpsth = reshape( csdblockpsth(temptloind, trialsoi, :), nnz(tloi)*nnz(trialsoi), numel(csdelectinds));
                for ci = 1:Nneu
                    temppsth = reshape( psth.(whichblock)(tloi, trialsoi, ci), [],1 );
                    stCSDvisprobe{typi}(ci,t,:) = mean(tempcsdpsth(temppsth,:),1);
                end
            end
            toc
        end
        %end
        save( stpsthCSDfn, 'lfpelecspacing', 'csdelectinds', 'vistrialtypes', 'sttrange', 'stCSDvisprobe', '-v7.3')
    end
    toc(sesclk)
end

%%
load(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{lfpprobe}), 'lfpelecvec')
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

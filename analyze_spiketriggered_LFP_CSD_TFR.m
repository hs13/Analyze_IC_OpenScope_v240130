clear all; close all; %%clc

%addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\SpectralEvents')) % for TFR

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E'};

for ises = 1:numel(nwbsessions)
    sesclk = tic;
    clearvars -except ises datadir nwbsessions probes
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));

    nwblfpfiles = nwbfiles(contains({nwbfiles.name}, 'probe'));
    if isempty(nwblfpfiles)
        fprintf('%d/%d %s does not have LFP recordings\n', ises, numel(nwbsessions), nwbsessions{ises})
        continue
    end

    % take filename  with shortest length or filename that does not contain probe
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    % nwbspikefile = string(nwbspikefile);
    disp(nwbspikefile)
    nwb = nwbRead(nwbspikefile); %, 'ignorecache');

    %% extract spike times
    unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this session
    unit_peakch = nwb.units.vectordata.get('peak_channel_id').data.load();
    unit_times_data = nwb.units.spike_times.data.load();
    unit_times_idx = nwb.units.spike_times_index.data.load();
    % unit_waveform = nwb.units.waveform_mean.data.load();
    unit_wfdur = nwb.units.vectordata.get('waveform_duration').data.load();

    Nneurons = length(unit_ids);

    % all(ismember(unit_peakch, electrode_id))

    spiketimes = cell(Nneurons, 1);
    last_idx = 0;
    for ii = 1:Nneurons
        unit_id = unit_ids(ii);

        %     assert(unit_trials_idx(i) == unit_times_idx(i), 'Expected unit boundaries to match between trials & spike_times jagged arrays')
        start_idx = last_idx + 1;
        end_idx = unit_times_idx(ii);

        spiketimes{ii} = unit_times_data(start_idx:end_idx);

        last_idx = end_idx;
    end

    load( [pathpp 'LFP1000Hz_probeC.mat'], 'lfptimeresamp' )
    load( [pathpp 'LFP_psth_probeC.mat'], 'Tres', 'vis' )
    stlen = round((lfptimeresamp(end))/Tres);
    spiketrain = false(stlen, Nneurons);
    for ii = 1:Nneurons
        spiketrain(floor(spiketimes{ii}/Tres)+1, ii) = true;
    end
    ststart = round((vis.ICwcfg1_presentations.start_time(1))/Tres);
    stend = round((vis.ICwcfg1_presentations.stop_time(end))/Tres);
    spiketrain = spiketrain(ststart:stend,:);

    lfpstart = floor((vis.ICwcfg1_presentations.start_time(1)-lfptimeresamp(1))/Tres)+1;
    lfpend = floor((vis.ICwcfg1_presentations.stop_time(end)-lfptimeresamp(1))/Tres)+1;
    if stend-ststart ~= lfpend-lfpstart
        error('mismatch between spike train length and lfp length')
    end


    %% spike triggered CSD and TFR
    for iprobe = 1:numel(probes)
        fprintf('%d/%d %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, probes{iprobe})
        probeclk = tic;
        load( sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{iprobe}) )
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}), 'neuoind')

        % 'lfpelecid', 'lfpelecvec', 'lfptimeresamp', 'lfpresamp'
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

        %% spike triggered CSD and TFR
        Twin = 250;
        trange = -Twin:Twin;
        spiketraintrimmed = spiketrain;
        spiketraintrimmed(1:Twin,:)=false;
        spiketraintrimmed(end-Twin+1:end,:)=false;
        stLFPprobe = NaN(Nneurons, length(trange), numel(lfpelecid));
        stCSDprobe = NaN(Nneurons, length(trange), numel(csdelectinds));
        stTFRprobe = NaN(Nneurons, length(trange), numel(fVec));
        for ii = 1:Nneurons
            tempstinds = find(spiketraintrimmed(:,ii))+trange;
            tempNspikes = size(tempstinds,1);
            tempstlfp = NaN( tempNspikes, length(trange), numel(lfpelecid) );
            for e = 1:numel(lfpelecid)
                templfp = lfpblock(e,:);
                tempstlfp(:,:,e) = templfp(tempstinds);
            end

            tempstcsd = NaN( tempNspikes, length(trange), numel(csdelectinds) );
            for e = 1:numel(csdelectinds)
                tempcsd = csdblock(e,:);
                tempstcsd(:,:,e) = tempcsd(tempstinds);
            end

            tempsttfr = NaN( tempNspikes, length(trange), numel(fVec) );
            for e = 1:numel(fVec)
                temptfr = TFR_L23_block(e,:);
                tempsttfr(:,:,e) = temptfr(tempstinds);
            end

            stLFPprobe(ii,:,:) = mean(tempstlfp,1);
            stCSDprobe(ii,:,:) = mean(tempstcsd,1);
            stTFRprobe(ii,:,:) = mean(tempsttfr,1);
        end
        save( sprintf('%sspiketriggered_LFP_CSD_TFR_probe%s.mat', pathpp, probes{iprobe}), ...
            'stLFPprobe', 'stCSDprobe', 'stTFRprobe', '-v7.3')
        toc(probeclk)
    end

    toc(sesclk)
end

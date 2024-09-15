clear all; close all; %%clc

%addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\SpectralEvents')) % for TFR

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E'};

for ises = 3:numel(nwbsessions)
    clearvars -except ises datadir nwbsessions probes
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));
    nwblfpfiles = nwbfiles(contains({nwbfiles.name}, 'probe'));
    if isempty(nwblfpfiles)
        fprintf('%d/%d %s does not have LFP recordings\n', ises, numel(nwbsessions), nwbsessions{ises})
        continue
    end

    %% load lfp data for each probe
    for iprobe = 1:numel(probes)
        probeclk = tic;
        fprintf('%d/%d %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, probes{iprobe})
        if exist(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{iprobe}), 'file')
            fprintf('LFP1000Hz_probe%s.mat already exits\n', probes{iprobe})
            continue
        end

        probeind = iprobe-1;
        fileind = find(contains({nwbfiles.name}, ['probe-' num2str(probeind)]));
        nwblfpfile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);

        nwblfp = nwbRead(nwblfpfile);

        % lfp = nwblfp.acquisition.get('probe_2_lfp_data');
        % lfpdata = lfp.data.load();

        lfp_probe = nwblfp.acquisition.get( ['probe_' num2str(probeind) '_lfp'] );
        lfp = lfp_probe.electricalseries.get( ['probe_' num2str(probeind) '_lfp_data'] );
        lfpelectrodes = lfp.electrodes.data.load();
        lfpdata = lfp.data.load(); % nprobes * ntimepoints
        lfptimestamps = ( lfp.timestamps.load() )'; % 1*ntimepoints, units in seconds

        %% figure out which electrode is where
        % Neuropixels 1.0 had 20um vertical spacing, Neuropixels 2.0 has 15um spacing
        % Neuropixels 1.0 was used in this dataset: vertical spacing 20um,
        % horizontal spacing 16um, closest electrodes are ~25um apart (sqrt(16^2+20^2))
        % relevant discussion here: https://community.brain-map.org/t/least-distance-between-2-recording-sites-on-neuropixels/1223/2
        % LFP: every four electrodes, i.e., 40 micrometer spacing

        lfpelecid = nwblfp.general_extracellular_ephys_electrodes.id.data.load();

        lfpelecvec = struct();
        lfpelecvecfields = nwblfp.general_extracellular_ephys_electrodes.colnames;
        for c= 1:numel(lfpelecvecfields)
            try
                lfpelecvec.(lfpelecvecfields{c}) = nwblfp.general_extracellular_ephys_electrodes.vectordata.get(lfpelecvecfields{c}).data.load();
            catch
                lfpelecvec.(lfpelecvecfields{c}) = nwblfp.general_extracellular_ephys_electrodes.vectordata.get(lfpelecvecfields{c}).data;%.load();
            end
        end

        % check electrode location (brain area) correspondence between lfp file and
        % spike file
        lfpeleclocation = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();
        load([pathpp 'info_electrodes.mat'])
        lfpelecind = zeros(size(lfpelecid));
        for e = 1:numel(lfpelecid)
            lfpelecind(e) = find(electrode_id==lfpelecid(e));
        end
        if ~isequal(electrode_location(lfpelecind), lfpeleclocation)
            error('electrode info on lfp nwb file does not match spike nwb file')
        end

        %
        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        if ~isequal(ctxelecbottom:ctxelectop, find(ctxelec)')
            warning('cortex electrodes are not consecutive -- check')
        end

        %{
% x, y, z are ccf coordinates. ref: https://alleninstitute.github.io/openscope_databook/visualization/visualize_neuropixels_probes.html
lfpelecx = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('x').data.load();
lfpelecy = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('y').data.load();
lfpelecz = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('z').data.load();
% electrode spacing: from 15 to 50 -- why the large range? expected uniform spacing...
lfpelecspacings = sqrt(diff(lfpelecx).^2+diff(lfpelecy).^2+diff(lfpelecz).^2);

% horizontal position and probe id is the same across all electrodes in the probe, 
% vertical position is spaced 40 apart
lfpprobe_probe_id = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('probe_id').data.load();
lfpprobe_horizontal_position = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('probe_horizontal_position').data.load();
lfpprobe_vertical_position = nwblfp.general_extracellular_ephys_electrodes.vectordata.get('probe_vertical_position').data.load();
        %}



        %% resample from 1250Hz to 1000Hz (same as spiking data)
        % pchip is the best, but results in occasional nan values
        load([pathpp, 'postprocessed_probe' probes{iprobe} '.mat'], 'Tres', 'vis', 'psthtli')

        Nelec = size(lfpdata,1);
        lfptimeresamp = lfptimestamps(1):Tres:lfptimestamps(end);
        lfpresamp = zeros(Nelec, length(lfptimeresamp));
        lfprsnear = zeros(Nelec, length(lfptimeresamp));
        % tic
        for ielec = 1:Nelec
            lfpresamp(ielec,:) = interp1(lfptimestamps, lfpdata(ielec,:), lfptimeresamp,'pchip');
            lfprsnear(ielec,:) = interp1(lfptimestamps, lfpdata(ielec,:), lfptimeresamp,'nearest');
        end
        % toc % ~30 seconds per probe
        lfpresamp(isnan(lfpresamp)) = lfprsnear(isnan(lfpresamp));

        if nnz(isnan(lfpresamp))>0
            error('nan values in lfpresamp: this will lead to faulty TFR')
        end

        save(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{iprobe}), ...
            'lfpelecid', 'lfpelecvec', 'lfptimeresamp', 'lfpresamp', '-v7.3')

        %{
% take a snippet of lfp to determine which interpolation/resampling\ approach to use
tstart = lfptimeresamp(1)+range(lfptimeresamp)*rand(1); tend = tstart+.05; % in seconds
tinds = find(lfptimestamps>=tstart & lfptimestamps<tend);
trsinds = find(lfptimeresamp>=lfptimestamps(tinds(1)) & lfptimeresamp<=lfptimestamps(tinds(end)));

        %{
yticks = 5*10^-4*(0:size(lfpdata,1)-1)';
figure;
plot(lfptimestamps(tinds),yticks+lfpdata(:,tinds));
set(gca, 'YTick', yticks, 'YTickLabel', lfpelecvec.location, 'YDir', 'normal')
        %}


snipelec = round(median(find(contains(lfpelecvec.location, '2/3'))));
lfpsnip = lfpdata(snipelec,tinds);
lfpsnipn = interp1(lfptimestamps(tinds), lfpsnip, lfptimeresamp(trsinds),'nearest');
lfpsnipp = interp1(lfptimestamps(tinds), lfpsnip, lfptimeresamp(trsinds),'pchip');
lfpsnipc = interp1(lfptimestamps(tinds), lfpsnip, lfptimeresamp(trsinds),'cubic');
lfpsnips = interp1(lfptimestamps(tinds), lfpsnip, lfptimeresamp(trsinds),'spline');

Fsorig = 1/median(diff(lfptimestaps));
lfpsnipr = resample(lfpsnip,1000,round(Fsorig));

% very similar, but upon visual inspection, pchip seems to be the best.
% cubic is second best.
figure
for ii = 1:5
    switch ii
        case 1
            templfprs = lfpsnipn;
        case 2
            templfprs = lfpsnipp;%
        case 3
            templfprs = lfpsnipc;
        case 4
            templfprs = lfpsnips;
        case 5
            templfprs = lfpsnipr;
    end
    subplot(2,3,ii)
    hold all
    plot(lfptimestamps(tinds), lfpsnip, 'k-')
    plot(lfptimeresamp(trsinds), templfprs(1:length(trsinds)), 'r.')
    xlim([lfptimestamps(tinds(1)) lfptimestamps(tinds(end))])
end
        %}

        %% align to vis stim
        visblocks = fieldnames(vis);

        lfpvispsth = struct();
        %tic
        for b = 1:numel(visblocks)
            if contains(visblocks{b}, 'spontaneous')
                lfpvispsth.(visblocks{b}) = cell(numel(vis.(visblocks{b}).start_time),1);
                for itrial = 1:numel(vis.(visblocks{b}).start_time)
                    trsinds = find( lfptimeresamp>=vis.(visblocks{b}).start_time(itrial) ...
                        & lfptimeresamp<=vis.(visblocks{b}).stop_time(itrial) );
                    lfpvispsth.(visblocks{b}){itrial} = lfpresamp(:,trsinds)';
                end
                continue
            end

            % actual start time could be up to 1ms after 0 index
            %psthtrialinds = floor(vis.(visblocks{b}).trialstart'/Tres)+1 + psthtli;
            Ntrials = numel(vis.(visblocks{b}).trialstart);
            lfppsthtrialinds = zeros(length(psthtli), Ntrials);
            for itrial = 1:Ntrials
                t0rsind = find(lfptimeresamp<=vis.(visblocks{b}).trialstart(itrial),1,'last');
                lfppsthtrialinds(:,itrial) = t0rsind + psthtli;
            end

            lfpvispsth.(visblocks{b}) = zeros(length(psthtli), vis.(visblocks{b}).numtrials, Nelec);
            for ii = 1:Nelec
                templfp = lfpresamp(ii,:)';
                lfpvispsth.(visblocks{b})(:,:,ii) = templfp(lfppsthtrialinds);
            end
            clear templfp lfppsthtrialinds
        end
        %toc % takes 1 min for 87 electrodes

        %{
% sanity check
whichblock = 'ICwcfg1_presentations';
blocktrialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
tt2p = [0 111 511];
%yticks = 5*10^-4*(0:size(lfpdata,1)-1)';
ytgain = 1*10^-4;
yticks = ytgain*(0:size(lfpdata,1)-1)';
yl = ytgain*[-2 size(lfpdata,1)+1];
figure
for itt = 1:numel(tt2p)
    trialsoi = blocktrialorder==tt2p(itt);
    subplot(1,3,itt)
plot(psthtli, yticks+squeeze(mean(lfpvispsth.(whichblock)(:,trialsoi,:),2))' );
xlim([-200 600])
ylim(yl)
set(gca, 'XGrid', 'on', 'YTick', yticks, 'YTickLabel', lfpelecvec.location, 'YDir', 'normal')
title(tt2p(itt))
end
        %}

        %% align to opto stim: note, there are a lot of NAN values in the lfp block
        load([pathpp, 'psth_opto_probe' probes{iprobe} '.mat'], 'opto', 'optopsthtli')

        Ntrials = numel(opto.optostarttime);
        lfpoptopsthtrialinds = zeros(length(optopsthtli), Ntrials);
        for itrial = 1:Ntrials
            t0rsind = find(lfptimeresamp<=opto.optostarttime(itrial),1,'last');
            lfpoptopsthtrialinds(:,itrial) = t0rsind + optopsthtli;
        end

        lfpoptopsth = zeros(length(optopsthtli), Ntrials, Nelec);
        % tic
        for ii = 1:Nelec
            templfp = lfpresamp(ii,:)';
            lfpoptopsth(:,:,ii) = templfp(lfpoptopsthtrialinds);
        end
        clear templfp psthtrialinds
        % toc % 4 seconds

        % nanprct=100*mean(isnan(lfpoptopsth(optopsthtli>=0&optopsthtli<1000, :,:)), 'all');
        % fprintf('%.1f%% time points are nan during opto stim\n', nanprct)
        if nnz(isnan(lfpoptopsth))>0
            error('nan values in lfpoptopsth: this will lead to faulty TFR')
        end

        %{
% sanity check
tt2p = [3 7 11];
%yticks = 5*10^-4*(0:size(lfpdata,1)-1)';
ytgain = 1*10^-4;
yticks = ytgain*(0:size(lfpdata,1)-1)';
yl = ytgain*[-2 size(lfpdata,1)+1];
figure
for itt = 1:numel(tt2p)
    trialsoi = opto.optotrials==tt2p(itt);
    subplot(1,3,itt)
plot(optopsthtli, yticks+squeeze(mean(lfpoptopsth(:,trialsoi,:),2))' );
xlim([-200 600])
ylim(yl)
set(gca, 'XGrid', 'on', 'YTick', yticks, 'YTickLabel', lfpelecvec.location, 'YDir', 'normal')
title( opto.stimlist(tt2p(itt)) )
end
        %}

        save(sprintf('%sLFP_psth_probe%s.mat', pathpp, probes{iprobe}), ...
            'Tres', 'vis', 'psthtli', 'lfpvispsth', ...
            'opto', 'optopsthtli', 'lfpoptopsth', '-v7.3')


        %% CSD: negative values mean sink (influx, depol), positive values mean source (outflux)
        % lfpelecspacing = 0.04; % 40micrometers, i.e., 0.04mm
        % csdelectinds = 2:Nelec-1;
        % csdresamp = -( lfpresamp(csdelectinds+1,:)-2*lfpresamp(csdelectinds,:)+lfpresamp(csdelectinds-1,:) )/(lfpelecspacing.^2);

        % smooth with a gaussian kernel first
        kerwinhalf = 5; kersigma = 2;
        kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
        kergauss = (kergauss/sum(kergauss));
        lfpconv = convn(lfpresamp, kergauss, 'same');

        lfpelecspacing = 0.04; % 40micrometers, i.e., 0.04mm
        csdelectinds = 2:Nelec-1;
        csdconv = -( lfpconv(csdelectinds+1,:)-2*lfpconv(csdelectinds,:)+lfpconv(csdelectinds-1,:) )/(lfpelecspacing.^2);

        csdvispsth = struct();
        % tic
        for b = 1:numel(visblocks)
            if contains(visblocks{b}, 'spontaneous')
                csdvispsth.(visblocks{b}) = cell(numel(vis.(visblocks{b}).start_time),1);
                for itrial = 1:numel(vis.(visblocks{b}).start_time)
                    t0rsind = find( lfptimeresamp>=vis.(visblocks{b}).start_time(itrial) ...
                        & lfptimeresamp<=vis.(visblocks{b}).stop_time(itrial) );
                    csdvispsth.(visblocks{b}){itrial} = csdconv(:,t0rsind)';
                end
                continue
            end

            % actual start time could be up to 1ms after 0 index
            %psthtrialinds = floor(vis.(visblocks{b}).trialstart'/Tres)+1 + psthtli;
            Ntrials = numel(vis.(visblocks{b}).trialstart);
            lfppsthtrialinds = zeros(length(psthtli), Ntrials);
            for itrial = 1:Ntrials
                t0rsind = find(lfptimeresamp<=vis.(visblocks{b}).trialstart(itrial),1,'last');
                lfppsthtrialinds(:,itrial) = t0rsind + psthtli;
            end

            csdvispsth.(visblocks{b}) = NaN(length(psthtli), vis.(visblocks{b}).numtrials, Nelec-2);
            for ii = 1:Nelec-2
                tempcsd = csdconv(ii,:)';
                csdvispsth.(visblocks{b})(:,:,ii) = tempcsd(lfppsthtrialinds);
            end
            clear tempcsd lfppsthtrialinds
        end
        % toc % takes 1 min for 87 electrodes

        Ntrials = numel(opto.optostarttime);
        lfpoptopsthtrialinds = zeros(length(optopsthtli), Ntrials);
        for itrial = 1:Ntrials
            t0rsind = find(lfptimeresamp<=opto.optostarttime(itrial),1,'last');
            lfpoptopsthtrialinds(:,itrial) = t0rsind + optopsthtli;
        end

        csdoptopsth = NaN(length(optopsthtli), Ntrials, Nelec-2);
        % tic
        for ii = 1:Nelec-2
            tempcsd = csdconv(ii,:)';
            csdoptopsth(:,:,ii) = tempcsd(lfpoptopsthtrialinds);
        end
        clear tempcsd psthtrialinds
        % toc

        %{
% compare with averaging first then calculating CSD
whichblock = 'ICwcfg1_presentations';
blocktrialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
tt2p = [0 111 511];
xl = [0 200];
yl = [ctxelecbottom ctxelectop]+1;
figure
for itt = 1:numel(tt2p)
    trialsoi = blocktrialorder==tt2p(itt);
    subplot(2,3,itt)
    imagesc(psthtli, csdelectinds, squeeze(mean( csdvispsth.(whichblock)(:,trialsoi,:),2))' );
    xlim(xl)
    ylim(yl)
    caxis([-0.025 0.025])
    set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
    title(tt2p(itt))
    colorbar

    lfpavg = squeeze(mean( convn(lfpvispsth.(whichblock)(:,trialsoi,:),kergauss', 'same') ,2))';
    csdavg = -( lfpavg(csdelectinds+1,:)-2*lfpavg(csdelectinds,:)+lfpavg(csdelectinds-1,:) )/(lfpelecspacing.^2);
    subplot(2,3,itt+3)
    imagesc(psthtli, csdelectinds, csdavg );
    xlim(xl)
    ylim(yl)
    caxis([-0.025 0.025])
    set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
    title(tt2p(itt))
    colorbar
end
colormap jet

% compare csdresamp vs csdconv
t0ind = prctile(trsinds,10); % one of the spontaneous blocks
tdur = 1000;
xl = [t0ind t0ind+tdur] - trsinds(1);
yl = [ctxelecbottom ctxelectop]+1;
cl = [-0.1 0.1];
figure
subplot(1,2,1)
imagesc(csdresamp(:,trsinds))
xlim(xl)
set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
colorbar
%cl=caxis;
caxis(cl)
ylim(yl)
subplot(1,2,2)
imagesc(csdconv(:,trsinds))
set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
colorbar
caxis(cl)
xlim(xl)
ylim(yl)
colormap jet
        %}

        save(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{iprobe}), ...
            'Tres', 'kergauss', 'lfpelecspacing', 'csdelectinds', ...
            'vis', 'psthtli', 'csdvispsth', ...
            'opto', 'optopsthtli', 'csdoptopsth', '-v7.3')

        %% TFR : for now, just choose one electrode in layer 2/3
        %elecL23 = round(median( find(contains(lfpelecvec.location, '2/3')) ));
        electrodesL23 = find(contains(lfpelecvec.location(csdelectinds), '2/3'));

        % CHOOSE THE ELECTRODE WITH THE STRONGEST SINK IN FEEDBACK PERIOD IN L2/3
        whichblock = 'ICwcfg1_presentations';
        blocktrialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
        trialsoi = ismember(blocktrialorder, [1105, 1100]);
        tloi = psthtli>=100 & psthtli<150;
        [mv,mi]=min( squeeze(mean( csdvispsth.(whichblock)(tloi,trialsoi,electrodesL23),[1,2])) );
        elecL23 = csdelectinds( electrodesL23(mi) ) ;

        %{
yl = [ctxelecbottom ctxelectop];
figure
tempcsd = squeeze(mean( csdvispsth.(whichblock)(:,trialsoi,:),2))';
imagesc(psthtli, csdelectinds, tempcsd)
hold on
templfp = squeeze(mean( lfpvispsth.(whichblock)(:,trialsoi,csdelectinds),2))';
lfpfac = 1.5/prctile(range(templfp,2),50);
plot(psthtli, csdelectinds'+lfpfac*templfp, 'k-')
scatter(125, elecL23,100,'w*', 'linewidth', 1)
set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
colorbar
%cl=caxis;
caxis([-0.02 0.02])
xlim([0 200])
ylim(yl)
colormap jet
        %}

        % 30 seconds and 7 GB for each electrode
        S = lfpresamp(elecL23,:)';
        % get rid of NaN values in S
        [TFR_L23,tVec,fVec] = spectralevents_ts2tfr(S,1:100,1/Tres,5);


        tfrvispsth = struct();
        %tic
        for b = 1:numel(visblocks)
            if contains(visblocks{b}, 'spontaneous')
                tfrvispsth.(visblocks{b}) = cell(numel(vis.(visblocks{b}).start_time),1);
                for itrial = 1:numel(vis.(visblocks{b}).start_time)
                    trsinds = find( lfptimeresamp>=vis.(visblocks{b}).start_time(itrial) ...
                        & lfptimeresamp<=vis.(visblocks{b}).stop_time(itrial) );
                    tfrvispsth.(visblocks{b}){itrial} = TFR_L23(:,trsinds)';
                end
                continue
            end

            % actual start time could be up to 1ms after 0 index
            %psthtrialinds = floor(vis.(visblocks{b}).trialstart'/Tres)+1 + psthtli;
            Ntrials = numel(vis.(visblocks{b}).trialstart);
            lfppsthtrialinds = zeros(length(psthtli), Ntrials);
            for itrial = 1:Ntrials
                t0rsind = find(lfptimeresamp<=vis.(visblocks{b}).trialstart(itrial),1,'last');
                lfppsthtrialinds(:,itrial) = t0rsind + psthtli;
            end

            tfrvispsth.(visblocks{b}) = zeros(length(psthtli), vis.(visblocks{b}).numtrials, length(fVec));
            for ii = 1:length(fVec)
                temptfr = TFR_L23(ii,:)';
                tfrvispsth.(visblocks{b})(:,:,ii) = temptfr(lfppsthtrialinds);
            end
            clear temptfr lfppsthtrialinds
        end

        Ntrials = numel(opto.optostarttime);
        lfpoptopsthtrialinds = zeros(length(optopsthtli), Ntrials);
        for itrial = 1:Ntrials
            t0rsind = find(lfptimeresamp<=opto.optostarttime(itrial),1,'last');
            lfpoptopsthtrialinds(:,itrial) = t0rsind + optopsthtli;
        end
        tfroptopsth = zeros(length(optopsthtli), Ntrials, length(fVec));
        % tic
        for ii = 1:length(fVec)
            temptfr = TFR_L23(ii,:)';
            tfroptopsth(:,:,ii) = temptfr(lfpoptopsthtrialinds);
        end
        clear temptfr psthtrialinds

        %{
tstart = 1800; tend = tstart+10; % in seconds
trsinds = find(lfptimeresamp>=tstart & lfptimeresamp<=tend );
[TFRsnip,~,~] = spectralevents_ts2tfr(lfpresamp(elecL23,trsinds)',fVec,1/Tres,5);
figure;
subplot(1,2,1)
imagesc(TFR_L23(:,trsinds))
colorbar
subplot(1,2,2)
imagesc(TFRsnip)
colorbar
colormap jet
corr(reshape(TFR_L23(:,trsinds),[],1), reshape(TFRsnip,[],1)) % 0.9934 correlation

% sanity check
whichblock = 'ICwcfg1_presentations';
blocktrialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
tt2p = [0 111 511];
figure
for itt = 1:numel(tt2p)
    trialsoi = blocktrialorder==tt2p(itt);
    subplot(1,3,itt)
imagesc(psthtli, fVec, squeeze(mean( log10(tfrvispsth.(whichblock)(:,trialsoi,:)),2))' );
xlim([-200 600])
clim([-10 -9])
set(gca, 'XGrid', 'on')
title(tt2p(itt))
colorbar
end


whichblock = 'ICwcfg1_presentations';
blocktrialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
tt2p = [0 111 511];
tloi = psthtli>=0 & psthtli<400;
figure
hold all
for itt = 1:numel(tt2p)
    trialsoi = blocktrialorder==tt2p(itt);
    plot(fVec, squeeze(mean( log10(tfrvispsth.(whichblock)(tloi,trialsoi,:)),[1, 2])) )
end
legend(tt2p)
        %}

        save(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{iprobe}), ...
            'Tres', 'elecL23', 'fVec', 'vis', 'psthtli', 'tfrvispsth', ...
            'opto', 'optopsthtli', 'tfroptopsth', '-v7.3')
toc(probeclk)
    end
end
% check here for camstim orientation: clockwise, i.e., 30deg is 1o'clock
% http://observatory.brain-map.org/visualcoding/stimulus/static_gratings

clear all

addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130')
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

%%
for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions datadir
    sesclk = tic;

    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));

    % take filename  with shortest length or filename that does not contain probe
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    % nwbspikefile = string(nwbspikefile);
    disp(nwbspikefile)
    nwb = nwbRead(nwbspikefile); %, 'ignorecache');

    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    if ~exist(pathpp, 'dir')
        mkdir(pathpp)
    end

    %% where are the units from? probe position and within-probe electrode
    % position? brain area labels?
    % There's a 'peak channel' field associated with each unit.
    % The formula for the id is 1000 *probe_val + the local electrode id.
    % This should be the same as the ids that are present in the id field of the
    % electrode section of the NWB (i.e. the 2304th value is 5383 because it is
    % the 383rd value of the 5th probe.)

    % nwb.general_extracellular_ephys
    electrode_probeid = nwb.general_extracellular_ephys_electrodes.vectordata.get('probe_id').data.load();
    electrode_localid = nwb.general_extracellular_ephys_electrodes.vectordata.get('local_index').data.load();
    electrode_id = 1000*electrode_probeid + electrode_localid;
    electrode_location = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();

    disp(unique(electrode_location)')
    save([pathpp 'info_electrodes.mat'], 'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')

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

    Tres = 0.001; % 1ms
    stlen = ceil((max(unit_times_data)+1)/Tres); % add 1s buffer/padding after the last spike timing

    disp([stlen, Nneurons])

    save([pathpp 'info_units.mat'], 'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',

    %% extract vis stim info
    % in IC blocks frame contains trialorder info
    % in RFCI and sizeCI blocks, Mask, orientation, x_position and  y_position

    % note from Ahad via Slack DM 220520
    % image_path_list variable that is currently not being entered into the nwb file.

    visblocks = nwb.intervals.keys;
    vis = struct();
    for b = 1:numel(visblocks)
        disp(visblocks{b})
        vis.(visblocks{b}).start_time = nwb.intervals.get(visblocks{b}).start_time.data.load();
        vis.(visblocks{b}).stop_time = nwb.intervals.get(visblocks{b}).stop_time.data.load();
        % vis.(visblocks{b}).tags = nwb.intervals.get(visblocks{b}).tags.data.load();
        % vis.(visblocks{b}).tags_index = nwb.intervals.get(visblocks{b}).tags_index.data.load();
        % vis.(visblocks{b}).timeseries = nwb.intervals.get(visblocks{b}).timeseries.data.load();
        % vis.(visblocks{b}).timeseries_index = nwb.intervals.get(visblocks{b}).timeseries_index.data.load();
        % vis.(visblocks{b}).id = nwb.intervals.get(visblocks{b}).id.data.load();
        viskeys =  nwb.intervals.get(visblocks{b}).vectordata.keys;
        for k = 1:numel(viskeys)
            vis.(visblocks{b}).(viskeys{k}) = nwb.intervals.get(visblocks{b}).vectordata.get(viskeys{k}).data.load();
        end

        % IC blocks
        if ismember('frame', viskeys)
            % expect 61, 31, 31, 31 for frame_firsttrial
            frame_firsttrial = find(vis.(visblocks{b}).frame~=0, 1, 'first');
            if mod(frame_firsttrial, 10) ~=1
                warning('first trial was blank')
                frame_firsttrial = 10*floor(frame_firsttrial/10)+1;
            end
            frame_lasttrial = find(vis.(visblocks{b}).frame~=0, 1, 'last');
            if mod(frame_lasttrial, 10) ~=9
                warning('last trial was blank')
                frame_lasttrial = 10*floor(frame_lasttrial/10)+9;
            end

            trialframeinds = frame_firsttrial:2:frame_lasttrial;
            vis.(visblocks{b}).trialframeinds = trialframeinds;
            vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time(trialframeinds);
            vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time(trialframeinds);
            vis.(visblocks{b}).numtrials = length(trialframeinds);
            vis.(visblocks{b}).trialorder = vis.(visblocks{b}).frame(trialframeinds);
            if contains(visblocks{b}, 'cfg1')
                vis.(visblocks{b}).trialtypedescription = {'Blank', 'X', 'TC1', 'IC1', 'LC1', 'TC2', 'LC2', 'IC2', ...
                    'IRE1', 'IRE2', 'TRE1', 'TRE2', 'XRE1', 'XRE2', ...
                    'InBR', 'InBL', 'InTL', 'InTR', 'OutBR', 'OutBL', 'OutTL', 'OutTL'};
            elseif contains(visblocks{b}, 'cfg0')
                vis.(visblocks{b}).trialtypedescription = {'Blank', 'X', 'TC1', 'IC1', 'LC1', 'TC2', 'LC2', 'IC2', ...
                    'IRE1', 'IRE2', 'TRE1', 'TRE2', 'XRE1', 'XRE2', ...
                    'InR', 'InB', 'InL', 'InT', 'OutR', 'OutB', 'OutL', 'OutT'};
            else
                error('unrecognized configuration')
            end
            vis.(visblocks{b}).ICtrialtypes = [0 101 105 106 107 109 110 111 ...
                506 511 1105 1109 1201 1299 ...
                1301 1302 1303 1304 1305 1306 1307 1308];

            disp([frame_firsttrial, frame_lasttrial vis.(visblocks{b}).numtrials])
            if contains(visblocks{b}, 'ICwcfg1') || (contains(nwbsessions{ises}, 'Placeholder') && b==1)
                expectedNtrials = 5300; % 12*400+10*50
            elseif contains(nwbsessions{ises}, 'Placeholder')
                expectedNtrials = 550; % 22*25
            else
                expectedNtrials = 22*30;
            end
            if ~( vis.(visblocks{b}).numtrials==expectedNtrials )
                error('check numtrials')
            end
        end

        % RFCIblocks
        % (+right,+up))
        % rfpos = [(0,0), (0,-203.3786), (203.3786/2**0.5,-203.3786/2**0.5), (203.3786,0), \
        %             (203.3786/2**0.5,203.3786/2**0.5), (0,203.3786), (-203.3786/2**0.5,203.3786/2**0.5), \
        %             (-203.3786,0), (-203.3786/2**0.5,-203.3786/2**0.5)]
        % '10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
        if ismember('y_position', viskeys)
            vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time;
            vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time;

            yx_position = [vis.(visblocks{b}).y_position vis.(visblocks{b}).x_position];
            vis.(visblocks{b}).sizepix = 203.3786;
            vis.(visblocks{b}).MaskDiaVisDeg = 16;
            % note, made the order of RFcentersrel match that in matlab (start
            % down, then move counterclockwise)
            % Y-DIRECTION FLIPPED FROM PSYCHTOOLBOX (in psychtoolbox, + is
            % down; in openscope dataset, -1 is down)
            vis.(visblocks{b}).RFcentersrel = [-1 0
                -1/sqrt(2) 1/sqrt(2)
                0 1
                1/sqrt(2) 1/sqrt(2)
                1 0
                1/sqrt(2) -1/sqrt(2)
                0 -1
                -1/sqrt(2) -1/sqrt(2)];
            vis.(visblocks{b}).RFcenters = vis.(visblocks{b}).sizepix * vis.(visblocks{b}).RFcentersrel;
            vis.(visblocks{b}).RFcentersVisDeg = vis.(visblocks{b}).MaskDiaVisDeg * vis.(visblocks{b}).RFcentersrel;
            if ~all(ismember(vis.(visblocks{b}).RFcenters, unique(yx_position, 'rows'), 'rows'))
                error('check RFcenters')
            end

            vis.(visblocks{b}).directions = unique(vis.(visblocks{b}).orientation);
            vis.(visblocks{b}).MaskList = unique(vis.(visblocks{b}).Mask);
            disp(vis.(visblocks{b}).MaskList)

            % vertical is zero, then rotates clockwise (45 is 1.5o'clock)
            vis.(visblocks{b}).numtrials = length(vis.(visblocks{b}).orientation);
            vis.(visblocks{b}).trialorder = zeros(vis.(visblocks{b}).numtrials, 1);
            for typi = 1:numel(vis.(visblocks{b}).directions)
                trialsoi = vis.(visblocks{b}).orientation==vis.(visblocks{b}).directions(typi);
                vis.(visblocks{b}).trialorder(trialsoi) = typi + vis.(visblocks{b}).trialorder(trialsoi);
            end
            for typi = 1:size(vis.(visblocks{b}).RFcenters,1)
                trialsoi = ismember(yx_position, vis.(visblocks{b}).RFcenters(typi,:), 'rows');
                vis.(visblocks{b}).trialorder(trialsoi) = 10*typi + vis.(visblocks{b}).trialorder(trialsoi);
            end
            for typi = 1:numel(vis.(visblocks{b}).MaskList)
                trialsoi = strcmp(vis.(visblocks{b}).Mask, vis.(visblocks{b}).MaskList(typi));
                tempss = strsplit(vis.(visblocks{b}).MaskList{typi}, '\');
                tempss = strsplit(tempss{end}, '.tif');
                maskno = str2num(tempss{1});
                vis.(visblocks{b}).trialorder(trialsoi) = maskno + vis.(visblocks{b}).trialorder(trialsoi);
            end
            vis.(visblocks{b}).trialtypedescription = ['10000s: classic 0 vs inverse 1,', ...
                ' 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'];
        end

        if contains(visblocks{b}, 'sizeCI')
            vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time;
            vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time;

            vis.(visblocks{b}).directions = unique(vis.(visblocks{b}).orientation);
            vis.(visblocks{b}).MaskList = unique(vis.(visblocks{b}).Mask);
            disp(vis.(visblocks{b}).MaskList)

            vis.(visblocks{b}).numtrials = length(vis.(visblocks{b}).orientation);
            vis.(visblocks{b}).trialorder = zeros(vis.(visblocks{b}).numtrials, 1);
            for typi = 1:numel(vis.(visblocks{b}).directions)
                trialsoi = vis.(visblocks{b}).orientation==vis.(visblocks{b}).directions(typi);
                vis.(visblocks{b}).trialorder(trialsoi) = typi + vis.(visblocks{b}).trialorder(trialsoi);
            end
            for typi = 1:numel(vis.(visblocks{b}).MaskList)
                trialsoi = strcmp(vis.(visblocks{b}).Mask, vis.(visblocks{b}).MaskList(typi));
                tempss = strsplit(vis.(visblocks{b}).MaskList{typi}, '\');
                tempss = strsplit(tempss{end}, '.tif');
                maskno = str2num(tempss{1});
                vis.(visblocks{b}).trialorder(trialsoi) = maskno + vis.(visblocks{b}).trialorder(trialsoi);
            end
            vis.(visblocks{b}).MaskDiaVisDeg = [0, 4, 8, 16, 32, 64];
            vis.(visblocks{b}).trialtypedescription = ['10000s: classic 0 vs inverse 1,', ...
                ' 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'];
        end

        % disp(unique(vis.(visblocks{b}).trialorder)')

    end

    % vistimestamps = nwb.processing.get('stimulus').nwbdatainterface.get('timestamps').timestamps.load();
    % visdata = nwb.processing.get('stimulus').nwbdatainterface.get('timestamps').data.load();
    % isequal(visdata, vistimestamps) % true
    % isequal(vistimestamps(1:ntrials), vis.(visblocks{b}).start_time) % true


    %% spontaneous block
    %{
visblockorder = [4 3 2 1 5 6];
blockstarts = zeros(size(visblocks));
blockends = zeros(size(visblocks));
for b = 1:numel(visblocks)
    blockstarts(b) = vis.(visblocks{b}).start_time(1);
    blockends(b) = vis.(visblocks{b}).stop_time(end);
end
% blockstarts(visblockorder(2:end))-blockends(visblockorder(1:end-1))
if ~isequal(vis.spontaneous_presentations.stop_time, blockstarts(visblockorder)')
    error('check blockstarts/visblockorder')
end
if ~isequal(vis.spontaneous_presentations.start_time(2:end), blockends(visblockorder(1:end-1))')
    error('check blockends/visblockorder')
end
    %}
    % start_time and stop_time are in seconds
    [mv,mi] = max(vis.spontaneous_presentations.stop_time - vis.spontaneous_presentations.start_time);
    fprintf('longest spontaneous block is %.0fs\n', mv)
    sponTstartind = floor(vis.spontaneous_presentations.start_time(mi)'/Tres)+1;
    sponTendind = floor(vis.spontaneous_presentations.stop_time(mi)'/Tres);
    % figure; hold all
    % plot(sponFR, meanFR, 'k.')
    % xl = xlim; yl = ylim; al = [min([xl yl]) max([xl yl])];
    % plot(al, al, 'r-')
    % axis([al al])

    %% psthall and Rall
    % A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
    % visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
    % visind = [6 5 1 2 4 3];

    probes = {'A', 'B', 'C', 'D', 'E', 'F'};

    psthtli = (-500:1000)';
    neucnt = 0;
    for iprobe = 1:numel(probes)
        neuoind = find(floor(unit_peakch/1000)==iprobe-1);
        neucnt = neucnt + numel(neuoind);
        fprintf('Probe %s: %d\n', probes{iprobe}, numel(neuoind) )

        if numel(neuoind)==0
            continue
        end

        probespiketrain = false(stlen, numel(neuoind));
        ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
        for ii = 1:numel(neuoind)
            ci = neuoind(ii);
            probespiketrain(floor(spiketimes{ci}/Tres)+1, ii) = true;
        end


        ppfn = sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe});
        if exist(ppfn, 'file')
            load(ppfn)
        else
            psth = struct();
            tic
            for b = 1:numel(visblocks)
                if contains(visblocks{b}, 'spontaneous')
                    continue
                end
                psthtrialinds = floor(vis.(visblocks{b}).trialstart'/Tres)+1 + psthtli;
                psth.(visblocks{b}) = false(length(psthtli), vis.(visblocks{b}).numtrials, numel(neuoind));
                for ii = 1:numel(neuoind)
                    tempST = probespiketrain(:,ii);
                    psth.(visblocks{b})(:,:,ii) = tempST(psthtrialinds);
                end
                clear tempST psthtrialinds
            end
            toc % takes 6 min for 2500 units

            save(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}), ...
                'neuoind', 'vis', 'Tres', 'psthtli', 'psth', '-v7.3')
        end

        % figure
        % for b = 1:numel(visblocks)
        %     subplot(2,3,b)
        %     plot(psthtli, 1000*squeeze(mean(psthall.(visblocks{b}), [2,3])))
        %     xlabel('Time (ms)')
        %     ylabel('Rate (Hz)')
        % end

        sponFRvec = mean(probespiketrain(sponTstartind:sponTendind, :),1)/Tres;
        meanFRvec = mean(probespiketrain,1)/Tres;

        % %%
        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        ICblocks = 1:4;
        ICsig = struct();
        for b = ICblocks
            disp(visblocks{b})
            tloi = psthtli>0 & psthtli<=400;
            tempR = squeeze(1000*mean(psth.(visblocks{b})(tloi,:,:), 1))';
            temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
            ICsig.(visblocks{b}) = analyzeStaticICtxi(tempR, temptrialorder);
        end
        % 'visblocknames', 'trialtypes', 'ICblocks',

        % RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning grating
        % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
        % durstim = vis.RFCI_presentations.stop_time-vis.RFCI_presentations.start_time;
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
        % durstim = vis.RFCI_presentations.start_time(2:end)-vis.RFCI_presentations.stop_time(1:end-1);
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])

        %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
        % 10000's: which type (classic 0 vs inverse 1), 1000's which ctrsizes,
        % 10-100's: which RFcenter, 1's: which direction
        tloi = psthtli>0 & psthtli<=1000;
        tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,1:4:end,:), 1))';
        temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
        RFCI = analyzeRFCI(tempR, temptrialorder, sponFRvec);

        tloi = psthtli>0 & psthtli<=250;
        tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,:,:), 1))';
        temptrialorder = vis.RFCI_presentations.trialorder;
        RFCIspin = analyzeRFCIspin(tempR, temptrialorder, sponFRvec);

        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        %szvec = [0, 4, 8, 16, 32, 64];
        % durstim = vis.sizeCI_presentations.stop_time-vis.sizeCI_presentations.start_time;
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
        % durstim = vis.sizeCI_presentations.start_time(2:end)-vis.sizeCI_presentations.stop_time(1:end-1);
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
        tloi = psthtli>0 & psthtli<=250;
        tempR = squeeze(1000*mean(psth.sizeCI_presentations(tloi,:,:), 1))';
        temptrialorder = vis.sizeCI_presentations.trialorder;
        [sizeCI, oriparams, ori4params] = analyzesizeCI(tempR, temptrialorder);

        save(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}), ...
            'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'RFCIspin', ...
            'sizeCI', 'oriparams', 'ori4params', '-v7.3')

    end
    %%
    if neucnt~=Nneurons
        error('not all neurons were accounted for')
    end


    toc(sesclk)
end

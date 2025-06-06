%% merge vis variables e.g., ICsigall
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

%%
for ises = 1:Nsessions
    clearvars -except ises Nsessions nwbsessions datadir
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    tic
    
    gazedistthresh = 20;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    load([pathpp 'postprocessed_probeC.mat'])%, 'psthtli', 'vis')
    load([pathpp 'visresponses_probeC.mat'])
    validfixgaze = exist([pathpp 'trackmouseeye.mat'], 'file');
    if validfixgaze
        load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh))
    end

    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    neuallloc = electrode_location(revmapelecid(unit_peakch+1));

    probes = {'A', 'B', 'C', 'D', 'E', 'F'};
    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
    Nneurons = numel(neuoind);
    
    Nneuronsall = length(unit_peakch);
    meanFRall = NaN(1, Nneuronsall);
    sponFRall = NaN(1, Nneuronsall);
    Rall = struct();
    for b = 1:numel(visblocks)
        Rall.(visblocks{b}) = NaN(vis.(visblocks{b}).numtrials, Nneuronsall);
    end
    
    ICblocks = {'ICwcfg1_presentations','ICwcfg0_presentations','ICkcfg1_presentations','ICkcfg0_presentations'};
    % 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
    allfields = fieldnames(ICsig.ICwcfg1_presentations);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(ICsig.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
    end
    ICsigfields = allfields(validfields);
    % {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
    %     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
    %     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
    %     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
    %     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
    %     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
    %     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
    %     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
    ICsigall = struct();
    for b = 1:numel(ICblocks)
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(ICsig.(ICblocks{b}).(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                ICsigall.(ICblocks{b}).(allfields{f}) = NaN(tempsizeall);
            else
                ICsigall.(ICblocks{b}).(allfields{f}) = ICsig.(ICblocks{b}).(allfields{f});
            end
        end
    end
    
    allfields = fieldnames(RFCI);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(RFCI.(allfields{f}),1)==Nneurons;
    end
    RFCIfields = allfields(validfields);
    % {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
    %     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
    RFCIall = struct();
    for f= 1:numel(allfields)
        if validfields(f)
            tempsize = size(RFCI.(allfields{f}));
            tempsizeall = tempsize;
            tempsizeall(tempsize==Nneurons) = Nneuronsall;
            RFCIall.(allfields{f}) = NaN(tempsizeall);
        else
            RFCIall.(allfields{f}) = RFCI.(allfields{f});
        end
    end
    
    allfields = fieldnames(RFCIspin);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(RFCIspin.(allfields{f}),1)==Nneurons;
    end
    RFCIspinfields = allfields(validfields);
    RFCIspinall = struct();
    for f= 1:numel(allfields)
        if validfields(f)
            tempsize = size(RFCIspin.(allfields{f}));
            tempsizeall = tempsize;
            tempsizeall(tempsize==Nneurons) = Nneuronsall;
            RFCIspinall.(allfields{f}) = NaN(tempsizeall);
        else
            RFCIspinall.(allfields{f}) = RFCIspin.(allfields{f});
        end
    end
    
    allfields = fieldnames(sizeCI);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(sizeCI.(allfields{f}),1)==Nneurons;
    end
    sizeCIfields = allfields(validfields);
    % sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
    sizeCIall = struct();
    for f= 1:numel(allfields)
        if validfields(f)
            tempsize = size(sizeCI.(allfields{f}));
            tempsizeall = tempsize;
            tempsizeall(tempsize==Nneurons) = Nneuronsall;
            sizeCIall.(allfields{f}) = NaN(tempsizeall);
        else
            sizeCIall.(allfields{f}) = sizeCI.(allfields{f});
        end
    end
    
    allfields = fieldnames(oriparams);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(oriparams.(allfields{f}),1)==Nneurons;
    end
    oriparamsfields = allfields(validfields);
    % oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
    oriparamsall = struct();
    for f= 1:numel(allfields)
        if validfields(f)
            tempsize = size(oriparams.(allfields{f}));
            tempsizeall = tempsize;
            tempsizeall(tempsize==Nneurons) = Nneuronsall;
            oriparamsall.(allfields{f}) = NaN(tempsizeall);
        else
            oriparamsall.(allfields{f}) = oriparams.(allfields{f});
        end
    end
    
    allfields = fieldnames(ori4params);
    validfields = true(size(allfields));
    for f = 1:numel(allfields)
        validfields(f) = size(ori4params.(allfields{f}),1)==Nneurons;
    end
    ori4paramsfields = allfields(validfields);
    % oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
    ori4paramsall = struct();
    for f= 1:numel(allfields)
        if validfields(f)
            tempsize = size(ori4params.(allfields{f}));
            tempsizeall = tempsize;
            tempsizeall(tempsize==Nneurons) = Nneuronsall;
            ori4paramsall.(allfields{f}) = NaN(tempsizeall);
        else
            ori4paramsall.(allfields{f}) = ori4params.(allfields{f});
        end
    end
    
    if validfixgaze
        ICblocks = {'ICwcfg1_presentations','ICwcfg0_presentations','ICkcfg1_presentations','ICkcfg0_presentations'};
        % 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
        allfields = fieldnames(ICsig_fixedgaze.ICwcfg1_presentations);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
        end
        ICsigall_fixedgaze = struct();
        for b = 1:numel(ICblocks)
            if isempty(fieldnames(ICsig_fixedgaze.(ICblocks{b})))
                ICsigall_fixedgaze.(ICblocks{b}) = struct();
            else
            for f= 1:numel(allfields)
                if validfields(f)
                    tempsize = size(ICsig_fixedgaze.(ICblocks{b}).(allfields{f}));
                    tempsizeall = tempsize;
                    tempsizeall(tempsize==Nneurons) = Nneuronsall;
                    ICsigall_fixedgaze.(ICblocks{b}).(allfields{f}) = NaN(tempsizeall);
                else
                    ICsigall_fixedgaze.(ICblocks{b}).(allfields{f}) = ICsig_fixedgaze.(ICblocks{b}).(allfields{f});
                end
            end
            end
        end
        
        allfields = fieldnames(RFCI_fixedgaze);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(RFCI_fixedgaze.(allfields{f}),1)==Nneurons;
        end
        % RFCIfields = allfields(validfields);
        % {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
        %     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
        RFCIall_fixedgaze = struct();
        if ~isempty(fieldnames(RFCI_fixedgaze))
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(RFCI_fixedgaze.(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                RFCIall_fixedgaze.(allfields{f}) = NaN(tempsizeall);
            else
                RFCIall_fixedgaze.(allfields{f}) = RFCI_fixedgaze.(allfields{f});
            end
        end
        end
        
        allfields = fieldnames(RFCIspin_fixedgaze);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(RFCIspin_fixedgaze.(allfields{f}),1)==Nneurons;
        end
        % RFCIspinfields = allfields(validfields);
        RFCIspinall_fixedgaze = struct();
        if ~isempty(fieldnames(RFCIspin_fixedgaze))
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(RFCIspin_fixedgaze.(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                RFCIspinall_fixedgaze.(allfields{f}) = NaN(tempsizeall);
            else
                RFCIspinall_fixedgaze.(allfields{f}) = RFCIspin_fixedgaze.(allfields{f});
            end
        end
        end
        
        allfields = fieldnames(sizeCI_fixedgaze);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(sizeCI_fixedgaze.(allfields{f}),1)==Nneurons;
        end
        % sizeCIfields = allfields(validfields);
        % sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
        sizeCIall_fixedgaze = struct();
        if ~isempty(fieldnames(sizeCI_fixedgaze))
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(sizeCI_fixedgaze.(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                sizeCIall_fixedgaze.(allfields{f}) = NaN(tempsizeall);
            else
                sizeCIall_fixedgaze.(allfields{f}) = sizeCI_fixedgaze.(allfields{f});
            end
        end
        end
        
        allfields = fieldnames(oriparams_fixedgaze);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(oriparams_fixedgaze.(allfields{f}),1)==Nneurons;
        end
        % oriparamsfields = allfields(validfields);
        % oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
        oriparamsall_fixedgaze = struct();
        if ~isempty(fieldnames(oriparams_fixedgaze))
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(oriparams_fixedgaze.(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                oriparamsall_fixedgaze.(allfields{f}) = NaN(tempsizeall);
            else
                oriparamsall_fixedgaze.(allfields{f}) = oriparams_fixedgaze.(allfields{f});
            end
        end
        end
        
        allfields = fieldnames(ori4params_fixedgaze);
        validfields = true(size(allfields));
        for f = 1:numel(allfields)
            validfields(f) = size(ori4params_fixedgaze.(allfields{f}),1)==Nneurons;
        end
        % ori4paramsfields = allfields(validfields);
        % oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
        ori4paramsall_fixedgaze = struct();
        if ~isempty(fieldnames(ori4params_fixedgaze))
        for f= 1:numel(allfields)
            if validfields(f)
                tempsize = size(ori4params_fixedgaze.(allfields{f}));
                tempsizeall = tempsize;
                tempsizeall(tempsize==Nneurons) = Nneuronsall;
                ori4paramsall_fixedgaze.(allfields{f}) = NaN(tempsizeall);
            else
                ori4paramsall_fixedgaze.(allfields{f}) = ori4params_fixedgaze.(allfields{f});
            end
        end
        end
    end
    
    neucheck = false(Nneuronsall,1);
    
    for iprobe = 1:numel(probes)
        if nnz(floor(unit_peakch/1000)==iprobe-1)==0
            continue
        end
        
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        if validfixgaze
            load(sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe}))
        end
        
        if any(neucheck(neuoind))
            error('some of the neurons seem to already have been accounted for, check')
        end
        neucheck(neuoind)=true;
        
        meanFRall(neuoind) = meanFRvec;
        sponFRall(neuoind) = sponFRvec;
        
        for b = 1:numel(visblocks)
            if ismember(visblocks{b}, ICblocks)
                tloi = psthtli>0 & psthtli<=400;
            else % RFCI, sizeCI
                tloi = psthtli>0 & psthtli<=250;
            end
            Rall.(visblocks{b})(:,neuoind) = 1000*squeeze(mean(psth.(visblocks{b})(tloi,:,:),1));
        end
        
        for b = 1:numel(ICblocks)
            for f= 1:numel(ICsigfields)
                ICsigall.(ICblocks{b}).(ICsigfields{f})(neuoind,:) = ICsig.(ICblocks{b}).(ICsigfields{f});
            end
        end
        
        for f= 1:numel(RFCIfields)
            RFCIall.(RFCIfields{f})(neuoind,:) = RFCI.(RFCIfields{f});
        end
        
        for f= 1:numel(RFCIspinfields)
            RFCIspinall.(RFCIspinfields{f})(neuoind,:) = RFCIspin.(RFCIspinfields{f});
        end
        
        for f= 1:numel(sizeCIfields)
            sizeCIall.(sizeCIfields{f})(neuoind,:) = sizeCI.(sizeCIfields{f});
        end
        
        for f= 1:numel(oriparamsfields)
            oriparamsall.(oriparamsfields{f})(neuoind,:) = oriparams.(oriparamsfields{f});
        end
        
        for f= 1:numel(ori4paramsfields)
            ori4paramsall.(ori4paramsfields{f})(neuoind,:) = ori4params.(ori4paramsfields{f});
        end
        
        if validfixgaze
            for b = 1:numel(ICblocks)
                if ~isempty(fieldnames(ICsig_fixedgaze.(ICblocks{b})))
                for f= 1:numel(ICsigfields)
                    ICsigall_fixedgaze.(ICblocks{b}).(ICsigfields{f})(neuoind,:) = ICsig_fixedgaze.(ICblocks{b}).(ICsigfields{f});
                end
                end
            end
            
            if ~isempty(fieldnames(RFCI_fixedgaze))
                for f= 1:numel(RFCIfields)
                    RFCIall_fixedgaze.(RFCIfields{f})(neuoind,:) = RFCI_fixedgaze.(RFCIfields{f});
                end
            end
            
            if ~isempty(fieldnames(RFCIspin_fixedgaze))
            for f= 1:numel(RFCIspinfields)
                RFCIspinall_fixedgaze.(RFCIspinfields{f})(neuoind,:) = RFCIspin_fixedgaze.(RFCIspinfields{f});
            end
            end
            
            if ~isempty(fieldnames(sizeCI_fixedgaze))
            for f= 1:numel(sizeCIfields)
                sizeCIall_fixedgaze.(sizeCIfields{f})(neuoind,:) = sizeCI_fixedgaze.(sizeCIfields{f});
            end
            end
            
            if ~isempty(fieldnames(oriparams_fixedgaze))
            for f= 1:numel(oriparamsfields)
                oriparamsall_fixedgaze.(oriparamsfields{f})(neuoind,:) = oriparams_fixedgaze.(oriparamsfields{f});
            end
            end
            
            if ~isempty(fieldnames(ori4params_fixedgaze))
            for f= 1:numel(ori4paramsfields)
                ori4paramsall_fixedgaze.(ori4paramsfields{f})(neuoind,:) = ori4params_fixedgaze.(ori4paramsfields{f});
            end
            end
        end
        
    end
    
    if ~all(neucheck)
        error('not all neurons in this recording session were accounted for, check')
    end
    
    save(sprintf('%spostprocessed.mat', pathpp), 'vis', 'Tres', 'Rall', 'neuallloc', '-v7.3')
    
    save(sprintf('%svisresponses.mat', pathpp), 'neuallloc', ...
        'meanFRall', 'sponFRall', 'ICtrialtypes', 'ICsigall', 'RFCIall', 'RFCIspinall', ...
        'sizeCIall', 'oriparamsall', 'ori4paramsall', '-v7.3')
    
    if validfixgaze
        save(sprintf('%svisresponses_fixedgaze%dpix.mat', pathpp, gazedistthresh), 'neuallloc', ...
            'meanFRall', 'sponFRall', 'ICtrialtypes', 'spinmaxdistmodecom', 'spinlikelyblink', ...
            'ICsigall_fixedgaze', 'RFCIall_fixedgaze', 'RFCIspinall_fixedgaze', ...
            'sizeCIall_fixedgaze', 'oriparamsall_fixedgaze', 'ori4paramsall_fixedgaze', '-v7.3')
    end
    toc
end


%% spontaneous block activity (added 241118)
clear all
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130')
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

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

    load([pathpp, 'postprocessed.mat'])

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
    if ~isequal( 1000*electrode_probeid + electrode_localid, ...
        nwb.general_extracellular_ephys_electrodes.id.data.load()) 
        error('check how electrode_id is calculated')
    end
    electrode_location = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();

    % disp(unique(electrode_location)')
    % save([pathpp 'info_electrodes.mat'], 'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')

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
    % save([pathpp 'info_units.mat'], 'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',

    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    neuallloc = electrode_location(revmapelecid(unit_peakch+1));

    % sanity check
    electrodeid = nwb.general_extracellular_ephys_electrodes.id.data.load();
    unit_location = cell(Nneurons,1);
    for ii = 1:Nneurons
        unit_location(ii) = electrode_location(electrodeid==unit_peakch(ii));
    end
    if ~isequal(neuallloc, unit_location)
        error('check how neuallloc is calculated')
    end
    
    %% spontaneous block
    % number of spontaneous chunks
    Nspon = numel(vis.spontaneous_presentations.start_time);
    spondurs = vis.spontaneous_presentations.stop_time-vis.spontaneous_presentations.start_time;
    sponstlens = 1+floor(spondurs/Tres);
    psthspon = cell(Nspon,1);
    for ispon = 1:Nspon

        sponspiketrain = false(sponstlens(ispon), Nneurons);
        for ci = 1:Nneurons
            tempvalspks = spiketimes{ci}>=vis.spontaneous_presentations.start_time(ispon) & ...
                spiketimes{ci}<vis.spontaneous_presentations.stop_time(ispon);
            tempst = spiketimes{ci}(tempvalspks) - vis.spontaneous_presentations.start_time(ispon); % spike timing in seconds
            tempstind = floor(tempst/Tres)+1;
            sponspiketrain(tempstind, ci) = true;
        end

        psthspon{ispon} = sponspiketrain;
    end

    save([pathpp, 'psth_spontaneous.mat'], 'neuallloc', 'vis', 'psthspon', '-v7.3')
    toc(sesclk)
end

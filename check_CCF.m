% addpath(genpath('S:\OpenScopeData\matnwb_HSLabDesktop'))
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\Analize_IC_OpenScope_v230821'))
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

%% probevisareas
probevisareas = cell(Nsessions,numel(probes));
for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions datadir probevisareas probes
    sesclk = tic;
    
    nwbfiles = cat(1, dir([datadir nwbsessions{ises} '\*.nwb']), dir([datadir nwbsessions{ises} '\*\*.nwb']));
    
    % take filename with shortest length or filename that does not contain probe
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
    
    %disp(unique(electrode_location)')
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
%     spiketrain = false(stlen, Nneurons);
%     ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
%     for ii = 1:Nneurons
%         spiketrain(floor(spiketimes{ii}/Tres)+1, ii) = true;
%     end
%     
%     disp(size(spiketrain))
    disp([stlen, Nneurons])
    
    save([pathpp 'info_units.mat'], 'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',
    
    %% probevisareas
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        %         tic
        %         load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        %         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
%         load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        %         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
        %             error('check neuoind')
        %         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        %         disp(probes{iprobe})
        %         disp(unique(neuloc(contains(neuloc, 'VIS')))')
        % probevisareas{ises, iprobe} = sprintf('%s ',unique(neuloc(contains(neuloc, 'VIS'))));
        tempcell = unique(neuloc(contains(neuloc, 'VIS')));
        probevisareas{ises, iprobe} = sprintf('%s ',tempcell{:});
    end
end
open probevisareas

%% report number of units in each area/session/probe
Nsessions = numel(nwbsessions);
% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};

neulocagg = cell(size(probes));
sesneuagg = cell(size(probes));
for ises = 1:Nsessions
    clearvars -except neulocagg sesneuagg probes visareas datadir nwbsessions Nsessions ises
    
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        
        probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
%         disp(unique(probelocs)')
        
        %probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
        neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
        %neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
        %neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
        
        sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
        %sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
    end
end

%%
aggneuloc = cat(1,neulocagg{:});
aggsesneu = cat(1,sesneuagg{:});

% [v,c]=uniquecnt(aggneuloc);
% disp([v c])

v = unique(aggneuloc);
c = zeros(size(v));
for ii = 1:numel(v)
    c(ii) = nnz(strcmp(aggneuloc, v{ii}));
end

areaunitsperses = zeros(length(v), Nsessions);
for ii = 1:length(v)
    neuoi = strcmp(aggneuloc, v(ii));
    hc=histcounts(aggsesneu(neuoi), 0.5:1:Nsessions+0.5);
    if ~( nnz(neuoi)==c(ii) && sum(hc)==c(ii) )
        error('mismatch between uniquecnt and strcmp -- check')
    end
    areaunitsperses(ii,:) = hc;
end

areaunitstab = table(v,c,areaunitsperses);
open areaunitstab

addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visctxareas = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};

lowpassopt = false;
whichneuarea = 'V1';

%% first, get psth of spectral events (vis stim aligned)
for ises = 1:numel(nwbsessions)
    clearvars -except ises datadir nwbsessions probes visareas visctxareas whichneuarea
    
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
    load(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{neuprobe}))
    
    eventvis = struct();
    eventvispsth = struct();
    
    % whichblock = 'spontaneous_presentations';
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
    
    %   TFR - time-frequency response (TFR) (frequency-by-time-trial) for a
    %       single subject/session.
    thrFOM = 6;
    % eventBand = [15 29];
    eventBand = [1 100];
    if contains(whichblock, 'spontaneous')
        [durtrial, itrial] = max(vis.(whichblock).stop_time-vis.(whichblock).start_time);
        TFR = permute(tfrvispsth.(whichblock){itrial}, [2 1]);
        classLabels = itrial;
        %tVec = 1:size(TFR,2);
    else
        TFR = permute(tfrvispsth.(whichblock), [3 1 2]);
        classLabels = reshape(trialorder,1,[]);
        %tVec = psthtli;
    end
    tVec = 1:size(TFR,2);
    if numel(classLabels) ~= size(TFR,3)
        error('TFR third dimension should be number of trials')
    end
    Ntrials = length(classLabels);
    
    % spectralEvents is 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
    % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
    evcols = {'trial_index', 'trial_class', 'maxima_frequency', ...
        'lowerbound_frequency', 'upperbound_frequency', 'frequency_span', ...
        'maxima_timing', 'event_onset_timing', 'event_offset_timing', ...
        'event_duration', 'maxima_power', 'maxima_normmedian_power'};
    tic
    [spectralEvents] = find_localmax_method_2(eventBand, thrFOM, tVec, fVec, ...
        TFR, classLabels);
    toc % 312s spontaneous period: 200s for for beta band, 406s for 1-100Hz
    
    % tic
    % [spectralEvents] = find_localmax_method_1(eventBand, thrFOM, tVec, fVec, ...
    %     TFR, classLabels);
    % toc % 312s spontaneous period: 126s for 1-100Hz
    % % ICwcfg1 block: 155s for beta band
    
    eventvis.(whichblock) = spectralEvents;
    
    % which event occurred when
    eventvispsth.(whichblock) = zeros(length(tVec), Ntrials);
    rowinds = spectralEvents(:,strcmp(evcols, 'maxima_timing'));
    colinds = spectralEvents(:,strcmp(evcols, 'trial_index'));
    rcinds = sub2ind([length(tVec), Ntrials], rowinds, colinds);
    eventvispsth.(whichblock)(rcinds) = 1:size(spectralEvents,1);
    
    %{
evfreq = spectralEvents(:,strcmp(evcols, 'maxima_frequency'));
betaeventinds = find(evfreq>=15 & evfreq<30);
betaeventpsth = ismember(eventvispsth.(whichblock), betaeventinds);

figure;
trialsoi = trialorder==106;
plot(psthtli, mean(betaeventpsth(:,trialsoi),2))

evfreq = spectralEvents(:,strcmp(evcols, 'maxima_frequency'));
evtime = spectralEvents(:,strcmp(evcols, 'maxima_timing'));
evoi = psthtli(evtime)>=0 & psthtli(evtime)<400 & evfreq>=15 & evfreq<30;
[v,c]=uniquecnt(spectralEvents(evoi,strcmp(evcols, 'trial_class')));
evhist = zeros(1,numel(vistrialtypes));
evhist( ismember(vistrialtypes, v) ) = c;
figure; hold all
plot(1:numel(vistrialtypes), evhist, 'o-')
set(gca, 'XTick', 1:numel(vistrialtypes), 'XTickLabel', vistrialtypes)
    %}
    
    save( sprintf('%s%s_eventvis%s_probe%s.mat', pathpp, whichneuarea, blockname, probes{neuprobe}), ...
        'eventBand', 'evcols', 'eventvis', 'eventvispsth', '-v7.3')
    
    
    ettrange = -400:400;
    evfreq = spectralEvents(:,strcmp(evcols, 'maxima_frequency'));
    evtime = spectralEvents(:,strcmp(evcols, 'maxima_timing'));
    evtrialtype = spectralEvents(:,strcmp(evcols, 'trial_class'));
    betaetvisspike = cell(size(vistrialtypes));
    for typi = 1:numel(vistrialtypes)
        tic
        trialsoi = trialorder==vistrialtypes(typi);
        evoind = find( evtrialtype==vistrialtypes(typi) & psthtli(evtime)>=0 & psthtli(evtime)<400 & evfreq>=15 & evfreq<30 );
        evpsthind = find(ismember(reshape(eventvispsth.(whichblock),[],1), evoind));
        evtind = evpsthind+ettrange;
        
        betaetvisspike{typi} = NaN(Nneu, length(ettrange));
        for ci = 1:Nneu
            temppsth = reshape(psth.(whichblock)(:,:,ci),[],1);
            betaetvisspike{typi}(ci,:) = mean( temppsth(evtind), 1);
        end
        toc
    end
    
    % all events in the block
    evfreq = spectralEvents(:,strcmp(evcols, 'maxima_frequency'));
    evtime = spectralEvents(:,strcmp(evcols, 'maxima_timing'));
    evoind = find( psthtli(evtime)>=0 & psthtli(evtime)<800 & evfreq>=15 & evfreq<30 );
    evpsthind = find(ismember(reshape(eventvispsth.(whichblock),[],1), evoind));
    evtind = evpsthind+ettrange;
    betaetspikeblock = NaN(Nneu, length(ettrange));
    for ci = 1:Nneu
        temppsth = reshape(psth.(whichblock)(:,:,ci),[],1);
        betaetspikeblock(ci,:) = mean( temppsth(evtind), 1);
    end
    clearvars evtind evpsthind
    
    %{
figure; plot(ettrange, etspikeblock)
xlim([-150 150])
kerwinhalf = 12; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));
figure; plot(ettrange, convn(etspikeblock, kergauss, 'same'))
xlim([-150 150])
figure; imagesc(ettrange, 1:Nneu, convn(etspikeblock, kergauss, 'same')); xlim([-150 150])
    %}
    
    save( sprintf('%s%s_betaetspike%s_probe%s.mat', pathpp, whichneuarea, blockname, probes{neuprobe}), ...
        'vistrialtypes', 'ettrange', 'betaetvisspike', 'betaetspikeblock', '-v7.3')
    
end

%% basic spectral event analyses
% are there clear "bands" of spectral event frequencies? (boundary between beta and gamma
% beta events are characterized by a prominent supragranuar sink -- do gamma events have a granular sink? and what about alpha events?

%% in continuous data, calculate event triggered

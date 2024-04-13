%{
node strength can be calculated with:
-CRF edge potential
-Pearson correlation
-Okun's population coupling
-mutual information between pairs of neurons

sensory sensitivity can be cauculated with:
-CRF AUC
-SP_ICvsRC
-mutual information between vis stim and neural responses
%}
% IMPORTANT NOT TO DO GENPATH
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
ses2anal = [1 6];
for ises = ses2anal
    clearvars -except datadir nwbsessions ses2anal ises
    %disp(nwbsessions{ises})
    fprintf('%d %s\n', ises, nwbsessions{ises})
    sesclk = tic;
    
    neuopt = 'ctx';
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed_probeC.mat'], 'vis')
    
    pathsv = [datadir 'CCG' filesep nwbsessions{ises} filesep];
    if ~exist(pathsv, 'dir')
        mkdir(pathsv)
    end
    
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location'
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data'
%     load([pathpp 'spiketimes.mat'])
    load([pathpp 'visresponses.mat'])
    % isequal(ststartend, [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1])
    Tres = 0.001; % 1ms
    
    %% extract spike times
    neuctx = contains(neuallloc, 'VIS');
    neuctxind = find(neuctx);
    
    Nneuctx = numel(neuctxind);
    neulocctx = neuloc(neuctxind);
    % isequal(neulocctx, find(contains(neuloc, 'VIS')))
    
    recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
    visarealabels = zeros(Nneuctx,1);
    for a = 1:numel(recarealabels)
        if strcmp(recarealabels{a}, 'VISp')
            neuinarea = contains(neulocctx, 'VISp') & ~contains(neulocctx, 'VISpm');
        else
            neuinarea = contains(neulocctx, recarealabels{a});
        end
        visarealabels(neuinarea) = a;
    end
    
    ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
    ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
    ICblockstartend = floor(ICblockstartend/Tres)+1;
    
    Nneurons = length(unit_ids);
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
    save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuallloc', 'neuctx', 'neuctxind', 'neulocctx', ...
        'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')
    
    %%    
    switch neuopt
        case 'ctx'
            neuctx = contains(neuallloc, 'VIS');
        case 'ctrctx'
            neuctx = contains(neuallloc, 'VIS') & RFCIall.RFindclassic==1 & RFCIall.Pkw_rfclassic<0.05;
        case 'ctxL23'
            neuctx = contains(neuallloc, 'VIS') & contains(neuallloc, '2/3');
    end
    neuctxind = find(neuctx) ;
    
    %%
    ctxspiketrain = false(Nneuctx, ststartend(end) );
    for ii = 1:Nneuctx
        ci = neuctxind(ii);
        ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
    end
    
    ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));
    % ctxspiketrain = ctxspiketrain(:,ststartend(1):ststartend(2));
    
    stlen = size(ctxspiketrain, 2);
    spkcntvec = sum(ctxspiketrain,2);
    sqrtspkcntmat = sqrt( spkcntvec * spkcntvec' );

    
    %%
    NFFT = 2^ceil(log2(stlen));
    CCGffttl = -NFFT/2:NFFT/2-1;
    Thalfwin = 100;
    CCGtli_fft = -Thalfwin:Thalfwin;
    fprintf('Number of neurons: %d\n', Nneuctx)
    
    tic
    ctxCCG_fft = NaN(Nneuctx, Nneuctx, length(CCGtli_fft));
    for ci = 1:Nneuctx
        for i100 = 1:ceil((Nneuctx-ci+1)/100)
            tempneuvec = ci+100*(i100-1): min([Nneuctx ci-1+100*(i100)]);
            %disp([tempneuvec(1) tempneuvec(end)])
            tempCCGfft = fftshift(ifft(fft(ctxspiketrain(tempneuvec,:)', NFFT) .* conj(fft(ctxspiketrain(ci,:)', NFFT))),1);
            ctxCCG_fft(ci,tempneuvec,:) = tempCCGfft(ismember(CCGffttl, CCGtli_fft),:)';
            ctxCCG_fft(tempneuvec,ci,:) = tempCCGfft(ismember(CCGffttl, CCGtli_fft),:)';
        end
        if mod(ci,10)==1
            toc
            fprintf('done with neuron #%d\n', ci)
        end
    end
    ctxCCGfft = ctxCCG_fft./sqrtspkcntmat;
    toc
    disp('calculated ctxCCG')
    
    switch neuopt
        case 'ctx'
            save([pathsv 'ctxCCGfft.mat'], 'stlen', 'spkcntvec', 'neuctx', 'neulocctx', 'visarealabels', 'CCGtli_fft', 'ctxCCGfft', '-v7.3')
        case 'ctrctx'
            neuctrctx = neuctx;
            neulocctrctx = neulocctx;
            ctrctxCCGfft = ctxCCGfft;
            ctrctxCCGweight = ctxCCGweight;
            save([pathsv 'ctrctxCCGfft.mat'], 'stlen', 'spkcntvec', 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli_fft', 'ctrctxCCGfft', '-v7.3')
        case 'ctxL23'
            neuctxL23 = neuctx;
            neulocctxL23 = neulocctx;
            ctxL23CCGfft = ctxCCGfft;
            ctxL23CCGweight = ctxCCGweight;
            save([pathsv 'ctxL23CCGfft.mat'], 'stlen', 'spkcntvec', 'neuctxL23', 'neulocctxL23', 'visarealabels', 'CCGtli_fft', 'ctxL23CCGfft', '-v7.3')
    end
    
    
    toc(sesclk)
end

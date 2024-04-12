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
    %%
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location'
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data'
    load([pathpp 'spiketimes.mat'])
    load([pathpp 'visresponses.mat'])
    % isequal(ststartend, [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1])
    Tres = 0.001; % 1ms
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    electrode_location = cellstr(electrode_location);
    neuloc = electrode_location(revmapelecid(unit_peakch+1));
    
    
    switch neuopt
        case 'ctx'
            neuctx = contains(neuloc, 'VIS');
        case 'ctrctx'
            neuctx = contains(neuloc, 'VIS') & RFCIall.RFindclassic==1 & RFCIall.Pkw_rfclassic<0.05;
        case 'ctxL23'
            neuctx = contains(neuloc, 'VIS') & contains(neuloc, '2/3');
    end
    neuctxind = find(neuctx) ;
    
    %%
    
    % neuctx = [];
    % neuctxind = [];
    % for iprobe = 1:6
    % tempneuoind = find(floor(unit_peakch/1000)==iprobe-1);
    % tempneuloc = electrode_location(revmapelecid(unit_peakch(tempneuoind)+1));
    % neuctx = cat(1, neuctx, contains(tempneuloc, 'VIS'));
    % neuctxind = cat(1, neuctxind, tempneuoind(contains(tempneuloc, 'VIS')));
    % end
    
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
    
    %%
    ctxspiketrain = false(Nneuctx, ststartend(end) );
    for ii = 1:Nneuctx
        ci = neuctxind(ii);
        ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
    end
    
    ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
    ICblockstartend = floor(ICblockstartend/Tres)+1;
    ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));
    % ctxspiketrain = ctxspiketrain(:,ststartend(1):ststartend(2));
    
    stlen = size(ctxspiketrain, 2);
    spkcntvec = sum(ctxspiketrain,2);
    sqrtspkcntmat = sqrt( spkcntvec * spkcntvec' );

    % save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuloc', 'neuctx', 'neuctxind', 'neulocctx', ...
    %     'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')
    
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

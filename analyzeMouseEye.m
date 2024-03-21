addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\Analize_IC_OpenScope_v240130'))
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

% Skipping session1 sub-619293 -- EyeTracking error in nwb
% Skipping session12 sub-637484 -- EyeTracking error in nwb

%%
for ises = 1:Nsessions
    clearvars -except ises Nsessions nwbsessions datadir
    sesclk = tic;
    
    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};

    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed_probeC.mat'], 'vis')
    
    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));    
    % take filename with shortest length or filename that does not contain probe
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    % nwbspikefile = string(nwbspikefile);
    disp(nwbspikefile)
    nwb = nwbRead(nwbspikefile); %, 'ignorecache');
    
    try
        TrackEye = nwb.acquisition.get('EyeTracking');
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    
    TrackEyeTimestamps = nwb.acquisition.get('EyeTracking').eye_tracking.timestamps.load();
    % figure; histogram(diff(TrackEyeTimestamps))
    
    likelyblink = nwb.acquisition.get('EyeTracking').likely_blink.data.load();
    
    % PupilTracking Eye-tracking data, representing pupil size
    % EyeTracking Eye-tracking data, representing direction of gaze
    % angle refers to phi, which is the rotation fit to the height
    % NaNs occur when DeepLabCut didn’t output enough points to fit an ellipse.
    % all positions (width, height, center_x, center_y) are in units of pixels
    % all areas are in units of pixels2
    eyefields = {'angle', 'area', 'area_raw', 'height', 'width', 'data'};
    eyetracking = struct();
    pupiltracking = struct();
    cornealreflection = struct();
    for e = 1:numel(eyefields)
        eyetracking.(eyefields{e}) = nwb.acquisition.get('EyeTracking').eye_tracking.(eyefields{e}).load();
        pupiltracking.(eyefields{e}) = nwb.acquisition.get('EyeTracking').pupil_tracking.(eyefields{e}).load();
        cornealreflection.(eyefields{e}) = nwb.acquisition.get('EyeTracking').corneal_reflection_tracking.(eyefields{e}).load();
    end
    
    % bin 3 pixels
    % no need to normalize since camera magnification settings are standardized
    binpix = 3;
    % figure; hold all
    % histogram2(pupiltracking.data(1,:), pupiltracking.data(2,:), 'binwidth', binpix, 'displaystyle', 'tile')
    [N,XEDGES,YEDGES] = histcounts2(pupiltracking.data(1,:), pupiltracking.data(2,:), 'binwidth', binpix);
    % 1st row corresponds to X, which correspond to rows of N
    xbinctrs = (XEDGES(1:end-1)+XEDGES(2:end))/2;
    ybinctrs = (YEDGES(1:end-1)+YEDGES(2:end))/2;
    [mv,mi]=max(N,[],'all', 'linear');
    [r,c]=find(N==max(N,[],'all'));
    if (c-1)*size(N,1)+r ~= mi
        error('check max point in histogram')
    end
    disp([xbinctrs(r) ybinctrs(c)])
    
    modecom = [xbinctrs(r) ybinctrs(c)];
    distmodecom = sqrt(sum((pupiltracking.data'-modecom).^2,2));
    % when distmodecom threshold is set to 10 pixels, 50% of the trials are preserved
    
    if ~isequal(size(distmodecom), size(likelyblink))
        error('timestamps ot shared between likely_blink, eye_tracking and pupil_tracking?')
    end
    
    % %% trials should start -0.5s before and end 1s after stim onset
    % 1/nanmedian(diff(TrackEyeTimestamps)) roughly 60 Hz frame rate
    trialdistmodecom = struct();
    trialmaxdistmodecom = struct();
    triallikelyblink = struct();
    trackeyetli = -30:60;
    tic
    for b = 1:numel(visblocks)
        if contains(visblocks{b}, 'spontaneous')
            continue
        end
        [r,c]=find(cumsum(TrackEyeTimestamps-vis.(visblocks{b}).trialstart'>0,1)==1);
        if ~isequal(c, (1:numel(vis.(visblocks{b}).trialstart))' )
            error('missing some trials')
        end
        trackeyetrialinds = (r-1)+trackeyetli;
        trackeyepsth = distmodecom(trackeyetrialinds);
        
        trialdistmodecom.(visblocks{b}).trackeyetli = trackeyetli;
        trialdistmodecom.(visblocks{b}).psthtrialinds = trackeyetrialinds;
        trialdistmodecom.(visblocks{b}).psth = trackeyepsth;
        
        likelyblinkpsth = likelyblink(trackeyetrialinds);
        
        % for IC blocks psthtli>0 & psthtli<=400
        % for RFCI blocks psthtli>0 & psthtli<=1000
        % for sizeCI blocks psthtli>0 & psthtli<=250
        if contains(visblocks{b}, 'IC')
            endframe = round(0.4*60);
        elseif contains(visblocks{b}, 'RFCI')
            endframe = round(1*60);
        elseif contains(visblocks{b}, 'sizeCI')
            endframe = round(0.25*60);
        else
            error('visblock not recognized')
        end
        
        trialmaxdistmodecom.(visblocks{b}) = max(trackeyepsth(:,trackeyetli>=0 & trackeyetli<endframe),[],2);
        triallikelyblink.(visblocks{b}) = any(likelyblinkpsth(:,trackeyetli>=0 & trackeyetli<endframe), 2);
    end
    toc % takes 6 min for 2500 units
    
    save([pathpp 'trackmouseeye.mat'], 'TrackEyeTimestamps', 'likelyblink', ...
        'eyetracking', 'pupiltracking', 'cornealreflection', 'modecom', 'distmodecom', ...
        'trialdistmodecom', 'trialmaxdistmodecom', 'triallikelyblink', '-v7.3')
    toc
end

%% aggregate across sessions to determine appropriate threshold

% figure; histogram(distmodecom, 'normalization', 'probability')
% xlabel('Distance (pixels)')
%
% mean(distmodecom<10)
% % 0.05*2*sqrt(nanmean(pupiltracking.area)/pi) correspond to ~3 pixels
% % pupil diameter is roughly 70 pixels, eye diameter is roughly 270 pixels
% % nanmedian(eyetracking.width) is ~155 pixels.
% width is defined as horizontal halfaxis of the ellipse fit of the pupil. 
% see https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/4e/be/4ebe2911-bd38-4230-86c8-01a86cfd758e/visual_behavior_2p_technical_whitepaper.pdf
% % this corresponds to mouse's gaze range, which is ~160 degrees
% % diameter of the mouse eye is eyewidth/cosd((180-gazerange)/2), i.e., ~160 pixels
% % circumference of eye = eyediameter*pi (~160*pi = ~500 pixels)
% % corresponds to 180 degrees
% % 10 pixels roughly correspond to 3.6 degrees (~50% timepoints)
% % 7 pixels roughly correspond to 2.5 degrees (~30% timepoints)
% % 22 pixels roughly correspond to 8 degrees

% figure; plot(eyetracking.area, pi*eyetracking.width.*eyetracking.height, '.')
gazerange = 160;

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
trialdistmodecomagg = struct();
trialmaxdistmodecomagg = struct();
triallikelyblinkagg = struct();
for b = 1:numel(visblocks)
    trialdistmodecomagg(Nsessions).(visblocks{b}) = [];
    trialmaxdistmodecomagg(Nsessions).(visblocks{b}) = [];
    triallikelyblinkagg(Nsessions).(visblocks{b}) = [];
end
pupilwidth = NaN(1, Nsessions);
eyewidth = NaN(1, Nsessions);
pix8visdeg = NaN(1, Nsessions);
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    try
        load([pathpp 'trackmouseeye.mat'])
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    pupilwidth(ises) = nanmedian(pupiltracking.width);
    eyewidth(ises) = nanmedian(eyetracking.width);
    tempeyediameter = eyewidth(ises)/cosd((180-gazerange)/2);
    tempeyecircumference = pi*tempeyediameter;
    pix8visdeg(ises) = tempeyecircumference*8/180;

    trialdistmodecomagg(ises) = trialdistmodecom;
    trialmaxdistmodecomagg(ises) = trialmaxdistmodecom;
    triallikelyblinkagg(ises) = triallikelyblink;
end

fs = 14;
figure
for ises = 1:Nsessions
    subplot(3,4,ises); hold all
    histogram(trialmaxdistmodecomagg(ises).ICwcfg1_presentations, 'Normalization', 'cdf')
    yl = ylim;
    plot(10*[1 1], yl, 'r--')
    %plot(7*[1 1], yl, 'r--')
    plot(pix8visdeg(ises)*[1 1], yl, 'r--')
    ylim(yl)
    xlim([0 50])
    set(gca, 'FontSize', fs)
    xlabel('camera pixels')
    ylabel('ICwcfg1 trial CDF')
    title(sprintf('%s: 8 vis deg ~%.0f pix\nwidth pupil %.0f eye %.0f', ...
        nwbsessions{ises}, pix8visdeg(ises), pupilwidth(ises), eyewidth(ises) ))
end

figure; histogram(pix8visdeg, 'binwidth', 1)
xlabel('pixels corresponding to 8 visual degrees')
ylabel('# Sessions')
title(sprintf('min %.1f mean %.1f median %.1f', min(pix8visdeg), nanmean(pix8visdeg), nanmedian(pix8visdeg)))

%% test gazedistthresh with RFCI blocks
% note by HS 240226
% 9/12 sessions are valid with gazedistthresh of 20 (roughly 8 visual degrees)

visagg = struct();
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
    if ises==1
        visagg=vis;
    else
        visagg=cat(1, visagg, vis);
    end
end

gazedistthresh = 20;
propfixedgazetrials = zeros(Nsessions,1);
propfixcrftrials = zeros(Nsessions,1);
figure
for ises = 1:Nsessions
    subplot(3,4,ises); hold all
    histogram(trialmaxdistmodecomagg(ises).RFCI_presentations, 'binwidth', 3, 'Normalization', 'cdf')
    yl = ylim;
    plot(gazedistthresh*[1 1], yl, 'r--')
    ylim(yl)
    xlim([0 50])
    set(gca, 'FontSize', fs)
    xlabel('camera pixels')
    ylabel('RFCI trial CDF')
    
    if ~isempty(trialmaxdistmodecomagg(ises).RFCI_presentations)
    temptrialorder = visagg(ises).RFCI_presentations.trialorder(1:4:end);
    tempgazedist = trialmaxdistmodecomagg(ises).RFCI_presentations(1:4:end);
    templikelyblink = triallikelyblinkagg(ises).RFCI_presentations(1:4:end);
    temptrialsfixedgaze = tempgazedist<gazedistthresh & ~templikelyblink;
    propfixedgazetrials(ises) = mean(temptrialsfixedgaze);
    propfixcrftrials(ises) = mean(temptrialsfixedgaze(floor(temptrialorder/10000) == 0));

    validRFCIfix = true;
    fixcrftrials = temptrialsfixedgaze & floor(temptrialorder/10000) == 0;% &
    if ~isequal( unique(floor(mod(temptrialorder(fixcrftrials), 1000) / 10)), (0:8)' )
        validRFCIfix = false;
    end
    else
        validRFCIfix = false;
    temptrialorder = visagg(ises).RFCI_presentations.trialorder(1:4:end);
    temptrialsfixedgaze = false(size(temptrialorder));
    end

    if validRFCIfix
        titcol = [0 0 1];
    else
        titcol = [0.5 0.5 0.5];
    end
    % title(sprintf('validRFCI %d %s: 8 vis deg ~%.0f pix\nwidth pupil %.0f eye %.0f', ...
    %     validRFCIfix, nwbsessions{ises}, pix8visdeg(ises), pupilwidth(ises), eyewidth(ises) ), 'Color', titcol)
    title(sprintf('validRFCI %d %s: 8 vis deg ~%.0f pix\nfixed-gaze trials %.0f%%', ...
        validRFCIfix, nwbsessions{ises}, pix8visdeg(ises), 100*mean(temptrialsfixedgaze) ), 'Color', titcol)
end

%% fixed gaze psthall and Rall
% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
% visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
% visind = [6 5 1 2 4 3];

for ises = 1:Nsessions
    clearvars -except ises nwbsessions Nsessions datadir
    sesclk = tic;
    fprintf('\nSession %d %s\n', ises, nwbsessions{ises})
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_units.mat'])
    try
        load([pathpp 'trackmouseeye.mat'])
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
    gazedistthresh = 20;
    
    trackeyetli = trialdistmodecom.RFCI_presentations.trackeyetli;
    trackeyetrialinds = trialdistmodecom.RFCI_presentations.psthtrialinds;
    trackeyepsth = trialdistmodecom.RFCI_presentations.psth;
    likelyblinkpsth = likelyblink(trackeyetrialinds);
    % for RFCIspin psthtli>0 & psthtli<=250
    endframe = round(0.25*60);
    spinmaxdistmodecom = max(trackeyepsth(:,trackeyetli>=0 & trackeyetli<endframe),[],2);
    spinlikelyblink = any(likelyblinkpsth(:,trackeyetli>=0 & trackeyetli<endframe), 2);
    
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
        
        vrfgfn= sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe});
        if exist( vrfgfn, 'file')
            fprintf('%s exists, skipping...\n', vrfgfn)
            continue
            %load(vrfgfn)
        end
        
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe})); %, 'meanFRvec', 'sponFRvec')
        
        %{
tloi = psthtli>0 & psthtli<=1000;
tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,1:4:end,:), 1))';
temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
RFCI = analyzeRFCI(tempR, temptrialorder, sponFRvec);
fprintf('ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
    nnz(RFCI.RFindclassic==1 & RFCI.pRFclassic<0.05), ...
    nnz(RFCI.RFindclassic==1 & RFCI.Pkw_rfclassic<0.05), ...
    nnz(RFCI.RFindclassic==1 & RFCI.sigmc_rfclassic), ...
    nnz(RFCI.RFindclassic==1 & RFCI.RFsigexclclassic), ...
    nnz(RFCI.RFindclassic==1 & RFCI.RFexclsigclassic) )

tloi = psthtli>0 & psthtli<=250;
tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,:,:), 1))';
temptrialorder = vis.RFCI_presentations.trialorder;
RFCIspin = analyzeRFCIspin(tempR, temptrialorder, sponFRvec);

        %{
% sanity check fails because of minor variations in duration of stim presentation
fields2check = {'Rrfclassic', 'Rrfinverse', 'RFindclassic', 'RFindinverse'};
checkspin = true;
for f = 1:numel(fields2check)
checkspin = checkspin && isequal(RFCI.(fields2check{f}), RFCIspin.(fields2check{f}));
end
if ~checkspin
    error('inconsistency between RFCI and RFCIspin: check')
end
        %}
        %{
% note, R4 and R1000 mostly match except when stim presentation deviated
% from exactly 250ms
tloi = psthtli>0 & psthtli<=250;
R250 = squeeze(1000*mean(psth.RFCI_presentations(tloi,:,:), 1))';
tloi = psthtli>0 & psthtli<=1000;
R1000 = squeeze(1000*mean(psth.RFCI_presentations(tloi,:,:), 1))';
R4 = squeeze(mean(reshape(tempR,numel(neuoind),4,[]),2));
figure; plot(R4, R1000(:,1:4:end), 'o')
        %}

fprintf('spin ctrCRF: pRFclassic* %d Pfried_rfclassic* %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
    nnz(RFCIspin.RFindclassic==1 & RFCIspin.pRFclassic<0.05), ...
    nnz(RFCIspin.RFindclassic==1 & RFCIspin.Pfried_rfclassic<0.05), ...
    nnz(RFCIspin.RFindclassic==1 & RFCIspin.RFsigexclclassic), ...
    nnz(RFCIspin.RFindclassic==1 & RFCIspin.RFexclsigclassic) )

save(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}), ...
    'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'RFCIspin', 'sizeCI', 'oriparams', '-v7.3')
        %}
        
        % %%
        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        ICblocks = 1:4;
        ICsig_fixedgaze = struct();
        for b = ICblocks
            disp([visblocks{b} ' fixed gaze'])
            tloi = psthtli>0 & psthtli<=400;
            tempR = squeeze(1000*mean(psth.(visblocks{b})(tloi,:,:), 1))';
            temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
            temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
            validICfix = all(ismember(ICtrialtypes, unique(temptrialorder(temptrialsfixedgaze))));
            if validICfix
                ICsig_fixedgaze.(visblocks{b}) = analyzeStaticICtxi(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze) );
            else
                ICsig_fixedgaze.(visblocks{b}) = struct();
                disp('skipped')
            end
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
        tempgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
        templikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
        temptrialsfixedgaze = tempgazedist<gazedistthresh & ~templikelyblink;
        tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,1:4:end,:), 1))';
        temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
        
        validRFCIfix = true;
        fixcrftrials = temptrialsfixedgaze & floor(temptrialorder/10000) == 0;% &
        if ~isequal( unique(floor(mod(temptrialorder(fixcrftrials), 1000) / 10)), (0:8)' )
            validRFCIfix = false;
        end
        fixirftrials = temptrialsfixedgaze & floor(temptrialorder/10000) == 1;% & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
        % if ~isequal( unique(floor(mod(temptrialorder(fixirftrials), 1000) / 10)), (0:8)' )
        %     validRFCIfix = false;
        % end
        if validRFCIfix
            RFCI_fixedgaze = analyzeRFCI(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze), sponFRvec);
            fprintf('fixed gaze ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
                nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.pRFclassic<0.05), ...
                nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.Pkw_rfclassic<0.05), ...
                nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.sigmc_rfclassic), ...
                nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.RFsigexclclassic), ...
                nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.RFexclsigclassic) )
        else
            RFCI_fixedgaze = struct();
            disp('fixedgaze RFCI block skipped')
        end
        
        tloi = psthtli>0 & psthtli<=250;
        temptrialsfixedgaze = spinmaxdistmodecom<gazedistthresh & ~spinlikelyblink;
        tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,:,:), 1))';
        temptrialorder = vis.RFCI_presentations.trialorder;
        
        validRFCIspinfix = isequal( unique(temptrialorder), unique(temptrialorder(temptrialsfixedgaze)) );
        
        if validRFCIspinfix
            RFCIspin_fixedgaze = analyzeRFCIspin(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze), sponFRvec);
            fprintf('fixed gaze spin ctrCRF: pRFclassic* %d Pfried_rfclassic* %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
                nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.pRFclassic<0.05), ...
                nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.Pfried_rfclassic<0.05), ...
                nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.RFsigexclclassic), ...
                nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.RFexclsigclassic) )
        else
            RFCIspin_fixedgaze = struct();
            disp('fixedgaze RFCIspin block skipped')
        end
        
        
        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        %szvec = [0, 4, 8, 16, 32, 64];
        % durstim = vis.sizeCI_presentations.stop_time-vis.sizeCI_presentations.start_time;
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
        % durstim = vis.sizeCI_presentations.start_time(2:end)-vis.sizeCI_presentations.stop_time(1:end-1);
        % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
        tloi = psthtli>0 & psthtli<=250;
        tempR = squeeze(1000*mean(psth.sizeCI_presentations(tloi,:,:), 1))';
        temptrialorder = vis.sizeCI_presentations.trialorder;
        temptrialsfixedgaze = trialmaxdistmodecom.sizeCI_presentations<gazedistthresh & ~triallikelyblink.sizeCI_presentations;
        [sizeCI_fixedgaze, oriparams_fixedgaze, ori4params_fixedgaze] = analyzesizeCI(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze) );
        
        
        save(vrfgfn, 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'spinmaxdistmodecom', 'spinlikelyblink', ...
            'ICsig_fixedgaze', 'RFCI_fixedgaze', 'RFCIspin_fixedgaze', ...
            'sizeCI_fixedgaze', 'oriparams_fixedgaze', 'ori4params_fixedgaze', '-v7.3')
        
    end
end

%% report number of ctrCRF neurons with different significance tests
for ises = 1:Nsessions
    %clearvars -except ises nwbsessions datadir
    fprintf('\nSession %d %s\n', ises, nwbsessions{ises})
    
    gazedistthresh = 20;
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed_probeC.mat'], 'vis')
    try
        load([pathpp 'trackmouseeye.mat'])
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    
    load(sprintf('%svisresponses_probeC.mat', pathpp))
    fprintf('ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
        nnz(RFCI.RFindclassic==1 & RFCI.pRFclassic<0.05), ...
        nnz(RFCI.RFindclassic==1 & RFCI.Pkw_rfclassic<0.05), ...
        nnz(RFCI.RFindclassic==1 & RFCI.sigmc_rfclassic), ...
        nnz(RFCI.RFindclassic==1 & RFCI.RFsigexclclassic), ...
        nnz(RFCI.RFindclassic==1 & RFCI.RFexclsigclassic) )
    
    fprintf('spin ctrCRF: pRFclassic* %d Pfried_rfclassic* %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
        nnz(RFCIspin.RFindclassic==1 & RFCIspin.pRFclassic<0.05), ...
        nnz(RFCIspin.RFindclassic==1 & RFCIspin.Pfried_rfclassic<0.05), ...
        nnz(RFCIspin.RFindclassic==1 & RFCIspin.RFsigexclclassic), ...
        nnz(RFCIspin.RFindclassic==1 & RFCIspin.RFexclsigclassic) )
    
    load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh))
    if ~isempty(fieldnames(RFCI_fixedgaze))
        fprintf('fixed gaze ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
            nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.pRFclassic<0.05), ...
            nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.Pkw_rfclassic<0.05), ...
            nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.sigmc_rfclassic), ...
            nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.RFsigexclclassic), ...
            nnz(RFCI_fixedgaze.RFindclassic==1 & RFCI_fixedgaze.RFexclsigclassic) )
    else
        disp('fixedgaze RFCI block skipped')
    end
    
    if ~isempty(fieldnames(RFCIspin_fixedgaze))
        fprintf('fixed gaze spin ctrCRF: pRFclassic* %d Pfried_rfclassic* %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
            nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.pRFclassic<0.05), ...
            nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.Pfried_rfclassic<0.05), ...
            nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.RFsigexclclassic), ...
            nnz(RFCIspin_fixedgaze.RFindclassic==1 & RFCIspin_fixedgaze.RFexclsigclassic) )
    else
        RFCIspin_fixedgaze = struct();
        disp('fixedgaze RFCIspin block skipped')
    end
    
end

%% report number of nonsigkwBI
% prediction: sigkwBI is a subset of sigkwBK
% non-sigkwBI should do inference,
% but non-sigkwBK should not be able to do inference
for ises = 1:Nsessions
    %clearvars -except ises nwbsessions datadir
    fprintf('Session %d %s\n', ises, nwbsessions{ises})
    
    gazedistthresh = 20;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    try
        load([pathpp 'trackmouseeye.mat'])
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    iprobe = 3;
    load([pathpp 'postprocessed_probeC.mat'], 'vis')
    load(sprintf('%svisresponses_probeC.mat', pathpp), 'ICsig')
    load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh), 'ICsig_fixedgaze')
    neuoind = find(floor(unit_peakch/1000)==iprobe-1);
    neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
    neuctx = contains(neuloc, 'VIS');
    whichvisarea = unique(neuloc(neuctx));
    
    ICblocks = fieldnames(ICsig);
    
    % expect sigkwBK & ~sigkwBI to be ICencoders/RCencoders
    % expect almost all of sigkwBI to be sigkwBK
    % expect almost all of nonsigkwBK to be nonsigkwBI
    for b = 1:numel(ICblocks)
        sigkwBI = ICsig.(ICblocks{b}).PkwBI(neuctx)<0.05;
        sigkwBK = ICsig.(ICblocks{b}).PkwBK(neuctx)<0.05;
        fprintf('%s Nneurons %d, sigkwBI %d, sigkwBK %d, sigkwBK&~sigkwBI %d, IC-encoder %d, RC-encoder %d\n', ...
            ICblocks{b}, ...
            nnz(neuctx), ...
            nnz(ICsig.(ICblocks{b}).PkwBI(neuctx)<0.05), ...
            nnz(ICsig.(ICblocks{b}).PkwBK(neuctx)<0.05), ...
            nnz(sigkwBK&~sigkwBI), ...
            nnz(ICsig.(ICblocks{b}).ICencoder(neuctx)), ...
            nnz(ICsig.(ICblocks{b}).RCencoder(neuctx)) )
        fprintf('%.2f%% sigkwBK/sigkwBI, %.2f%% nonsigkwBI/nonsigkwBK\n', ...
            100*mean(sigkwBK(sigkwBI)), 100*mean(~sigkwBI(~sigkwBK)) )
        
        if ~isempty(fieldnames(ICsig_fixedgaze.(ICblocks{b})))
            sigkwBI = ICsig_fixedgaze.(ICblocks{b}).PkwBI(neuctx)<0.05;
            sigkwBK = ICsig_fixedgaze.(ICblocks{b}).PkwBK(neuctx)<0.05;
            fprintf('fixed-gaze Nneurons %d, sigkwBI %d, sigkwBK %d, sigkwBK&~sigkwBI %d, IC-encoder %d, RC-encoder %d\n', ...
                nnz(neuctx), ...
                nnz(ICsig_fixedgaze.(ICblocks{b}).PkwBI(neuctx)<0.05), ...
                nnz(ICsig_fixedgaze.(ICblocks{b}).PkwBK(neuctx)<0.05), ...
                nnz(sigkwBK&~sigkwBI), ...
                nnz(ICsig_fixedgaze.(ICblocks{b}).ICencoder(neuctx)), ...
                nnz(ICsig_fixedgaze.(ICblocks{b}).RCencoder(neuctx)) )
            fprintf('%.2f%% sigkwBK/sigkwBI, %.2f%% nonsigkwBI/nonsigkwBK\n', ...
                100*mean(sigkwBK(sigkwBI)), 100*mean(~sigkwBI(~sigkwBK)) )
        else
            disp('fixedgaze IC block skipped')
        end
    end
end

clear all; close all; clc
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\Analize_IC_OpenScope_v240130'))
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

%% 240317 reran analyzeRFCI for all trials and fixed gaze to add one-tail stats

for ises = 1:Nsessions
    clearvars -except ises nwbsessions Nsessions datadir
    sesclk = tic;
    fprintf('\nSession %d %s\n', ises, nwbsessions{ises})

    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
    probes = {'A', 'B', 'C', 'D', 'E', 'F'};

    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_units.mat'])
    if exist([pathpp 'trackmouseeye.mat'], 'file')
        load([pathpp 'trackmouseeye.mat'])
        validRFCIfix = true;
        gazedistthresh = 20;
    else
        validRFCIfix = false;
    end
    load([pathpp 'postprocessed.mat'])
    load([pathpp 'visresponses.mat'])


    % load([pathpp 'postprocessed_probeC.mat'])
    % load([pathpp 'visresponses_probeC.mat'])
    % psthtli = (-500:1000)';
    % tloi = psthtli>0 & psthtli<=1000;
    % tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,1:4:end,:), 1))';
    % tempRall = squeeze(mean(reshape(Rall.RFCI_presentations(:,neuoind),4,length(temptrialorder),numel(neuoind)),1))';
    % isequal(tempR, tempRall) % not exactly the same
    % figure;imagesc(tempRall-tempR)
    % figure;plot(tempR(:,30), tempRall(:,30),'o')

    Nneurons = length(meanFRall);
    temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
    tempR = squeeze(mean(reshape(Rall.RFCI_presentations,4,length(temptrialorder),Nneurons),1))';

    RFCIall = analyzeRFCI(tempR, temptrialorder, sponFRall);

    % RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning grating
    % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
    % durstim = vis.RFCI_presentations.stop_time-vis.RFCI_presentations.start_time;
    % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])
    % durstim = vis.RFCI_presentations.start_time(2:end)-vis.RFCI_presentations.stop_time(1:end-1);
    % disp([mean(durstim) median(durstim) min(durstim) max(durstim)])

    %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
    % 10000's: which type (classic 0 vs inverse 1), 1000's which ctrsizes,
    % 10-100's: which RFcenter, 1's: which direction

    % temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
    % tloi = psthtli>0 & psthtli<=1000;
    % tempR = squeeze(1000*mean(psth.RFCI_presentations(tloi,1:4:end,:), 1))';

    if exist([pathpp 'trackmouseeye.mat'], 'file')
        tempgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
        templikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
        temptrialsfixedgaze = tempgazedist<gazedistthresh & ~templikelyblink;

        % validRFCIfix = true;
        fixcrftrials = temptrialsfixedgaze & floor(temptrialorder/10000) == 0;% &
        if ~isequal( unique(floor(mod(temptrialorder(fixcrftrials), 1000) / 10)), (0:8)' )
            validRFCIfix = false;
        end
        fixirftrials = temptrialsfixedgaze & floor(temptrialorder/10000) == 1;% & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
        % if ~isequal( unique(floor(mod(temptrialorder(fixirftrials), 1000) / 10)), (0:8)' )
        %     validRFCIfix = false;
        % end
    end

    if validRFCIfix
        RFCIall_fixedgaze = analyzeRFCI(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze), sponFRall);
        fprintf('fixed gaze ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
            nnz(RFCIall_fixedgaze.RFindclassic==1 & RFCIall_fixedgaze.pRFclassic<0.05), ...
            nnz(RFCIall_fixedgaze.RFindclassic==1 & RFCIall_fixedgaze.Pkw_rfclassic<0.05), ...
            nnz(RFCIall_fixedgaze.RFindclassic==1 & RFCIall_fixedgaze.sigmc_rfclassic), ...
            nnz(RFCIall_fixedgaze.RFindclassic==1 & RFCIall_fixedgaze.RFsigexclclassic), ...
            nnz(RFCIall_fixedgaze.RFindclassic==1 & RFCIall_fixedgaze.RFexclsigclassic) )
    else
        RFCIall_fixedgaze = struct();
        disp('fixedgaze RFCI block skipped')
    end

    save([pathpp 'RFCIall_onetail.mat'], 'RFCIall', 'RFCIall_fixedgaze')
end

%%
%% report percent trials left in modified fixed gaze condition, where only the first 150 ms out of 250 ms of each orientation is analyzed
% in the RFCI block, each trial is 1s spinning grating, 0.25s at each orientation
% with gazedistthresh of 10pix (~4visual degrees), this leaves ~30 percent
% trials in most sessions
for ises = 1:Nsessions
    clearvars -except ises Nsessions nwbsessions datadir
    sesclk = tic;
    fprintf('\nSession %d %s\n', ises, nwbsessions{ises})

    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
    probes = {'A', 'B', 'C', 'D', 'E', 'F'};
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed.mat'])
    load([pathpp 'visresponses.mat'])
    if exist([pathpp 'trackmouseeye.mat'], 'file')
        load([pathpp 'trackmouseeye.mat'])
    else
        continue
    end
trackeyetli = -30:60;
    % 1/nanmedian(diff(TrackEyeTimestamps)) roughly 60 Hz frame rate
% eyecamframerate = 1/nanmedian(diff(TrackEyeTimestamps));
eyecamframerate = 60;
trackeyetl = trackeyetli/eyecamframerate;

% %% get trialpupildata
load([pathpp 'trackmouseeye.mat'])
trialpupildata = struct();
for b = 1:numel(visblocks)
    if contains(visblocks{b}, 'spontaneous')
        continue
    end
    [r,c]=find(cumsum(TrackEyeTimestamps-vis.(visblocks{b}).trialstart'>0,1)==1);
    if ~isequal(c, (1:numel(vis.(visblocks{b}).trialstart))' )
        error('missing some trials')
    end
    % figure; plot(TrackEyeTimestamps(r),vis.(visblocks{b}).trialstart, 'o')
    % max(abs(TrackEyeTimestamps(r)-vis.(visblocks{b}).trialstart)) % 0.0167s, i.e., 1/60s
    trackeyetrialinds = (r-1)+trackeyetli;
%     trackeyetl = trackeyetli/eyecamframerate;

    tempdata = squeeze(pupiltracking.data(1,:));
    trialpupildata.(visblocks{b}).x = tempdata(trackeyetrialinds);
    tempdata = squeeze(pupiltracking.data(2,:));
    trialpupildata.(visblocks{b}).y = tempdata(trackeyetrialinds);
end


% %% get frontRFCItrialsfix10pix
b = find(strcmp(visblocks, 'RFCI_presentations'));
[r,c]=find(cumsum(TrackEyeTimestamps-vis.(visblocks{b}).trialstart'>0,1)==1);
if ~isequal(c, (1:numel(vis.(visblocks{b}).trialstart))' )
    error('missing some trials')
end
trackeyetrialinds = (r-1)+trackeyetli;
trackeyepsth = distmodecom(trackeyetrialinds);
likelyblinkpsth = likelyblink(trackeyetrialinds);

% 0.25seconds at each orientation
winfront = 0.15; % sec,
RFCIfronttloi = (trackeyetl>=0 & trackeyetl<winfront) | (trackeyetl>=0.25 & trackeyetl<0.25+winfront) | ...
    (trackeyetl>=0.25*2 & trackeyetl<0.25*2+winfront) | (trackeyetl>=0.25*3 & trackeyetl<0.25*3+winfront);
trialmaxdistmodecom.RFCIfront10pix = max(trackeyepsth(:,RFCIfronttloi),[],2);
triallikelyblink.RFCIfront10pix = any(likelyblinkpsth(:,RFCIfronttloi), 2);

frontRFCIgazedist = trialmaxdistmodecom.RFCIfront10pix(1:4:end);
frontRFCIlikelyblink = triallikelyblink.RFCIfront10pix(1:4:end);
frontRFCItrialsfix10pix = frontRFCIgazedist<10 & ~frontRFCIlikelyblink;

RFCIgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
RFCIlikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
RFCItrialsfixedgaze = RFCIgazedist<20 & ~RFCIlikelyblink;

disp([mean(RFCItrialsfixedgaze), mean(frontRFCItrialsfix10pix)])
% win 0.125: 0.3722    0.2056
% win 0.15: 0.3722    0.2056

% %% calculate Rall_RFCIfront
Nneurons = length(sponFRall);
temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
tempR = squeeze(mean(reshape(Rall.RFCI_presentations,4,length(temptrialorder),Nneurons),1))';

psthtl = psthtli/1000;
RFCIfrontpsthtli = (psthtl>=0 & psthtl<winfront) | (psthtl>=0.25 & psthtl<0.25+winfront) | ...
    (psthtl>=0.25*2 & psthtl<0.25*2+winfront) | (psthtl>=0.25*3 & psthtl<0.25*3+winfront);
Rall_RFCIfront10pix = NaN(length(temptrialorder), Nneurons);
for iprobe = 1:numel(probes)
    load([pathpp 'postprocessed_probe' probes{iprobe} '.mat'])
Rall_RFCIfront10pix(:,neuoind) = squeeze(mean(1000*psth.RFCI_presentations(RFCIfrontpsthtli,1:4:end,:),1));
end

% % figure; plot(RFCIgazedist, frontRFCIgazedist,'o')
% % figure; plot(psthtli, squeeze(mean(psth.RFCI_presentations(:,1:4:end,:),[2,3])) )
% % set(gca, 'XTick', -500:50:1000, 'XGrid', 'on')
% load([pathpp 'postprocessed_probeC.mat'])
% tempR1 = squeeze(mean(1000*psth.RFCI_presentations(psthtli>0 & psthtli<=1000,1:4:end,:),1))';
% % figure; plot(tempR(neuoind,:),tempR1,'o') % almost the same, not exactly, when winfront=0.25
% %figure; plot(tempR, Rall_RFCIfront', 'o') % almost the same, not exactly, when winfront=0.25
% figure; plot(Rall_RFCIfront(:,neuoind)',tempR1,'o') % exactly the same when winfront=0.25
% RR125 = diag(corr(tempR', Rall_RFCIfront125, 'Type','Spearman'));
% RR150 = diag(corr(tempR', Rall_RFCIfront150, 'Type','Spearman'));
% disp([nanmean(RR125) prctile(RR125,[25 50 75])]) % 0.7855    0.7333    0.7913    0.8447
% disp([nanmean(RR150) prctile(RR150,[25 50 75])]) % 0.8395    0.7985    0.8461    0.8889
% figure; hold all; histogram(RR150); histogram(RR125)

% RESUME HERE
% validRFCIfix = true;
fixcrftrials = frontRFCItrialsfix10pix & floor(temptrialorder/10000) == 0;% &
if ~isequal( unique(floor(mod(temptrialorder(fixcrftrials), 1000) / 10)), (0:8)' )
    validRFCIfix = false;
end
fixirftrials = frontRFCItrialsfix10pix & floor(temptrialorder/10000) == 1;% & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
% if ~isequal( unique(floor(mod(temptrialorder(fixirftrials), 1000) / 10)), (0:8)' )
%     validRFCIfix = false;
% end

if validRFCIfix
    % NOTE input to analyzeRFCI tempR is transposed from Rall
    RFCIall_frontfix10pix = analyzeRFCI(Rall_RFCIfront10pix(frontRFCItrialsfix10pix,:)', temptrialorder(frontRFCItrialsfix10pix), sponFRall);
    fprintf('fixed gaze ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
        nnz(RFCIall_frontfix10pix.RFindclassic==1 & RFCIall_frontfix10pix.pRFclassic<0.05), ...
        nnz(RFCIall_frontfix10pix.RFindclassic==1 & RFCIall_frontfix10pix.Pkw_rfclassic<0.05), ...
        nnz(RFCIall_frontfix10pix.RFindclassic==1 & RFCIall_frontfix10pix.sigmc_rfclassic), ...
        nnz(RFCIall_frontfix10pix.RFindclassic==1 & RFCIall_frontfix10pix.RFsigexclclassic), ...
        nnz(RFCIall_frontfix10pix.RFindclassic==1 & RFCIall_frontfix10pix.RFexclsigclassic) )
else
    RFCIall_frontfix10pix = struct();
    disp('fixedgaze RFCI block skipped')
end


tempgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
templikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
temptrialsfixedgaze = tempgazedist<20 & ~templikelyblink;

% validRFCIfix = true;
nonfixcrftrials = ~temptrialsfixedgaze & floor(temptrialorder/10000) == 0;% &
if ~isequal( unique(floor(mod(temptrialorder(nonfixcrftrials), 1000) / 10)), (0:8)' )
    validRFCInonfix = false;
end
nonfixirftrials = ~temptrialsfixedgaze & floor(temptrialorder/10000) == 1;% & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
% if ~isequal( unique(floor(mod(temptrialorder(fixirftrials), 1000) / 10)), (0:8)' )
%     validRFCIfix = false;
% end

if validRFCInonfix
    RFCIall_nonfix20pix = analyzeRFCI(tempR(:,~temptrialsfixedgaze), temptrialorder(~temptrialsfixedgaze), sponFRall);
    fprintf('non-fixed gaze ctrCRF: pRFclassic* %d Pkw_rfclassic* %d sigmc_rfclassic %d RFsigexclclassic %d RFexclsigclassic %d\n', ...
        nnz(RFCIall_nonfix20pix.RFindclassic==1 & RFCIall_nonfix20pix.pRFclassic<0.05), ...
        nnz(RFCIall_nonfix20pix.RFindclassic==1 & RFCIall_nonfix20pix.Pkw_rfclassic<0.05), ...
        nnz(RFCIall_nonfix20pix.RFindclassic==1 & RFCIall_nonfix20pix.sigmc_rfclassic), ...
        nnz(RFCIall_nonfix20pix.RFindclassic==1 & RFCIall_nonfix20pix.RFsigexclclassic), ...
        nnz(RFCIall_nonfix20pix.RFindclassic==1 & RFCIall_nonfix20pix.RFexclsigclassic) )
else
    RFCIall_nonfix20pix = struct();
    disp('fixedgaze RFCI block skipped')
end

save([pathpp 'RFCIall_frontfix.mat'], 'trialpupildata', 'winfront', 'frontRFCItrialsfix10pix', ...
    'trialmaxdistmodecom', 'triallikelyblink', 'Rall_RFCIfront10pix', ...
    'RFCIall_frontfix10pix', 'RFCIall_nonfix20pix')
    
toc(sesclk)
end


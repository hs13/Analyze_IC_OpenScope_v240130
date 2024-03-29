% SHOULD GAZE BE DEFINED BASED ON EYETRACKING POSITION OR PUPILTRACKING POSITION?
% NOTE CURRENTLY I DEFINED GAZE BASED ON PUPIL POSITION
% Siegle et al Neuropixels platfor paper:
% Across 50 mice with processed eye-tracking videos, we used 
% the gaze_mapping module of the AllenSDK to translate pupil position into 
% screen coordinates (in units of degrees). On average, 95% of gaze locations 
% fell within 6.4 ± 2.1° of the mean, with a maximum of 13.6°.


% <4 vis deg (stricter criterion) for fixed gaze and replicate Fig1 results (R2C1.1)
% 
% eye position on different trial types (esp. IC vs LC) (R1C1)
% pupil area on IC vs LC vs RE trials (R1C2)
% 
% Perhaps show that receptive field position is not different when using all trials vs fixed-gaze trials
% Alternatively, show that exclusively center responsive neurons defined with all trials do not respond to grating patches in other RF positions …

datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/00248_v240130/';
nwbdir = dir([datadir 'postprocessed']);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);
ises = 7;

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'trackmouseeye.mat'])
load([pathpp 'postprocessed.mat'], 'vis')

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};


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
    
    modecom = [xbinctrs(r) ybinctrs(c)];
%     disp(modecom)
%     figure; plot(pupiltracking.data(1,:),pupiltracking.data(2,:),'.')
%     hold on; scatter(modecom(1), modecom(2), 100,'rx', 'linewidth', 2)
%     figure; histogram2(pupiltracking.data(1,:),pupiltracking.data(2,:),'displaystyle', 'tile')
%     hold on; scatter(modecom(1), modecom(2), 100,'rx', 'linewidth', 2)
    
    distmodecom = sqrt(sum((pupiltracking.data'-modecom).^2,2));
    
figure; plot(eyetracking.data(1,:),eyetracking.data(2,:),'.'); axis equal
axis([300 400 220 320])
figure; 
subplot(2,2,1)
histogram2(eyetracking.data(1,:),eyetracking.data(2,:),'displaystyle','tile');
hold on; scatter(modecom(1),modecom(2), 100,'rx', 'linewidth', 2)
axis([300 400 220 320])
title('eye position')
subplot(2,2,2)
histogram2(eyetracking.data(1,distmodecom<20),eyetracking.data(2,distmodecom<20),'displaystyle','tile');
hold on; scatter(modecom(1),modecom(2), 100,'rx', 'linewidth', 2)
axis([300 400 220 320])
title('eye position distmodecom<20')
subplot(2,2,3)
histogram2(pupiltracking.data(1,:),pupiltracking.data(2,:),'displaystyle','tile');
hold on; scatter(modecom(1),modecom(2), 100,'rx', 'linewidth', 2)
axis([300 400 220 320])
title('pupil position')
subplot(2,2,4)
histogram2(pupiltracking.data(1,distmodecom<20),pupiltracking.data(2,distmodecom<20),'displaystyle','tile');
hold on; scatter(modecom(1),modecom(2), 100,'rx', 'linewidth', 2)
axis([300 400 220 320])
title('pupil position distmodecom<20')

%% trials should start -0.5s before and end 1s after stim onset
% 1/nanmedian(diff(TrackEyeTimestamps)) roughly 60 Hz frame rate
eyecamframerate = 60;
trackeyetli = -30:60;
trialdistmodecom = struct();
likelyblinkpsth = struct();
trialpupildata = struct();
tic
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
    trackeyetl = trackeyetli/eyecamframerate;

    trackeyepsth = distmodecom(trackeyetrialinds);

    trialdistmodecom.(visblocks{b}).trackeyetli = trackeyetli;
    trialdistmodecom.(visblocks{b}).psthtrialinds = trackeyetrialinds;
    trialdistmodecom.(visblocks{b}).psth = trackeyepsth;

    likelyblinkpsth.(visblocks{b}) = likelyblink(trackeyetrialinds);
    
    tempdata = squeeze(pupiltracking.data(1,:));
    trialpupildata.(visblocks{b}).x = tempdata(trackeyetrialinds);
    tempdata = squeeze(pupiltracking.data(2,:));
    trialpupildata.(visblocks{b}).y = tempdata(trackeyetrialinds);
end

%{
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
%}

%% eye position on IC vs LC vs IRE trials
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];

whichblock = 'ICwcfg1_presentations';
tloi = trackeyetli>0 & trackeyetli<=0.4*eyecamframerate;

trialsoi = vis.(whichblock).trialorder
trialpupildata.(whichblock).x(:,tloi)
trialpupildata.(whichblock).y(:,tloi)

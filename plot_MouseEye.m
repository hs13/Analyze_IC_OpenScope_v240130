% <4 vis deg (stricter criterion) for fixed gaze and replicate Fig1 results (R2C1.1)
% 
% eye position on different trial types (esp. IC vs LC) (R1C1)
% pupil area on IC vs LC vs RE trials (R1C2)
% 
% Perhaps show that receptive field position is not different when using all trials vs fixed-gaze trials
% Alternatively, show that exclusively center responsive neurons defined with all trials do not respond to grating patches in other RF positions …

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'trackmouseeye.mat'])

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};

% %% trials should start -0.5s before and end 1s after stim onset
% 1/nanmedian(diff(TrackEyeTimestamps)) roughly 60 Hz frame rate
trackeyetli = -30:60;
trialdistmodecom = struct();
trialdistmodecom = struct();
trialdistmodecom = struct();
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
    eyecamframerate = 60;
    trackeyetl = trackeyetli/eyecamframerate;

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

end
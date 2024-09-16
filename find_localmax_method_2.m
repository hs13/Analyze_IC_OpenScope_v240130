% code snippet from SpectralEvents/spectralevents_find.m

%   2) Find spectral events by first thresholding
%      entire normalize TFR (over all frequencies), then finding local 
%      maxima. Discard those of lesser magnitude in each suprathreshold 
%      region, respectively, s.t. only the greatest local maximum in each 
%      region survives (when more than one local maxima in a region have 
%      the same greatest value, their respective event timing, freq. 
%      location, and boundaries at full-width half-max are calculated 
%      separately and averaged). This method does not allow for overlapping
%      events to occur in a given suprathreshold region and does not 
%      guarantee the presence of within-band, suprathreshold activity in 
%      any given trial will render an event.

% Inputs:
%   findMethod - integer value specifying which event-finding method to use
%       (1, 2, or 3). Note that the method specifies how much overlap
%       exists between events. Use 1 to replicate the method used in 
%       et al. eLife 2017.
%   eventBand - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
%       which above-threshold spectral power events are classified.
%   thrFOM - factors of median threshold; positive real number used to
%       threshold local maxima and classify events (see Shin et al. eLife 
%       2017 for discussion concerning this value).
%   tVec - time vector (s) over which the time-frequency response (TFR) is 
%       calcuated.
%   fVec - frequency vector (Hz) over which the time-frequency response 
%       (TFR) is calcuated.
%   TFR - time-frequency response (TFR) (frequency-by-time-trial) for a
%       single subject/session.
%   classLabels - numeric or logical 1-row array of trial classification 
%       labels; associates each trial of the given subject/session to an 
%       experimental condition/outcome/state (e.g., hit or miss, detect or 
%       non-detect, attend-to or attend away).

function [spectralEvents] = find_localmax_method_2(eventBand, thrFOM, tVec, fVec, TFR, classLabels)
% 2nd event-finding method: Find spectral events by first thresholding
% entire normalize TFR (over all frequencies), then finding local
% maxima. This method does not allow for overlapping events to occur in
% a given suprathreshold region and does not guarantee the presence of
% within-band, suprathreshold activity in any given trial will render
% an event.

% spectralEvents: 12 column matrix for storing local max event metrics: trial
% index, hit/miss, maxima frequency, lowerbound frequency, upperbound
% frequency, frequency span, maxima timing, event onset timing, event
% offset timing, event duration, maxima power, maxima/median power
spectralEvents = [];

% Retrieve local maxima in normalized TFR using imregionalmax,
% discard those of lesser (un-normalized) magnitude in each suprathreshold
% region, respectively, and characterize event boundaries (at half max)
for ti=1:numTrials
    TFR_ST = squeeze(TFR(:,:,ti))./medianpower; %Suprathreshold TFR: first isolate 2D TFR matrix and normalize
    TFR_ST(TFR_ST<thrFOM) = 0; %Set infrathreshold values to zero
    
    % Find all local maxima in suprathreshold TFR
    TFR_LM = TFR_ST.*imregionalmax(TFR_ST); %Threshold TFR at each respective local maximum
    numTotalPeaks = nnz(TFR_LM);
    
    % Escape this iteration when this trial contains no suprathreshold
    % local maxima
    if numTotalPeaks==0
        continue
    end
    
    % Find max peak in each respective suprathreshold region
    [~,regions,numReg,~] = bwboundaries(TFR_ST>=thrFOM); %Separate suprathreshold regions
    evPeakF = cell(1,numReg);
    evPeakT = cell(1,numReg);
    evPeakpower = nan(numReg,1);
    for reg_i=1:numReg
        region = zeros(size(TFR_ST)); %Initialize a blank image that will contain a single region
        region(regions==reg_i) = 1; %Set elements (pixels) in region to the value 1
        TFR_reg = TFR_LM.*region; %Regional local maxima
        [peakF_reg,peakT_reg] = find(TFR_reg); %Indices of regional local maxima
        peakpower_reg = TFR(find(TFR_reg)+(ti-1)*flength*tlength); %Power values at regional local maxima
        maxPeakpower = max(peakpower_reg);
        maxPeak_inds = find(peakpower_reg==maxPeakpower); %Indices of all instances where local maxima have the max peak power
        evPeakF{reg_i} = peakF_reg(maxPeak_inds); %Select TFR indices at max regional peak
        evPeakT{reg_i} = peakT_reg(maxPeak_inds); %Select TFR indices at max regional peak
        evPeakpower(reg_i) = maxPeakpower(1);
    end
    
    % Find local maxima lowerbound, upperbound, and full width at half max
    % for both frequency and time
    evBndsF = nan(numReg,3);
    evBndsT = nan(numReg,3);
    evPeakF_inds = nan(numReg,1);
    evPeakT_inds = nan(numReg,1);
    evPeakpower_norm = nan(numReg,1);
    for reg_i=1:numReg
        numRegPeaks = numel(evPeakF{reg_i});
        peakF = evPeakF{reg_i};
        peakT = evPeakT{reg_i};
        peakpower = evPeakpower(reg_i);
        
        Ffwhm = nan(numRegPeaks,3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
        Tfwhm = nan(numRegPeaks,3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
        peakpower_norm = nan(numRegPeaks,1); %Vector for storing the normalized power at each regional peak
        for lmi=1:numRegPeaks
            lmF_underthr = find(squeeze(TFR(:,peakT(lmi),ti) < peakpower/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
            if ~isempty(find(lmF_underthr < peakF(lmi), 1)) && ~isempty(find(lmF_underthr > peakF(lmi), 1))
                Ffwhm(lmi,1) = fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
                Ffwhm(lmi,2) = fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
                Ffwhm(lmi,3) = Ffwhm(lmi,2)-Ffwhm(lmi,1)+ min(diff(fVec));
            elseif isempty(find(lmF_underthr < peakF(lmi),1)) && ~isempty(find(lmF_underthr > peakF(lmi),1))
                Ffwhm(lmi,1) = fVec(1);
                Ffwhm(lmi,2) = fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
                Ffwhm(lmi,3) = 2*(Ffwhm(lmi,2)-fVec(peakF(lmi)))+ min(diff(fVec));
            elseif ~isempty(find(lmF_underthr < peakF(lmi),1)) && isempty(find(lmF_underthr > peakF(lmi),1))
                Ffwhm(lmi,1) = fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
                Ffwhm(lmi,2) = fVec(end);
                Ffwhm(lmi,3) = 2*(fVec(peakF(lmi))-Ffwhm(lmi,1))+ min(diff(fVec));
            else
                Ffwhm(lmi,1) = fVec(1);
                Ffwhm(lmi,2) = fVec(end);
                Ffwhm(lmi,3) = 2*(fVec(end)-fVec(1)+min(diff(fVec)));
            end
            
            lmT_underthr = find(squeeze(TFR(peakF(lmi),:,ti) < peakpower/2)); %Indices of TFR times of < half max power at the freq of a given local peak
            if ~isempty(find(lmT_underthr < peakT(lmi), 1)) && ~isempty(find(lmT_underthr > peakT(lmi), 1))
                Tfwhm(lmi,1) = tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
                Tfwhm(lmi,2) = tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
                Tfwhm(lmi,3) = Tfwhm(lmi,2)-Tfwhm(lmi,1)+ min(diff(tVec));
            elseif isempty(find(lmT_underthr < peakT(lmi),1)) && ~isempty(find(lmT_underthr > peakT(lmi),1))
                Tfwhm(lmi,1) = tVec(1);
                Tfwhm(lmi,2) = tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
                Tfwhm(lmi,3) = 2*(Tfwhm(lmi,2)-tVec(peakT(lmi)))+ min(diff(tVec));
            elseif ~isempty(find(lmT_underthr < peakT(lmi),1)) && isempty(find(lmT_underthr > peakT(lmi),1))
                Tfwhm(lmi,1) = tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
                Tfwhm(lmi,2) = tVec(end);
                Tfwhm(lmi,3) = 2*(tVec(peakT(lmi))-Tfwhm(lmi,1))+ min(diff(tVec));
            else
                Tfwhm(lmi,1) = tVec(1);
                Tfwhm(lmi,2) = tVec(end);
                Tfwhm(lmi,3) = 2*(tVec(end)-tVec(1)+min(diff(tVec)));
            end
            
            peakpower_norm(lmi) = TFR_ST(peakF(lmi),peakT(lmi));
        end
        
        evBndsF(reg_i,:) = mean(Ffwhm,1);
        evBndsT(reg_i,:) = mean(Tfwhm,1);
        evPeakF_inds(reg_i) = round(mean(peakF));
        evPeakT_inds(reg_i) = round(mean(peakT));
        evPeakpower_norm(reg_i) = mean(peakpower_norm);
    end
    
    % 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
    % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
    spectralEvents = [spectralEvents; ti*ones(size(evPeakpower)) classLabels(ti)*ones(size(evPeakpower))...
        fVec(evPeakF_inds)' evBndsF tVec(evPeakT_inds)' evBndsT evPeakpower evPeakpower_norm];
end

% Pick out maxima within the frequency band of interest
spectralEvents = spectralEvents((spectralEvents(:,3)>=eventBand(1) & spectralEvents(:,3)<=eventBand(2)),:); %Select local maxima
end
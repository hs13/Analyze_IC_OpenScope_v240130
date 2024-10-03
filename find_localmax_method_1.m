% code snippet from SpectralEvents/spectralevents_find.m

%   1) (Primary event detection method in Shin et al. eLife 2017): Find 
%      spectral events by first retrieving all local maxima in 
%      un-normalized TFR using imregionalmax, then selecting suprathreshold
%      peaks within the frequency band of interest. This method allows for 
%      multiple, overlapping events to occur in a given suprathreshold 
%      region and does not guarantee the presence of within-band, 
%      suprathreshold activity in any given trial will render an event.


function [spectralEvents] = find_localmax_method_1(eventBand, thrFOM, tVec, fVec, TFR, classLabels)
% 1st event-finding method (primary event detection method in Shin et
% al. eLife 2017): Find spectral events by first retrieving all local
% maxima in un-normalized TFR using imregionalmax, then selecting
% suprathreshold peaks within the frequency band of interest. This
% method allows for multiple, overlapping events to occur in a given
% suprathreshold region and does not guarantee the presence of
% within-band, suprathreshold activity in any given trial will render
% an event.

% Initialize general data parameters
eventBand_inds = fVec>=eventBand(1) & fVec<=eventBand(2); %Logical vector representing indices of freq vector within eventBand
if size(eventBand_inds,1)~=length(eventBand_inds)
    eventBand_inds = eventBand_inds'; %Transpose so that the dimensions correspond with the frequency-domain dimension of the TFR
end
flength = size(TFR,1); %Number of elements in discrete frequency spectrum
tlength = size(TFR,2); %Number of points in time
numTrials = size(TFR,3); %Number of trials
classes = unique(classLabels);
medianpower = median(reshape(TFR, size(TFR,1), size(TFR,2)*size(TFR,3)), 2); %Median power at each frequency across all trials
thr = thrFOM*medianpower; %Spectral event threshold for each frequency value

% Validate consistency of parameter dimensions
if flength~=length(fVec) || tlength~=length(tVec) || numTrials~=length(classLabels)
  error('Mismatch in input parameter dimensions!')
end

% spectralEvents: 12 column matrix for storing local max event metrics: trial
% index, hit/miss, maxima frequency, lowerbound frequency, upperbound
% frequency, frequency span, maxima timing, event onset timing, event
% offset timing, event duration, maxima power, maxima/median power
spectralEvents = [];

% Finds_localmax: stores peak frequecy at each local max (columns) for each
% trial (rows)
Finds_localmax = [];

% Retrieve all local maxima in TFR using imregionalmax
for ti=1:numTrials
    [peakF,peakT] = find(imregionalmax(squeeze(TFR(:,:,ti)))); %Indices of max local power
    peakpower = TFR(find(imregionalmax(squeeze(TFR(:,:,ti))))+(ti-1)*flength*tlength); %Power values at local maxima (vector; compiles across frequencies and time)
    
    % Find local maxima lowerbound, upperbound, and full width at half max
    % for both frequency and time
    Ffwhm = NaN(numel(peakpower),3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
    Tfwhm = NaN(numel(peakpower),3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
    for lmi=1:numel(peakpower)
        lmF_underthr = find(squeeze(TFR(:,peakT(lmi),ti) < peakpower(lmi)/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
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
        
        lmT_underthr = find(squeeze(TFR(peakF(lmi),:,ti) < peakpower(lmi)/2)); %Indices of TFR times of < half max power at the freq of a given local peak
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
    end
    
    % 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
    % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
    spectralEvents = [spectralEvents; ti*ones(size(peakF)) classLabels(ti)*ones(size(peakF)) fVec(peakF)' Ffwhm tVec(peakT)' Tfwhm peakpower peakpower./medianpower(peakF)];
    Finds_localmax = [Finds_localmax; peakF];
end

% Pick out maxima above threshold and within the frequency band of interest
spectralEvents = spectralEvents((spectralEvents(:,3)>=eventBand(1) & spectralEvents(:,3)<=eventBand(2) & spectralEvents(:,11)>=thr(Finds_localmax)),:); %Select local maxima
end

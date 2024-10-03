%% basic spectral event analyses
% are there clear "bands" of spectral event frequencies? (boundary between beta and gamma
% beta events are characterized by a prominent supragranuar sink -- do gamma events have a granular sink? and what about alpha events?

 
%% first, get psth of spectral events (vis stim aligned)


whichblock = 'spontaneous_presentations';
[durtrial, itrial] = max(vis.(whichblock).stop_time-vis.(whichblock).start_time);
lfpsnip = lfpvispsth.(whichblock){itrial};
lfpconv = convn(lfpsnip, reshape(kergauss,[],1), 'same');

% whichblock = 'spontaneous_presentations';
% Nelec = numel(lfpelecvec.location);
% figure; hold all
% plot(median(reshape(lfpvispsth.ICwcfg1_presentations,[],Nelec),1), 'k-', 'linewidth', 1)
% for itrial = 1:numel(lfpvispsth.(whichblock))
% lfpsnip = lfpvispsth.(whichblock){itrial};
% plot(median(lfpsnip,1))
% end

%   TFR - time-frequency response (TFR) (frequency-by-time-trial) for a
%       single subject/session.
thrFOM = 6;
eventBand = [15 29];
tVec = lfpvistimes.(whichblock){itrial};
if contains(whichblock, 'spontaneous')
    TFR = permute(tfrvispsth.(whichblock){itrial}, [2 1]);
    CSD = permute(csdvispsth.(whichblock){itrial}, [2 1]);
    classLabels = itrial;
else
    TFR = permute(tfrvispsth.(whichblock), [3 1 2]);
    CSD = permute(csdvispsth.(whichblock), [3 1 2]);
    classLabels = reshape(vis.(whichblock).trialorder,1,[]);
end

% spectralEvents is 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
% 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
evcols = {'trial_index', 'trial_class', 'maxima_frequency', ...
    'lowerbound_frequency', 'upperbound_frequency', 'frequency_span', ...
    'maxima_timing', 'event_onset_timing', 'event_offset_timing', ...
    'event_duration', 'maxima_power', 'maxima_normmedian_power'};
tic
[spectralEvents] = find_localmax_method_2(eventBand, thrFOM, tVec, fVec, ...
    TFR, classLabels);
toc % 200s for 312s period

%% in continuous data, calculate event triggered

addpath('C:\Users\USER\GitHub\helperfunctions')


pathpp = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/IllusionOpenScope/sub-620333/';
pathpp = 'S:\OpenScopeData\00248_v240130\postprocessed\sub-620333\';
load([pathpp 'LFP1000Hz_probeC.mat'])
load([pathpp 'LFP_psth_probeC.mat'])
load([pathpp 'LFP_CSD_probeC.mat'])
load([pathpp 'LFP_TFR_L23_probeC.mat'])

visblocks = fieldnames(vis);
lfpvispsthtrialinds = struct();
lfpvistimes = struct();
for b = 1:numel(visblocks)
    if contains(visblocks{b}, 'spontaneous')
        lfpvispsthtrialinds.(visblocks{b}) = cell(numel(vis.(visblocks{b}).start_time),1);
        for itrial = 1:numel(vis.(visblocks{b}).start_time)
            trsinds = find( lfptimeresamp>=vis.(visblocks{b}).start_time(itrial) ...
                & lfptimeresamp<=vis.(visblocks{b}).stop_time(itrial) );
            lfpvispsthtrialinds.(visblocks{b}){itrial} = trsinds;
            lfpvistimes.(visblocks{b}){itrial} = lfptimeresamp(trsinds);
        end
        continue
    end

    % actual start time could be up to 1ms after 0 index
    %psthtrialinds = floor(vis.(visblocks{b}).trialstart'/Tres)+1 + psthtli;
    Ntrials = numel(vis.(visblocks{b}).trialstart);
    lfppsthtrialinds = zeros(length(psthtli), Ntrials);
    for itrial = 1:Ntrials
        t0rsind = find(lfptimeresamp<=vis.(visblocks{b}).trialstart(itrial),1,'last');
        lfppsthtrialinds(:,itrial) = t0rsind + psthtli;
    end
    lfpvispsthtrialinds.(visblocks{b}) = lfppsthtrialinds;
    lfpvistimes.(visblocks{b}) = lfptimeresamp(lfppsthtrialinds);

    clear lfppsthtrialinds
end

Ntrials = numel(opto.optostarttime);
lfpoptopsthtrialinds = zeros(length(optopsthtli), Ntrials);
lfpoptotimes = zeros(length(optopsthtli), Ntrials);
for itrial = 1:Ntrials
    t0rsind = find(lfptimeresamp<=opto.optostarttime(itrial),1,'last');
    lfpoptopsthtrialinds(:,itrial) = t0rsind + optopsthtli;
end
lfpoptotimes = lfptimeresamp(lfpoptopsthtrialinds);

%%
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
eventBand = [7 100];
%eventBand = [15 29];
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
toc % 200s for 312s period, 406s for 1-100Hz

tic
[spectralEvents1] = find_localmax_method_1(eventBand, thrFOM, tVec, fVec, ...
    TFR, classLabels);
toc % 126s for 1-100Hz

figure; 
subplot(2,1,1)
histogram(spectralEvents1(:,strcmp(evcols,'maxima_frequency')), 0:99)
subplot(2,1,2)
histogram(spectralEvents(:,strcmp(evcols,'maxima_frequency')), 0:99)

figure; 
scatter(spectralEvents1(:,strcmp(evcols,'maxima_frequency')), ...
    spectralEvents1(:,strcmp(evcols,'maxima_normmedian_power')), 'o')

% tempx = spectralEvents(:, strcmp(evcols, 'maxima_frequency'));
% tempy = spectralEvents(:, strcmp(evcols, 'maxima_power'));
% figure; scatter(tempx, tempy, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.3)

figure; hold all
imagesc(tVec, fVec, TFR); 
tempx = spectralEvents(:, strcmp(evcols, 'maxima_timing'));
tempy = spectralEvents(:, strcmp(evcols, 'maxima_frequency'));
scatter(tempx,tempy, 'w*', 'LineWidth',1)
colorbar
colormap jet
t0 = tVec(randi(length(tVec)-10000, 1));
xlim(t0+[0 5]); ylim([7 80]); 
caxis(prctile(TFR(:), [2.5 97.5])); colorbar

figure; histogram(spectralEvents(:, strcmp(evcols, 'event_duration')))

%% top 50 events CSD
Nevents = size(spectralEvents,1);
[sv,si]= sort(spectralEvents(:, strcmp(evcols, 'maxima_power')),'descend');
topeventinds = si(1:50);
Ttopevents = spectralEvents(topeventinds, strcmp(evcols, 'maxima_timing'));

[indevent, tindevent] = find(cumsum(tVec>=spectralEvents(:, strcmp(evcols, 'maxima_timing')), 2)==1 );
if ~isequal(tVec(tindevent)', spectralEvents(:, strcmp(evcols, 'maxima_timing')) )
    error('check tindevent')
end

wvdur = round(1000/min(eventBand));
tevrange = -wvdur:wvdur;
evtroughtime = zeros(Nevents,1);
evtroughtind = zeros(Nevents,1);
evpeaktime = zeros(Nevents,1);
evpeaktind = zeros(Nevents,1);
for e = 1:Nevents
    temptinds = tindevent(e) + tevrange;
    [mv,mi]=min( lfpconv(temptinds,elecL23) );
    evtroughtind(e) = temptinds(mi);
    evtroughtime(e) = tVec(temptinds(mi));

    [mv,mi]=max( lfpconv(temptinds,elecL23) );
    evpeaktind(e) = temptinds(mi);
    evpeaktime(e) = tVec(temptinds(mi));
end
%figure; plot(tindevent, evtroughind,'.')

trange = -250:250;
yl = [7 80];
figure
for ii = 1:15
    tempevind = si(ii);
    temptinds = tindevent(tempevind) + trange;
subplot(3,5,ii)
hold all
imagesc(tVec(temptinds), fVec, TFR(:,temptinds))
tempvec = lfpconv(temptinds,elecL23)';
tempvec = tempvec-min(tempvec);
tempvecnorm = tempvec*0.9*range(yl)/range(tempvec)+yl(1)+0.05*range(yl);
plot(tVec(temptinds), tempvecnorm, 'w-', 'linewidth',1)
temptroughind = find(trange==0)+evtroughtind(tempevind)-tindevent(tempevind);
scatter(evtroughtime(tempevind), tempvecnorm(temptroughind), 30,'wv', 'markerfacecolor','m')
temppeakind = find(trange==0)+evpeaktind(tempevind)-tindevent(tempevind);
scatter(evpeaktime(tempevind), tempvecnorm(temppeakind), 30,'w^', 'markerfacecolor','m')
ylim(yl)
colorbar
colormap jet
end

% get CSD 200 ms before and after
% (Neletrodes-2) * time * Nevents
eventCSD = NaN(length(csdelectinds), length(trange), Nevents);
eventLFP = NaN(Nelec, length(trange), Nevents);
for e = 1:Nevents
    % eventCSD(:,:,e) = CSD(:, evtroughtind(e)+trange );
    % eventLFP(:,:,e) = lfpconv(evtroughtind(e)+trange, :)';

    eventCSD(:,:,e) = CSD(:, tindevent(e)+trange );
    eventLFP(:,:,e) = lfpconv(tindevent(e)+trange, :)';
end

ctxelec = contains(lfpelecvec.location, 'VIS');
ctxelectop = find(ctxelec, 1, 'last');
ctxelecbottom = find(ctxelec, 1, 'first');
if ~isequal(ctxelecbottom:ctxelectop, find(ctxelec)')
    warning('cortex electrodes are not consecutive -- check')
end

% top 50 beta events average CSD and LFP
figure('Position',[100 100 300 300])
hold all
yl = [ctxelecbottom ctxelectop]+.5;
imagesc(trange, csdelectinds, squeeze(mean(eventCSD(:,:,topeventinds),3)))
tempvec = squeeze( mean(eventLFP(elecL23,:,topeventinds), 3) );
tempvec = tempvec-min(tempvec);
tempvecnorm = tempvec*0.9*range(yl)/range(tempvec)+yl(1)+0.05*range(yl);
plot(trange, tempvecnorm, 'k-', 'linewidth',2)
scatter(0,elecL23,20,'w*', 'linewidth', 2)
xlim([-100 100])
ylim(yl)
set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpelecvec.location(csdelectinds), 'YDir', 'normal')
colorbar
colormap jet
caxis(0.04*[-1 1])

%% average TFR of top 50 elecL23 troughs 
Ntroughs = 50;
[PKS,LOCS] = findpeaks(-CSD(csdelectinds==elecL23,:));
[sv,si]=sort(PKS,'descend');

csdtroughtind = LOCS(si(1:Ntroughs));
csdtrough = -sv(1:Ntroughs);
csdtroughtime = tVec(csdtroughtind);

trange = -250:250;
csdtroughTFR = NaN(length(fVec), length(trange), Ntroughs);
csdtroughLFP = NaN(Nelec, length(trange), Ntroughs);
for e = 1:Ntroughs
    csdtroughTFR(:,:,e) = TFR(:, csdtroughtind(e)+trange );
    csdtroughLFP(:,:,e) = lfpconv(csdtroughtind(e)+trange, :)';
end

figure
imagesc(trange, fVec, squeeze(mean(csdtroughTFR,3)))
set(gca, 'YDir', 'normal')
ylim([8 80])
caxis([0 2*10^-8])
colormap jet
colorbar

figure; hold all
shadedErrorBar(fVec, squeeze(mean(TFR,2)), squeeze(std(TFR,0,2)), {'k-', 'linewidth', 2},1 )
temptfr = reshape(csdtroughTFR,length(fVec),[]);
shadedErrorBar(fVec, squeeze(mean(temptfr,2)), squeeze(std(temptfr,0,2)), {'r-', 'linewidth', 2},1 )
xlim([7 80])

temptfr = reshape(csdtroughTFR,length(fVec),[]);
figure; plot(fVec, squeeze(mean(temptfr,2))./squeeze(mean(TFR,2)) )

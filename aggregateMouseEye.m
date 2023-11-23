datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

%%
whichblock = 'ICwcfg1_presentations';

    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'trackmouseeye.mat'])

% find timestamps that correspond to each trial on IC blocks

trackeyetli = trialdistmodecom.(whichblock).trackeyetli;
trackeyetrialinds = trialdistmodecom.(whichblock).psthtrialinds;
% trackeyepsth = trialdistmodecom.(whichblock).psth;

% sanity check
%{
figure; hold all
plot(vis.(whichblock).trialstart, TrackEyeTimestamps(trackeyetrialinds(:, trackeyetli==0)), 'o')
plot(vis.(whichblock).trialstart, vis.(whichblock).trialstart, '-')

% roughly 60 Hz frame rate, max difference should be around that
figure
histogram(TrackEyeTimestamps(trackeyetrialinds(:, trackeyetli==0)) - vis.(whichblock).trialstart)
% vis trial start is always slightly later

max(abs(TrackEyeTimestamps(trackeyetrialinds(:, trackeyetli==0)) - vis.(whichblock).trialstart))

figure; histogram(diff(TrackEyeTimestamps))
%}

ifi = nanmean(diff(TrackEyeTimestamps)); % inter-frame interval

gazepsthfields = {'x', 'y'};
gazepsth = struct();
for f= 1:numel(gazepsthfields)
tempdata = pupiltracking.data(f,:);
gazepsth.(gazepsthfields{f}) = tempdata(trackeyetrialinds);
end

ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
trialtypes = unique(vis.(whichblock).trialorder);
ttoi = trialtypes(ismember(ICtrialtypes, [106 107 110 111]));

% eye position
figure
for typi = 1:numel(ttoi)
    trialsoi = vis.(whichblock).trialorder==ttoi(typi);
    tloi = trackeyetli>0 & trackeyetli<=round(0.4/ifi);
    subplot(2,2,typi)
    scatter(gazepsth.x(trialsoi,tloi), gazepsth.y(trialsoi,tloi), 5,'ko', 'filled', 'markeredgecolor', 'none', 'markerfacealpha', 0.1)
    axis([330 390 160 220 ])
end

% eye trajectory
tloi = trackeyetli>0 & trackeyetli<=round(0.4/ifi);
x0 = gazepsth.x(:,tloi)-gazepsth.x(:,trackeyetli==0);
y0 = gazepsth.y(:,tloi)-gazepsth.y(:,trackeyetli==0);
figure
for typi = 1:numel(ttoi)
    trialsoi = vis.(whichblock).trialorder==ttoi(typi);
    subplot(2,2,typi)
    hold all
    % plot(x0(trialsoi,:), y0(trialsoi,:))
    plot(nanmean(x0(trialsoi,:),1), nanmean(y0(trialsoi,:),1), 'k-', 'Linewidth', 2)
    scatter(nanmean(x0(trialsoi,:),1), nanmean(y0(trialsoi,:),1), 50, 1:nnz(tloi), 'filled')
    axis([-0.3 0.3 -0.3 0.3])
end

% eye position at the end of trial
figure
for typi = 1:numel(ttoi)
    trialsoi = vis.(whichblock).trialorder==ttoi(typi);
tloi = trackeyetli>0 & trackeyetli<=round(0.4/ifi);
    ttrialend = find(tloi,1,'last');
    subplot(2,2,typi)
    scatter(gazepsth.x(trialsoi,ttrialend), gazepsth.y(trialsoi,ttrialend), 5,'ko', 'filled', 'markeredgecolor', 'none', 'markerfacealpha', 0.1)
    axis([330 390 160 220 ])
end

% fit 2d gaussian to eye position




% compare pupil dilation on different trials
pupilareapsth = pupiltracking.area(trackeyetrialinds);
tloi = trackeyetli>0 & trackeyetli<=round(0.4/ifi);
trialendpupilarea = pupilareapsth(:,find(tloi,1,'last') );
trialmaxpupilarea = max(pupilareapsth(:,tloi), [],2);
trialmeanpupilarea = mean(pupilareapsth(:,tloi), 2);

trialsoi = ismember(vis.(whichblock).trialorder, ttoi);
[p,tbl,stats]=kruskalwallis(trialendpupilarea(trialsoi), vis.(whichblock).trialorder(trialsoi));
disp(p)
figure; multcompare(stats)

[p,tbl,stats]=kruskalwallis(trialmaxpupilarea(trialsoi), vis.(whichblock).trialorder(trialsoi));
disp(p)

[p,tbl,stats]=kruskalwallis(trialmeanpupilarea(trialsoi), vis.(whichblock).trialorder(trialsoi));
disp(p)

ttoi = trialtypes(ismember(ICtrialtypes, [106 111 1105 1109]));
trialsoi = ismember(vis.(whichblock).trialorder, ttoi);
[p,tbl,stats]=kruskalwallis(trialendpupilarea(trialsoi), vis.(whichblock).trialorder(trialsoi));
disp(p)

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
figure;
hold all
for b = 1:numel(visblocks)
    temptrackeyetli = trialdistmodecom.(visblocks{b}).trackeyetli;
    temptrackeyetrialinds = trialdistmodecom.(visblocks{b}).psthtrialinds;
    temppupilareapsth = pupiltracking.area(temptrackeyetrialinds);
    plot(ifi*temptrackeyetli, nanmean(temppupilareapsth,1)) % highest at the beginning of the trial, constricts during visual presentation due to luminance
end
set(gca, 'XTick', -0.4:0.4:1.2, 'XGrid', 'on')
legend(visblocks)

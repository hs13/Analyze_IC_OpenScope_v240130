ises = 12;
iprobe = 3;

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
% visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
%     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
% probes = {'C', 'D'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];

load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
% 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'

spklatencybins = [1 5 10 25]; % bin size for smoothing histogram
spklatencyT0s = [0 5 10]; % minimum spike latency
load(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{iprobe}))

open spklatency
%%
b = 4;
tt=3; ibin = 3;
whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));

ci = 45;
ii = 5;

trialsoi = vis.(visblocks{b}).trialorder==vistrialtypes.(visblocks{b})(ii);
vistrialrep.(visblocks{b})(ii) = nnz(trialsoi);

temptloi = psthtli>spklatencyT0s(tt) & psthtli<=250;
temptli = psthtli(temptloi);

% temppsth is Ntimepoints X Ntrials
temppsth = squeeze(psth.(visblocks{b})(temptloi,trialsoi,ci));
if nnz(temppsth)==0
    warning('%d spikes for this trial type', nnz(temppsth))
end

temppsthpre = squeeze(psth.(visblocks{b})(psthtli>-250 & psthtli<=0,trialsoi,ci));
[rpre,cpre]=find(temppsthpre & cumsum(temppsthpre,1)==1);
hcpre = histcounts(rpre, 0.5:1:nnz(temptloi)+0.5);

% find the first spike on each trial
[r,c]=find(temppsth & cumsum(temppsth,1)==1);
% spklatency.bin1ms.(visblocks{b})(ci,ii) = mode(r);

% % alternative way to find the first spike on each trial
% %[~,mi] = max(temppsth,[],1);
% %isequal(mi(c)', r)

hc = histcounts(r, 0.5:1:nnz(temptloi)+0.5);

figure; plot(temptli, hc)

figure; 
for ibin = 1:numel(spklatencybins)
hcsm = smooth( hc, spklatencybins(ibin) );
[mv,mi] = max(hcsm);
if spklatencybins(ibin)==1
    if mi ~= mode(r)
        error('check spike latency calculation for 1ms bin')
    end
end

hbw = 2*ceil(spklatencybins(ibin)/2)-1;
hbe = -0.5:hbw:nnz(temptloi)-hbw+0.5;
hcbin = NaN(hbw, length(hbe)-1);
for offset = 1:hbw
hcbin(offset,:) = histcounts(r, 'BinEdges', offset+hbe);
end

subplot(2,2,ibin)
hold all
plot(temptli, hcsm, 'b-')
plot(temptli(1)+(hbw-1)/2+(0:numel(hcbin)-1), hcbin(:)/hbw, 'r--')
title(sprintf('bin size %dms', spklatencybins(ibin)))
end

x=20.5

% NOTE, SMOOTHING IS THE SAME AS HISTOGRAM WITH CORRESPONDING BIN SIZE ONLY
% WHEN BIN SIZE IS AN ODD NUMBER (even numbers are reduced by 1)

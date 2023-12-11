datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

spklatencybins = [1 5 10 25]; % bin size for smoothing histogram
spklatencyT0s = [0 5 10]; % minimum spike latency

neulocagg = cell(numel(probes),Nsessions);
spklatencyagg = struct();
spklatencyprobagg = struct();
spklatencyadjagg = struct();
spklatencyadjprobagg = struct();
for iprobe = 1:numel(probes)
    for b = 1:numel(visblocks)
        for tt = 1:numel(spklatencyT0s)
            for ibin = 1:numel(spklatencybins)
                whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
                spklatencyagg.(whichbin).(visblocks{b}) = cell(numel(probes), Nsessions);
                spklatencyprobagg.(whichbin).(visblocks{b}) = cell(numel(probes), Nsessions);
                spklatencyadjagg.(whichbin).(visblocks{b}) = cell(numel(probes), Nsessions);
                spklatencyadjprobagg.(whichbin).(visblocks{b}) = cell(numel(probes), Nsessions);
            end
        end
    end
end

for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'

    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);

    for iprobe = 1:numel(probes)
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}), 'neuoind')
        load(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{iprobe}))
        

        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));

        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        neulocagg{iprobe, ises} = neuloc;

        for b = 1:numel(visblocks)
            for tt = 1:numel(spklatencyT0s)
                for ibin = 1:numel(spklatencybins)
                    whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
                    spklatency.(whichbin).(visblocks{b})(spklatencyprob.(whichbin).(visblocks{b})==0)= NaN;
                    spklatencyadj.(whichbin).(visblocks{b})(spklatencyadjprob.(whichbin).(visblocks{b})==0)= NaN;

                    spklatencyagg.(whichbin).(visblocks{b}){iprobe, ises} = ...
                        spklatency.(whichbin).(visblocks{b}) ;
                    spklatencyprobagg.(whichbin).(visblocks{b}){iprobe, ises} = ...
                        spklatencyprob.(whichbin).(visblocks{b}) ;
                    spklatencyadjagg.(whichbin).(visblocks{b}){iprobe, ises} = ...
                        spklatencyadj.(whichbin).(visblocks{b}) ;
                    spklatencyadjprobagg.(whichbin).(visblocks{b}){iprobe, ises} = ...
                        spklatencyadjprob.(whichbin).(visblocks{b}) ;
                end
            end
        end

    end
end

%%
neulocall = cat(1,neulocagg{:});
neuV1 = contains(neulocall, 'VISp') & ~contains(neulocall, 'VISpm');
neuLM = contains(neulocall, 'VISl');

ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
load('S:\OpenScopeData\00248_v230821\postprocessed\openscope_popavg_all.mat')

fs=12;

%% determine spklatencyprob threshold: 0.01
typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
for tt = 1:numel(spklatencyT0s)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
        subplot(numel(spklatencyT0s),numel(spklatencybins),numel(spklatencybins)*(tt-1)+ibin)

        templatencyprob = cat(1,spklatencyadjprobagg.(whichbin).(whichvisblock){:});
        histogram(templatencyprob(neuV1, typi) )%, 'normalization', 'cdf')
        title(whichbin)
    end
end

%% on IC and REl trials, compare V1 IC-encoder vs V1 inducer-encoder latency
bw=10;
adj = false;
if adj
    tempspklatencyprobagg = spklatencyadjprobagg;
    tempspklatencyagg = spklatencyadjagg;
    tempspklatencydesc = 'adjusted spike latency';
else
    tempspklatencyprobagg = spklatencyprobagg;
    tempspklatencyagg = spklatencyagg;
    tempspklatencydesc = 'spike latency';
end

figure;
annotation('textbox', [0.05 0.9 0.9 0.1], 'string', [tempspklatencydesc ' on I_C trials: V1 IC-encoder (red) vs V1 segment responders (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.9 0.9 0.1], 'string', [tempspklatencydesc ' on T_R_E trials: V1 IC-encoder (red) vs V1 segment responders (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

%%
figure;
annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [tempspklatencydesc ' on IC trials: LM IC-encoder (red) vs V1 IC-encoder (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuLM & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuLM & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [tempspklatencydesc ' on T_R_E trials: LM IC-encoder (red) vs V1 IC-encoder (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuLM & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuLM & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

%% LM IC-responsive (red) vs V1 IC-responsive (blue)
figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' on IC trials: LM IC-responsive (red) vs V1 IC-responsive (blue)'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuLM & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuLM & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' on I_R_E trials: LM IC-responsive (red) vs V1 IC-responsive (blue)'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuLM & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuLM & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' on T_R_E trials: LM IC-responsive (red) vs V1 IC-responsive (blue)'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuLM & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuLM & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

%% I_C (red) vs T_R_E (blue) trials
figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' V1 IC responsive neurons: I_C (red) vs T_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' pacman-encoder: I_C (red) vs T_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indenc13==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indenc24==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indenc13==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indenc24==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' segment responder: I_C (red) vs T_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' IC-encoder: I_C (red) vs T_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==1105), ...
                    templatency(neuoind2, ICtrialtypes==1109));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

%% I_C (red) vs I_R_E (blue) trials
figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' V1 IC responsive neurons: I_C (red) vs I_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICresp1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICresp2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' segment responder: I_C (red) vs I_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).indin1==1 | ICsigall.(whichvisblock).indin3==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).indin2==1 | ICsigall.(whichvisblock).indin4==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

figure;
annotation('textbox', [0.05 0.91 0.9 0.1], 'string', [tempspklatencydesc ' IC-encoder: I_C (red) vs I_R_E (blue) trials'], 'FontSize', fs, 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    plus0msbin = sprintf('plus0ms_bin%dms', spklatencybins(ibin));
    plus0mslatency = cat(1,tempspklatencyagg.(plus0msbin).(whichvisblock){:});
    plus5msbin = sprintf('plus10ms_bin%dms', spklatencybins(ibin));
    plus5mslatency = cat(1,tempspklatencyagg.(plus5msbin).(whichvisblock){:});
    templatency = plus5mslatency;
    templatency(plus0mslatency~=plus5mslatency)=NaN;
    templatency(templatency==250) = NaN;

    subplot(2,2,ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==506), ...
                    templatency(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = neuV1 & (ICsigall.(whichvisblock).ICencoder1==1);
                neuoind2 = neuV1 & (ICsigall.(whichvisblock).ICencoder2==1);
                tempspklatency = cat(1, templatency(neuoind1, ICtrialtypes==106), ...
                    templatency(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
    end
    set(gca, 'FontSize', fs)
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = signrank(tempspklatency1, tempspklatency2);
    title(sprintf('bin%dms p=%.4f', spklatencybins(ibin), p))
end

%% on IC and REl trials, compare ctrRF vs non-ctrRF latency
adj = false;
if adj
    tempspklatencyprobagg = spklatencyadjprobagg;
    tempspklatencyagg = spklatencyadjagg;
    tempspklatencydesc = 'adjusted spike lstency';
else
    tempspklatencyprobagg = spklatencyprobagg;
    tempspklatencyagg = spklatencyagg;
    tempspklatencydesc = 'spike lstency';
end

iprobe = find(strcmp(probes, 'C'));
typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
for tt = 1:numel(spklatencyT0s)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
        templatency = cat(1,tempspklatencyagg.(whichbin).(whichvisblock){:});
        subplot(numel(spklatencyT0s),numel(spklatencybins),numel(spklatencybins)*(tt-1)+ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind = neuV1 & (RFCIall.RFindclassic~=1 & RFCIall.pRFclassic<0.05);
            case 2
                neuoind = neuV1 & (RFCIall.RFindclassic==1 & RFCIall.pRFclassic<0.05);
        end
        %     neuoind = neuoind(tempspklatencyprobagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, typi)>spklatencyprobthr);
        histogram(templatency(neuoind, typi), 'binwidth', 5, 'normalization', 'probability')
    end
    title(whichbin)
end
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' on IC trials: center-RF (red) vs non-center-RF (blue)'], 'edgecolor', 'none')
for tt = 1:numel(spklatencyT0s)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
        templatency = cat(1,tempspklatencyagg.(whichbin).(whichvisblock){:});
        subplot(numel(spklatencyT0s),numel(spklatencybins),numel(spklatencybins)*(tt-1)+ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind = neuV1 & (RFCIall.RFindclassic~=1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==106), ...
                    templatency(neuoind, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind = neuV1 & (RFCIall.RFindclassic==1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==106), ...
                    templatency(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' non-center-RF neurons: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for tt = 1:numel(spklatencyT0s)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
        templatency = cat(1,tempspklatencyagg.(whichbin).(whichvisblock){:});
        subplot(numel(spklatencyT0s),numel(spklatencybins),numel(spklatencybins)*(tt-1)+ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind = neuV1 & (RFCIall.RFindclassic~=1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==506), ...
                    templatency(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind = neuV1 & (RFCIall.RFindclassic~=1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==106), ...
                    templatency(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
    end
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' center-RF neurons: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for tt = 1:numel(spklatencyT0s)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
        templatency = cat(1,tempspklatencyagg.(whichbin).(whichvisblock){:});
        subplot(numel(spklatencyT0s),numel(spklatencybins),numel(spklatencybins)*(tt-1)+ibin)
    hold all
    for ii = 1:2
        switch ii
            case 1
                neuoind = neuV1 & (RFCIall.RFindclassic==1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==506), ...
                    templatency(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind = neuV1 & (RFCIall.RFindclassic==1 & RFCIall.pRFclassic<0.05);
                tempspklatency = cat(1, templatency(neuoind, ICtrialtypes==106), ...
                    templatency(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
        histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
    end
end
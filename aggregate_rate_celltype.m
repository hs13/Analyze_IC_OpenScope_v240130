addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath('C:\Users\USER\GitHub\helperfunctions')
datadir = 'S:\OpenScopeData\00248_v240130\';

nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions datadir

    %% spike waveform

    % sesid = 'sub-620333';
    sesid = nwbsessions{ises};
    pathpp = [datadir 'postprocessed\' sesid '\'];

    load([pathpp 'info_units.mat'])
    load([pathpp 'qc_units.mat'])
    load([pathpp 'psth_spontaneous.mat'])

    neuV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
    neuRS = unit_wfdur>0.4;
    neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
    % Siegle et al's single unit filter criteria: isi_violations < 0.5 & amplitude_cutoff < 0.1 & presence_ratio > 0.9
    neufiltxtra = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9);

    nwbfiles = cat(1, dir([datadir sesid filesep '*.nwb']), dir([datadir sesid filesep '*' filesep '*.nwb']));
    % take filename  with shortest length or filename that does not contain probe
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    nwb = nwbRead(nwbspikefile); %, 'ignorecache');

    unit_waveform_mean = nwb.units.waveform_mean.data.load();

    unit_wfhalfwidth = nwb.units.vectordata.get('waveform_halfwidth').data.load();
    unit_waveform_mean_index = double( nwb.units.vectordata.get('waveform_mean_index').data.load() );

    Nelecperprobe = size(unit_waveform_mean,2)/length(unit_wfdur);
    if Nelecperprobe~=384
        error('electrodes per probe vary by session?')
    end

    % note that unit_peakch uses 0-based indexing, must add 1 to convert to matlab indexing
    unit_wfmeanind = Nelecperprobe*(0:length(unit_wfdur)-1)' + mod(double(unit_peakch+1), 1000);
    if isequal(unit_waveform_mean_index, unit_wfmeanind)
        disp('waveform_mean_index was set correctly')
    elseif isequal(unit_waveform_mean_index, Nelecperprobe*(1:length(unit_wfdur))')
        disp('waveform_mean_index was set incorrectly')
    else
        error('check waveform_mean_index')
    end

    unit_wfmean = unit_waveform_mean(:,unit_wfmeanind)';

    [troughV, troughT]=min(unit_wfmean,[],2);
    unit_wfmeannorm = unit_wfmean./(-troughV);
    % timeline aligned to trough
    unit_wftl = -troughT + (1:size(unit_wfmean,2));

    save([pathpp 'units_waveform.mat'], 'neuallloc', 'unit_peakch', 'unit_wfdur', ...
        'unit_isi_violations', 'unit_amplitude_cutoff', 'unit_presence_ratio', 'unit_amplitude', 'unit_wfhalfwidth', ...
        'unit_waveform_mean', 'unit_wfmeanind', 'unit_wfmean', 'unit_wfmeannorm', 'unit_wftl')


    %% plot spike waveforms
    %{
figure; plot(unit_wfdur, unit_wfhalfwidth, '.')
figure;histogram(unit_wfdur, 'binwidth', 0.05)
figure;histogram(unit_wfhalfwidth, 'binwidth', 0.02)

% spike waveform: RS vs FS
xl = [1 size(unit_wfmean,2)];
figure; 
subplot(2,3,1)
hold all
plot(unit_wfmean(neuV1, :)', 'k-')
plot(unit_wfmean(neuV1 & neufilt, :)', 'b-')
plot(unit_wfmean(neuV1 & neufiltxtra, :)', 'c-')
xlim(xl)
title(sprintf('V1 units n=%d', nnz(neuV1)))
subplot(2,3,2)
hold all
plot(unit_wfmean(neuV1 & neufilt, :)', 'b-')
xlim(xl)
title(sprintf('V1 soft-filtered units n=%d', nnz(neuV1 & neufilt)))
subplot(2,3,3)
hold all
plot(unit_wfmean(neuV1 & neufiltxtra, :)', 'c-')
xlim(xl)
title(sprintf('V1 Siegle-filtered units n=%d', nnz(neuV1 & neufiltxtra)))
subplot(2,3,4)
hold all
plot(unit_wfmean(neuV1 & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wfmean(neuV1 & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 RS+FS units n=%d + %d', nnz(neuV1 & neuRS), nnz(neuV1 & ~neuRS)))
subplot(2,3,5)
hold all
plot(unit_wfmean(neuV1 & neufilt & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wfmean(neuV1 & neufilt & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 soft-filtered RS+FS units n=%d + %d', nnz(neuV1 & neufilt & neuRS), nnz(neuV1 & neufilt & ~neuRS)))
subplot(2,3,6)
hold all
plot(unit_wfmean(neuV1 & neufiltxtra & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wfmean(neuV1 & neufiltxtra & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 Siegle-filtered RS+FS units n=%d + %d', nnz(neuV1 & neufiltxtra & neuRS), nnz(neuV1 & neufiltxtra & ~neuRS)))

% normalized spike waveform: RS vs FS
xl=[-30 70];
figure
subplot(2,3,1)
hold all
plot(unit_wftl(neuV1, :)', unit_wfmeannorm(neuV1, :)', 'k-')
plot(unit_wftl(neuV1 & neufilt, :)', unit_wfmeannorm(neuV1 & neufilt, :)', 'b-')
plot(unit_wftl(neuV1 & neufiltxtra, :)', unit_wfmeannorm(neuV1 & neufiltxtra, :)', 'c-')
xlim(xl)
title(sprintf('V1 units n=%d', nnz(neuV1)))
subplot(2,3,2)
hold all
plot(unit_wftl(neuV1 & neufilt, :)', unit_wfmeannorm(neuV1 & neufilt, :)', '-', 'Color', [0 0 0 0.5])
xlim(xl)
title(sprintf('V1 soft-filtered units n=%d', nnz(neuV1 & neufilt)))
subplot(2,3,3)
hold all
plot(unit_wftl(neuV1 & neufiltxtra, :)', unit_wfmeannorm(neuV1 & neufiltxtra, :)', '-', 'Color', [0 0 0 0.5])
xlim(xl)
title(sprintf('V1 Siegle-filtered units n=%d', nnz(neuV1 & neufiltxtra)))
subplot(2,3,4)
hold all
plot(unit_wftl(neuV1 & neuRS, :)', unit_wfmeannorm(neuV1 & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wftl(neuV1 & ~neuRS, :)', unit_wfmeannorm(neuV1 & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 RS+FS units n=%d + %d', nnz(neuV1 & neuRS), nnz(neuV1 & ~neuRS)))
subplot(2,3,5)
hold all
plot(unit_wftl(neuV1 & neufilt & neuRS, :)', unit_wfmeannorm(neuV1 & neufilt & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wftl(neuV1 & neufilt & ~neuRS, :)', unit_wfmeannorm(neuV1 & neufilt & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 soft-filtered RS+FS units n=%d + %d', nnz(neuV1 & neufilt & neuRS), nnz(neuV1 & neufilt & ~neuRS)))
subplot(2,3,6)
hold all
plot(unit_wftl(neuV1 & neufiltxtra & neuRS, :)', unit_wfmeannorm(neuV1 & neufiltxtra & neuRS, :)', '-', 'Color', [0 0 1 0.5])
plot(unit_wftl(neuV1 & neufiltxtra & ~neuRS, :)', unit_wfmeannorm(neuV1 & neufiltxtra & ~neuRS, :)', '-', 'Color', [1 0 0 0.5])
xlim(xl)
title(sprintf('V1 Siegle-filtered RS+FS units n=%d + %d', nnz(neuV1 & neufiltxtra & neuRS), nnz(neuV1 & neufiltxtra & ~neuRS)))
    %}
    %% opto-tagged neurons: SALT p-value + opto-evoked psth sharp peak within 10ms
    load([pathpp 'psth_opto_probeC.mat'])
    probes = {'A', 'B', 'C', 'D', 'E', 'F'};
    unit_saltI = NaN(size(unit_wfdur));
    unit_saltp = NaN(size(unit_wfdur));

    optopsthall = false(length(optopsthtli), length(opto.optotrials), length(unit_wfdur));

    unit_optopulsepsthavg = NaN(length(optopsthtli), length(unit_wfdur));
    unit_optopulsepsthavgbs = NaN(length(optopsthtli), length(unit_wfdur));

    unit_optobasestd = NaN(size(unit_wfdur));
    unit_optopulsepsthpk = NaN(size(unit_wfdur));
    unit_optopulsepsthTpk = NaN(size(unit_wfdur));
    unit_optopulsepsthpkstdx = NaN(size(unit_wfdur));

    unitsaccountedfor = false(size(unit_wfdur));
    for iprobe = 1:numel(probes)
        if ~exist([pathpp 'psth_opto_probe' probes{iprobe} '.mat'])
            fprintf('%s %s Probe %s does not exist!\n', sesid,  opto.genotype, probes{iprobe})
            continue
        else
            fprintf('%s %s Probe %s\n', sesid,  opto.genotype,  probes{iprobe})
        end

        opto_old = opto;
        load([pathpp 'psth_opto_probe' probes{iprobe} '.mat'])
        if ~isequal(opto, opto_old)
            error('why is opto trial info different among sessions?')
        end

        Nneuinprobe = numel(probeunits_saltp);
        [~,neusortsaltp]= sort(probeunits_saltp);
        Nsigp05 = nnz(probeunits_saltp<0.05);
        Nsigp01 = nnz(probeunits_saltp<0.01);

        % ttoind = find(~contains(opto.stimlist, 'cosine'));
        % trialsoi = ismember( opto.optotrials, ttoind);
        optopulsepsthavg = 1000*squeeze(mean(optopsth(:,salttrials,:),2));
        % tlibase = optopsthtli>=-500 & optopsthtli<0;
        tlibase = ismember(optopsthtli, saltbasetli);
        optopulsepsthavgbs = optopulsepsthavg - mean(optopulsepsthavg(tlibase,:),1);

        % additional criterion: there must be a sharp peak within 10ms
        % tli10ms = optopsthtli>=0 & optopsthtli<10;
        tli10ms = ismember(optopsthtli, salttesttli);
        tl10ms = optopsthtli(tli10ms);
        optobasestd = reshape( 1000*std( reshape(optopsth(tlibase,:,:), nnz(tlibase)*size(optopsth,2), size(optopsth,3)), 0,1), [],1);
        optopulsepsthpk = NaN(Nneuinprobe,1);
        optopulsepsthTpk = NaN(Nneuinprobe,1);
        for ci = 1:size(optopulsepsthavg,2)
            [pks,locs]=findpeaks(optopulsepsthavg(tli10ms, ci));
            if ~isempty(pks)
                [mv,mi]=max(pks);
                maxpk = pks(mi);
                maxpkloc = locs(mi);
                optopulsepsthpk(ci) = maxpk;
                optopulsepsthTpk(ci) = tl10ms( maxpkloc );
            end
        end

        optopulsepsthpkstdx = optopulsepsthpk./optobasestd;
        neusigp01pk = find(probeunits_saltp<0.01 & optopulsepsthpkstdx>=1);
        fprintf('non-opto driven units peak %.2f X std\n', max(optopulsepsthpkstdx(probeunits_saltp>=0.05)))
        fprintf('opto driven units p<0.05 n=%d, p<0.01 n=%d, p<0.01 & peak>1*std n=%d\n', nnz(probeunits_saltp<0.05), nnz(probeunits_saltp<0.01), numel(neusigp01pk))

        unit_saltI(neuoind) = probeunits_saltI;
        unit_saltp(neuoind) = probeunits_saltp;
        optopsthall(:,:,neuoind) = optopsth;

        unit_optopulsepsthavg(:,neuoind) = optopulsepsthavg;
        unit_optopulsepsthavgbs(:,neuoind) = optopulsepsthavgbs;
        unit_optobasestd(neuoind) = optobasestd;
        unit_optopulsepsthpk(neuoind) = optopulsepsthpk;
        unit_optopulsepsthTpk(neuoind) = optopulsepsthTpk;

        if all(~unitsaccountedfor(neuoind))
            unitsaccountedfor(neuoind) = true;
        else
            error('why is there an overlap between units in different probes?')
        end
    end
    if ~all(unitsaccountedfor)
        error('failed to account for all units w.r.t optotagging')
    end

    save([pathpp 'psth_opto.mat'], 'opto', 'optopsthtli', 'optopsthall', 'saltbasetli', 'salttesttli', 'salttrials', ...
        'unit_saltI', 'unit_saltp', 'unit_optopulsepsthavg', 'unit_optopulsepsthavgbs', ...
        'unit_optobasestd', 'unit_optopulsepsthpk', 'unit_optopulsepsthTpk', 'unit_optopulsepsthpkstdx', '-v7.3')


    %% plot opto evoked spiking activity
    %{
% example neuron: spike rastergram on each of the 12 opto conditions
ci = neusortsaltp(33);
optostimHz = [1 1 5 10 20 30 40 50 60 80 1 1];
figure
for icond = 1:numel(opto.stimlist)
    trialsoi = opto.optotrials==icond;
    [r,c]=find( squeeze(optopsth(:,trialsoi,ci)) );
subplot(3,4,icond)
hold all
plot(optopsthtli(r), c, 'k.')
set(gca, 'XTick', 0:1000/optostimHz(icond):1000, 'XGrid', 'on')
ylim(0.5+[0 nnz(trialsoi)])
title(opto.stimlist{icond}, 'interpreter', 'none')
end

% heatmap of baseline subtracted mean opto-evoked psth: rows are neurons, columns are time
yl = 0.5+[0 Nneuinprobe];
figure; hold all
imagesc(optopsthtli, 1:Nneuinprobe, optopulsepsthavgbs(:,neusortsaltp)')
plot([0 0], yl, 'k--', 'linewidth', 1)
plot(5+[0 0], yl, 'k--', 'linewidth', 1)
plot([optopsthtli(1) optopsthtli(end)], Nsigp05*[1 1]+0.5, 'k-', 'linewidth', 1)
plot([optopsthtli(1) optopsthtli(end)], Nsigp01*[1 1]+0.5, 'k--', 'linewidth', 1)
set(gca, 'YDir', 'reverse')
xlim([-50 50])
ylim(yl)
caxis(20*[-1 1])
colormap redblue

figure; hold all
scatter(probeunits_saltp, optopulsepsthpkstdx, 10, 'o', 'filled')
plot([0 1], [1 1], 'k-')
plot([0 1], 0.5*[1 1], 'k-')

% spike rastergram of opto-driven units
[sv,si]=sort(probeunits_saltp(neusigp01pk));
neu2plot = neusigp01pk(si);
[~,trialsorted] = sort(opto.optotrials);
yl = 0.5+[0 numel(trialsorted)];
figure;
for ii = 1:min([numel(neu2plot), 15])
    ci = neu2plot(ii);
    [r,c]=find( squeeze(optopsth(:,trialsorted,ci)) );
    subplot(3,5,ii)
    hold all
    plot(optopsthtli(r), c, 'k.')
    plot([0 0], yl, 'c-')
    xlim([-50 50])
    set(gca, 'XTick', 0:100:1000, 'XGrid', 'on', 'YDir', 'reverse')
    ylim(yl)
    title(sprintf('peak %.2f*std p=%.4f', optopulsepsthpkstdx(ci), probeunits_saltp(ci)))
end
    %}
end

%{
1. Re: 강동호 - spontaneous firing rate distribution for E & I cells
2. Re: 나동현 - firing rate distribution for E & P & S & V cells
fit with gamma? log-normal? poisson?
check spike waveform of E & I cells
check optotagging analysis

PCA before and after z-scoring -- compare PC1 coefficient
%}
%%

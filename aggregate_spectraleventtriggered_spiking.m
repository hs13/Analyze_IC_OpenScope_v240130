addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visctxareas = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};

whichneuarea = 'V1';
whichblock = 'ICwcfg1_presentations';

whichblock = 'ICwcfg1_presentations';
blocksplit = strsplit(whichblock, '_');
blockname = blocksplit{1};

%% spike triggered CSD during visual presentation period of each trial type
betaetvisspikeagg = cell(1, Nsessions);
betaetspikeblockagg  = cell(1, Nsessions);
% elecL23visagg = zeros(1, Nsessions);
neuprobelocagg  = cell(1, Nsessions);
Ravgblockagg  = cell(1, Nsessions);
ICsigprobeagg = cell(1, Nsessions);
for ises = 1:Nsessions
    tic
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];

    load([pathpp 'info_electrodes.mat'])
    if strcmp(whichneuarea, 'V1')
        neuareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
    else
        neuarea = visctxareas{strcmp(visareas, whichneuarea)};
        neuareaprobeinds = 1+electrode_probeid( contains(electrode_location, neuarea) );
    end
    if numel(unique(neuareaprobeinds))>1
        warning('this session has more than one %s probe', neuarea)
    end
    neuprobe = mode(neuareaprobeinds);

    fprintf('%d/%d %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, probes{neuprobe})
    if ~exist(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{neuprobe}), 'file')
        fprintf('LFP_TFR_L23_probe%s.mat does not exit!!!\n', probes{neuprobe})
        continue
    end

    load([pathpp 'postprocessed.mat']);%, 'neuallloc')
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{neuprobe} ), 'neuoind')
    neuprobelocagg{ises} = neuallloc(neuoind);
    blocktt = unique(vis.(whichblock).trialorder);
    Ravg = NaN(numel(neuoind), numel(blocktt));
    for typi = 1:numel(blocktt)
        trialsoi = vis.(whichblock).trialorder==blocktt(typi);
        Ravg(:,typi) = mean(Rall.(whichblock)(trialsoi,neuoind),1)';
    end
    Ravgblockagg{ises} = Ravg;

    load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{neuprobe} ))
    ICsigprobeagg{ises} = ICsig;

    load(sprintf('%s%s_betaetspike%s_probe%s.mat', pathpp, whichneuarea, whichblock, probes{neuprobe}))
    betaetvisspikeagg{ises} = betaetvisspike;
    betaetspikeblockagg{ises} = betaetspikeblock;
    toc
end

%%
neuprobelocall = cat(1,neuprobelocagg{:});
neuV1 = contains(neuprobelocall, 'VISp') & ~contains(neuprobelocall, 'VISpm');

temp = cat(1, ICsigprobeagg{:});
ICsigprobe = cat(1, temp.(whichblock));

V1ICenc = neuV1 & cat(1,ICsigprobe(:).ICencoder);
V1indenc = neuV1 & cat(1,ICsigprobe(:).inducerencoder);

betaetvisspikeacc = cat(1,betaetvisspikeagg{:});
betaetspikeblockacc = cat(1,betaetspikeblockagg{:});

Ravgblockacc = cat(1,Ravgblockagg{:});

%%
kerwinhalf = 12; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));

neucol = [0 0.7 0; 1 0 1; 0 0 1; 0 0 0];

figure
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', 'Beta event triggered spiking activity', 'edgecolor', 'none')
for iic = 1:4
        subplot(2,2,iic)
        hold all
    for neuopt = 1:4

switch iic
    case 1
        typi = vistrialtypes==106;
        ttdesc = 'I_C_1';
        switch neuopt
            case 1
                neutitle = 'IC1-encoder';
                neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder1);
            case 2
                neutitle = 'BR+TL';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin1)|cat(1,ICsigprobe(:).indin3));
            case 3
                neutitle = 'IC1-resp';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp1) & ~cat(1,ICsigprobe(:).ICencoder1));
            case 4
                neutitle = 'non-resp';
                neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp2));
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end
    case 2
        typi = vistrialtypes==111;
        ttdesc = 'I_C_2';
        switch neuopt
            case 1
                neutitle = 'IC2-encoder';
                neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder2);
            case 2
                neutitle = 'BL+TR';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin2)|cat(1,ICsigprobe(:).indin4));
            case 3
                neutitle = 'IC2-resp';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp2) & ~cat(1,ICsigprobe(:).ICencoder2));
            case 4
                neutitle = 'non-resp';
                neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp1));
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end
    case 3
        typi = vistrialtypes==0;
        ttdesc = 'Blank';
        switch neuopt
            case 1
                neutitle = 'IC-encoder';
                neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder);
            case 2
                neutitle = 'seg. resp.';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin1)|cat(1,ICsigprobe(:).indin2) ...
                    |cat(1,ICsigprobe(:).indin3)|cat(1,ICsigprobe(:).indin4));
            case 3
                neutitle = 'IC-resp';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp1) | cat(1,ICsigprobe(:).ICresp2) & ~cat(1,ICsigprobe(:).ICencoder));
            case 4
                neutitle = 'non-resp';
                neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp1)) & (~cat(1,ICsigprobe(:).ICresp2));
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end
    case 4
        typi = NaN;
        ttdesc = 'ICwcfg1 block';
        switch neuopt
            case 1
                neutitle = 'IC-encoder';
                neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder);
            case 2
                neutitle = 'seg. resp.';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin1)|cat(1,ICsigprobe(:).indin2) ...
                    |cat(1,ICsigprobe(:).indin3)|cat(1,ICsigprobe(:).indin4));
            case 3
                neutitle = 'IC-resp';
                neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp1) | cat(1,ICsigprobe(:).ICresp2) & ~cat(1,ICsigprobe(:).ICencoder));
            case 4
                neutitle = 'non-resp';
                neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp1)) & (~cat(1,ICsigprobe(:).ICresp2));
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end
    otherwise
        error('iic %d not recognized', iic)
end
if isnan(typi)
    betaettl = allettrange;
tempbetaetbs = 1000*betaetspikeblockacc-Ravgblockacc(:,vistrialtypes==0);
xl = 175*[-1 1];
else
    betaettl = ettrange;
tempbetaetspike = cat(1,betaetvisspikeacc{:,typi});
tempbetaetbs = 1000*tempbetaetspike-Ravgblockacc(:,typi);
xl = 375*[-1 1];
end

% plot(ettrange, convn(tempbetaetspike(neuingroup,:), kergauss, 'same'))
shadedErrorBar(betaettl, mean(convn(tempbetaetbs(neuingroup,:), kergauss, 'same'),1), ...
    std(convn(tempbetaetbs(neuingroup,:), kergauss, 'same'),0,1)/sqrt(nnz(neuingroup)), {'Color', neucol(neuopt,:), 'LineWidth', 1},1)

    end
    plot([ettrange(1) ettrange(end)], [0 0], 'k-')
    xlabel('Time (ms)')
    ylabel('Firing Rate (Hz)')
title([ttdesc ' trials'])
    xlim(xl)
end


%% merge IC1 and IC2 trials
fs = 16;
kerwinhalf = 12; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));

neucol = [0 0.7 0; 1 0 1; 0 0 1; 0 0 0];
neu2plt = [2 1];
neulegs ={'IC-encoders', 'segment responders', 'IC-responsive', 'non-responsive'};

figure('Position', [100 100 400 360])
%annotation('textbox', [0.1 0.92 0.9 0.1], 'string', 'Beta event triggered spiking activity', 'edgecolor', 'none')
hold all
for neuopt = neu2plt
    tempbetaetbsmerge = [];
    for iic = 1:2
        switch iic
            case 1
                typi = vistrialtypes==106;
                ttdesc = 'I_C_1';
                switch neuopt
                    case 1
                        neutitle = 'IC1-encoder';
                        neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder1);
                    case 2
                        neutitle = 'BR+TL';
                        neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin1)|cat(1,ICsigprobe(:).indin3));
                    case 3
                        neutitle = 'IC1-resp';
                        neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp1) & ~cat(1,ICsigprobe(:).ICencoder1));
                    case 4
                        neutitle = 'non-resp';
                        neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp2));
                    otherwise
                        error('neuopt %.0f not recognized', neuopt)
                end
            case 2
                typi = vistrialtypes==111;
                ttdesc = 'I_C_2';
                switch neuopt
                    case 1
                        neutitle = 'IC2-encoder';
                        neuingroup = neuV1 & cat(1,ICsigprobe(:).ICencoder2);
                    case 2
                        neutitle = 'BL+TR';
                        neuingroup = neuV1 & (cat(1,ICsigprobe(:).indin2)|cat(1,ICsigprobe(:).indin4));
                    case 3
                        neutitle = 'IC2-resp';a
                        neuingroup = neuV1 & (cat(1,ICsigprobe(:).ICresp2) & ~cat(1,ICsigprobe(:).ICencoder2));
                    case 4
                        neutitle = 'non-resp';
                        neuingroup = neuV1 & (~cat(1,ICsigprobe(:).ICresp1));
                    otherwise
                        error('neuopt %.0f not recognized', neuopt)
                end
            otherwise
                error('iic %d not recognized', iic)
        end
        if isnan(typi)
            betaettl = allettrange;
            tempbetaetbs = 1000*betaetspikeblockacc-Ravgblockacc(:,vistrialtypes==0);
            xl = 175*[-1 1];
        else
            betaettl = ettrange;
            tempbetaetspike = cat(1,betaetvisspikeacc{:,typi});
            tempbetaetbs = 1000*tempbetaetspike-Ravgblockacc(:,typi);
            xl = 375*[-1 1];
        end
        tempbetaetbsmerge = cat(1, tempbetaetbsmerge, tempbetaetbs(neuingroup,:));
    end
    % plot(ettrange, convn(tempbetaetspike(neuingroup,:), kergauss, 'same'))
    shadedErrorBar(betaettl, mean(convn(tempbetaetbsmerge, kergauss, 'same'),1), ...
        std(convn(tempbetaetbsmerge, kergauss, 'same'),0,1)/sqrt(size(tempbetaetbsmerge,1)), {'Color', neucol(neuopt,:), 'LineWidth', 2},1)

end
set(gca, 'FontSize', fs)
yl = ylim;
cnt = 0;
for neuopt = neu2plt
    cnt = cnt+1;
    text(xl(2),yl(1)+0.08*range(yl)*(cnt), neulegs{neuopt}, 'FontSize', fs, 'Color', neucol(neuopt,:), 'HorizontalAlignment', 'right', 'VerticalAlignment','middle')
end
ylim(yl)
plot([ettrange(1) ettrange(end)], [0 0], 'k-')
plot([0 0], yl, 'k--')
xlabel('Time (ms)', 'FontSize', fs)
ylabel('Firing Rate (Hz)', 'FontSize', fs)
title('Beta_ Event Triggered Spiking Activity', 'FontSize', fs, 'FontWeight', 'normal')% on I_C Trials')
xlim(xl)

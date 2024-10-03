clear all; close all; clc

%addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath(genpath('C:\Users\USER\GitHub\SpectralEvents')) % for TFR
addpath 'C:\Users\USER\GitHub\helperfunctions'

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

probes = {'C', 'D', 'E'};
neuindV1agg = cell(1, numel(nwbsessions));
stLFPagg = cell(numel(probes), numel(nwbsessions));
stCSDagg = cell(numel(probes), numel(nwbsessions));
stTFRagg = cell(numel(probes), numel(nwbsessions));
lfpelecvecagg  = cell(numel(probes), numel(nwbsessions));
ctxelecindsagg = cell(numel(probes), numel(nwbsessions));
elecL23agg = zeros(numel(probes), numel(nwbsessions));
tic
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    for iprobe = 1:numel(probes)
        clearvars('neuindV1', 'stLFPprobe', 'stCSDprobe', 'stTFRprobe')
        if ~exist(sprintf('%sV1spiketriggered_LFP_CSD_TFR_probe%s.mat', pathpp, probes{iprobe}), 'file')
            fprintf('V1spiketriggered_LFP_CSD_TFR_probe%s.mat does not exit!!!\n', probes{iprobe})
            continue
        end
        load( sprintf('%sV1spiketriggered_LFP_CSD_TFR_probe%s.mat', pathpp, probes{iprobe}) )
        if isempty(neuindV1agg{ises})
            neuindV1agg{ises} = neuindV1;
        else
            if ~isequal(neuindV1agg{ises}, neuindV1)
                error('neuindV1 should be the same across probes within the same session')
            end
        end
        stLFPagg{iprobe,ises} = stLFPprobe;
        stCSDagg{iprobe,ises} = stCSDprobe;
        stTFRagg{iprobe,ises} = stTFRprobe;

        load(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{iprobe}), 'lfpelecvec')
        ctxelec = contains(lfpelecvec.location, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        if ~isequal(ctxelecbottom:ctxelectop, find(ctxelec)')
            warning('%d/%d %s cortex electrodes are not consecutive -- check', ises, Nsessions, nwbsessions{ises})
        end
        lfpelecvecagg{iprobe,ises} = lfpelecvec;
        ctxelecindsagg{iprobe,ises} = ctxelec;

        load(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{iprobe}), 'elecL23', 'fVec')
        if ~isempty(elecL23)
            elecL23agg(iprobe,ises) = elecL23;
        else
            elecL23agg(iprobe,ises) = 0;
        end
    end
end
toc

cellfun(@nnz, ctxelecindsagg)

%% low pass filter spike triggered average LFP before calculating CSD
% how much does this differ from low pass filtering LFP before calculating
% spike triggered average LFP?

% y = lowpass(x,wpass) : If x is a matrix, the function filters each column independently.
lfpelecspacing=0.04;

stavglpCSDagg = cell(numel(probes), numel(nwbsessions));
for ises = 1:Nsessions
    for iprobe = 1:numel(probes)
        stavglplfp = NaN(size(stLFPagg{iprobe,ises}));
        tic
        Nneu = size(stLFPagg{iprobe,ises},1);
        Nelec = size(stLFPagg{iprobe,ises},3);
        for ci = 1:Nneu
        tempstlfp = squeeze(stLFPagg{iprobe,ises}(ci,:,:));
        if all(isnan(tempstlfp(:)))
            continue
        end
        stavglplfp(ci,:,:) = lowpass(tempstlfp,100,1000);
        end
        csdelectinds = 2:Nelec-1;
        stavglpcsd = -( stavglplfp(:,:,csdelectinds+1)-2*stavglplfp(:,:,csdelectinds)+stavglplfp(:,:,csdelectinds-1) )/(lfpelecspacing.^2);
        toc
        stavglpCSDagg{iprobe,ises} = stavglpcsd;
    end
end

% figure; hold all
% plot(trange, squeeze(stLFPagg{iprobe,ises}(ci,:,elecL23agg(iprobe,ises))) )
% plot(trange, squeeze(stavglplfp(ci,:,elecL23agg(iprobe,ises))) )

%%
load('S:\OpenScopeData\00248_v240130\postprocessed\openscope_popavg_all.mat')

neuindV1all = false(size(sesneuall));
for ises = 1:Nsessions
    tempsesneuind = find(sesneuall==ises);
    neuindV1all( tempsesneuind(neuindV1agg{ises}) ) = true;
end

neuV1 = contains(neulocall,'VISp') & ~contains(neulocall,'VISpm');
neuRS = unit_wfdur_all>0.4;

if ~isequal(neuindV1all, ismember(sesneuall, find(~cellfun(@isempty, neuindV1agg))) & neuV1 ) 
    error('check neuindV1agg')
end
if ~all( ismember(neuindV1all, neuV1) )
    error('check neuindV1agg')
end

% [v,c]=uniquecnt( sesneuall(ICsigall.ICwcfg1_presentations.ICencoder & neuV1) );
% disp([v,c])

%% IC-encoder vs segment responder CSD averaged in each session
iprobe = find(strcmp(probes, 'E'));

trange = -250:250;
figure
for ises = 1:Nsessions
    if isempty(neuindV1agg{ises}) || nnz(contains(lfpelecvecagg{iprobe,ises}.location, 'VIS'))<=1
        continue
    end
    for neuopt = 1:2
        switch neuopt
            case 1
                neutitle = 'IC-encoder';
                neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
            case 2
                neutitle = 'segment responder';
                neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
                    ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end

        lfpeleclocation = lfpelecvecagg{iprobe,ises}.location;
        Nelec = numel(lfpeleclocation);
        csdelectinds = 2:Nelec-1;
        sesneuoi = neuingroup(sesneuall==ises);
        sesneuindV1 = sesneuoi(neuindV1agg{ises});

        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        yl = [ctxelecbottom ctxelectop]+.5;

        subplot(4,6,2*(ises-1)+neuopt)
        hold all
        imagesc(trange, csdelectinds, squeeze(mean(stCSDagg{iprobe,ises}(sesneuindV1,:,:),1))' )
        scatter(0, elecL23agg(iprobe,ises), 50, 'w*', 'linewidth', 1)
        set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpeleclocation(csdelectinds), 'YDir', 'normal')
        title(sprintf('Session%d %s n=%d', ises, neutitle, nnz(sesneuindV1)))
        ylim(yl)
        xlim([trange(1) trange(end)])
        if neuopt==1
            cl = caxis;
            cl = range(cl)/2 *[-1 1];
        end
        cl = 0.015*[-1 1];
        caxis(cl)
        colorbar
    end
end
colormap jet

%% IC-encoder vs segment responder LFP averaged in each session
trange = -250:250;
figure
for ises = 1:Nsessions
    if isempty(neuindV1agg{ises})
        continue
    end

    for neuopt = 1:2
        switch neuopt
            case 1
                neutitle = 'IC-encoder';
                neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
            case 2
                neutitle = 'segment responder';
                neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
                    ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end

        lfpeleclocation = lfpelecvecagg{iprobe,ises}.location;
        Nelec = numel(lfpeleclocation);
        csdelectinds = 2:Nelec-1;
        sesneuoi = neuingroup(sesneuall==ises);
        sesneuindV1 = sesneuoi(neuindV1agg{ises});

        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        yl = [ctxelecbottom ctxelectop]+.5;


        subplot(4,6,2*(ises-1)+neuopt)
        hold all
        imagesc(trange, 1:Nelec, squeeze(mean(stLFPagg{iprobe,ises}(sesneuindV1,:,:),1))' )
        scatter(0, elecL23agg(iprobe,ises), 50, 'w*', 'linewidth', 1)
        set(gca, 'XGrid', 'on', 'YTick', 1:Nelec, 'YTickLabel', lfpeleclocation, 'YDir', 'normal')
        title(sprintf('Session%d %s n=%d', ises, neutitle, nnz(sesneuindV1)))
        ylim(yl)
        xlim([trange(1) trange(end)])
        %cl = caxis;
        %cl = range(cl)/2 *[-1 1];
        cl = 5*10.^-5*[-1 1];
        caxis(cl)
        colorbar
    end
end
colormap jet

%% IC-encoder vs segment responder TFR averaged in each session
trange = -250:250;
figure
for ises = 1:Nsessions
    if isempty(neuindV1agg{ises})
        continue
    end

    for neuopt = 1:2
        switch neuopt
            case 1
                neutitle = 'IC-encoder';
                neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
            case 2
                neutitle = 'segment responder';
                neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
                    ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end

        sesneuoi = neuingroup(sesneuall==ises);
        sesneuindV1 = sesneuoi(neuindV1agg{ises});

        subplot(4,6,2*(ises-1)+neuopt)
        hold all
        imagesc(trange, fVec, squeeze(mean(stTFRagg{iprobe,ises}(sesneuindV1,:,:),1))' )
        set(gca, 'XGrid', 'on', 'YTick', fVec, 'YTickLabel', fVec, 'YDir', 'normal')
        title(sprintf('Session%d %s n=%d', ises, neutitle, nnz(sesneuindV1)))
        ylim([7 80])
        xlim([trange(1) trange(end)])
        if neuopt==1
        cl = caxis;
        cl(1)=0; cl(2)=0.8*cl(2);
        end
        %cl = [0 3*10^-9];
        caxis(cl)
        colorbar
    end
end
colormap jet

%%
trange = -250:250;
stCSDprobeL23 = NaN(nnz(neuindV1all),1);
stCSDpre50probeL23 = NaN(nnz(neuindV1all),1);
for ises = 1:Nsessions
    if isempty(stCSDagg{iprobe,ises}) | elecL23agg(iprobe,ises)==0
        continue
    end
    tempneuses = sesneuall(neuindV1all)==ises;
    stCSDprobeL23(tempneuses) = stCSDagg{iprobe,ises}(:, trange==0,elecL23agg(iprobe,ises) );
    stCSDpre50probeL23(tempneuses) = sum(stCSDagg{iprobe,ises}(:, trange<=0 & trange>-50,elecL23agg(iprobe,ises) ), 2);
end

figure
hold all
for neuopt = 1:2
    switch neuopt
        case 1
            neutitle = 'IC-encoder';
            neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
        case 2
            neutitle = 'segment responder';
            neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
                ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
        otherwise
            error('neuopt %.0f not recognized', neuopt)
    end
%histogram( stCSDprobeL23(neuingroup(neuindV1all)) , 'normalization', 'pdf', 'binwidth', 0.002)
histogram( stCSDpre50probeL23(neuingroup(neuindV1all)) , 'binwidth', 0.1, 'normalization', 'pdf')
end

figure; histogram( stCSDpre50probeL23)

neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
ICL23sink = stCSDpre50probeL23(neuingroup(neuindV1all));

neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
    ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
segL23sink = stCSDpre50probeL23(neuingroup(neuindV1all));
ranksum(ICL23sink, segL23sink)



%%
fs=14;
VISlayers = {'1', '2/3', '4', '5', '6a', '6b'};
xcols = jet(numel(VISlayers));

figure
for neuopt = 1:2
    switch neuopt
        case 1
            neutitle = 'IC-encoder';
            neuingroup = ICsigall.ICwcfg1_presentations.ICencoder;
        case 2
            neutitle = 'segment responder';
            neuingroup = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
                ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
        otherwise
            error('neuopt %.0f not recognized', neuopt)
    end
    subplot(1,2,neuopt)
    hold all
    for ises = 1:Nsessions
        if isempty(stCSDagg{iprobe,ises})
            continue
        end
        sesneuoi = neuingroup(sesneuall==ises);
        sesneuindV1 = sesneuoi(neuindV1agg{ises});

        lfpeleclocation = lfpelecvecagg{iprobe,ises}.location;
        Nelec = numel(lfpeleclocation);
        csdelectinds = 2:Nelec-1;

        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelecinds = find(ctxelec);
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        [C,ya,yc]=unique(lfpeleclocation(ctxelec));

        xvec = -(find(ctxelec)-ctxelecbottom)/(ctxelectop-ctxelecbottom);
        for il = 1:numel(VISlayers)
            xoi = contains(lfpeleclocation(ctxelec), VISlayers{il});
            plot( xvec(xoi), ...
                squeeze(nanmean(stCSDagg{iprobe,ises}(sesneuindV1, trange<=0 & trange>-100,ctxelecinds(xoi)), 2)), ...
                'o', 'Color', [xcols(il,:) 0.2])%, 'MarkerFaceColor', [xcols(il,:)])
            text(0, yl(2)-0.06*range(yl)*(il-1), VISlayers{il}, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', xcols(il,:), 'FontSize', fs)
        end
    end
    if neuopt==1
        % yl = 0.05*[-1 1];
        % yl = 0.3*[-1 1];
        yl = ylim;
        yl = abs(max(yl))*[-1 1];
    end
    set(gca, 'FontSize',fs)
    ylabel('stCSD -50~0ms from spike')
    xlabel('electrode')
    %legend(VISplayers)
    ylim(yl)
    title(neutitle)
end


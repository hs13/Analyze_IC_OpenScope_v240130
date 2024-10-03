datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visctxareas = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};

lowpassopt = true;
whichneuarea = 'V1';
lfpareas = {'V1', 'LM'};%, 'AL'};%
whichblock = 'ICwcfg1_presentations';

%% spike triggered CSD during visual presentation period of each trial type
stCSDvisagg = cell(numel(lfpareas), Nsessions);
lfpelecvecvisagg  = cell(numel(lfpareas), Nsessions);
ctxelecindsvisagg = cell(numel(lfpareas), Nsessions);
elecL23visagg = zeros(numel(lfpareas), Nsessions);
ICsigV1agg = cell(1, Nsessions);

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


    load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{neuprobe} ))
    ICsigV1agg{ises} = ICsig;

    blocksplit = strsplit(whichblock);
    blockname = blocksplit{1};

    for a = 1:numel(lfpareas)
        lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
        if strcmp(lfpareas{a}, 'V1')
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, 'VISp') & ~contains(electrode_location, 'VISpm') );
        else
            lfparea = visctxareas{strcmp(visareas, lfpareas{a})};
            lfpareaprobeinds = 1+electrode_probeid( contains(electrode_location, lfparea) );
        end
        if numel(unique(lfpareaprobeinds))>1
            warning('this session has more than one %s probe', lfpareas{a})
        end
        lfpprobe = mode(lfpareaprobeinds);

    fprintf('%d/%d %s %s Probe%s\n', ises, numel(nwbsessions), nwbsessions{ises}, lfpareas{a}, probes{lfpprobe})
    if ~exist(sprintf('%sLFP_CSD_probe%s.mat', pathpp, probes{lfpprobe}), 'file')
        fprintf('LFP_CSD_probe%s.mat does not exit!!!\n', probes{lfpprobe})
        continue
    end

    load(sprintf('%sLFP1000Hz_probe%s.mat', pathpp, probes{lfpprobe}), 'lfpelecvec')
    ctxelec = contains(lfpelecvec.location, 'VIS');
    ctxelectop = find(ctxelec, 1, 'last');
    ctxelecbottom = find(ctxelec, 1, 'first');
    if ~isequal(ctxelecbottom:ctxelectop, find(ctxelec)')
        warning('%d/%d %s cortex electrodes are not consecutive -- check', ises, Nsessions, nwbsessions{ises})
    end
    lfpelecvecvisagg{a,ises} = lfpelecvec;
    ctxelecindsvisagg{a,ises} = ctxelec;

    load(sprintf('%sLFP_TFR_L23_probe%s.mat', pathpp, probes{lfpprobe}), 'elecL23', 'fVec')
    if ~isempty(elecL23)
        elecL23visagg(a,ises) = elecL23;
    else
        elecL23visagg(a,ises) = 0;
    end

    if lowpassopt
        stpsthCSDfn = sprintf('%s%s_stCSD%s_lowpassprobe%s.mat', pathpp, whichneuarea, blockname, probes{lfpprobe});
    else
        stpsthCSDfn = sprintf('%s%s_stCSD%s_probe%s.mat', pathpp, whichneuarea, blockname, probes{lfpprobe});
    end
    %     'lfpelecspacing', 'csdelectinds', 'vistrialtypes', 'sttrange', 'stCSDvisprobe'
    load(stpsthCSDfn)
    stCSDvisagg{a,ises} = stCSDvisprobe;
    end

    toc
end


%% IC1 trials: IC-encoder vs segment responder CSD averaged in each session
a = find(strcmp(lfpareas,'LM'));

typi = vistrialtypes==106;
figure
for ises = 1:Nsessions
    if isempty(stCSDvisagg{ises})
        continue
    end
    for neuopt = 1:2
        switch neuopt
            case 1
                neutitle = 'IC1-encoder';
                neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.ICencoder1;
            case 2
                neutitle = 'BR+TL';
                neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indin1 | ...
                    ICsigV1agg{ises}.ICwcfg1_presentations.indin3;
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end

        lfpeleclocation = lfpelecvecvisagg{a,ises}.location;
        Nelec = numel(lfpeleclocation);
        csdelectinds = 2:Nelec-1;

        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        yl = [ctxelecbottom ctxelectop]+.5;

        subplot(4,6,2*(ises-1)+neuopt)
        hold all
        imagesc(sttrange, csdelectinds, squeeze(mean(stCSDvisagg{a,ises}{typi}(neuingroup,:,:),1))' )
        scatter(0, elecL23visagg(a,ises), 50, 'w*', 'linewidth', 1)
        set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpeleclocation(csdelectinds), 'YDir', 'normal')
        title(sprintf('Session%d %s n=%d %s CSD', ises, neutitle, nnz(neuingroup), lfpareas{a}))
        ylim(yl)
        xlim([sttrange(1) sttrange(end)])
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

% %% IC2 trials: IC-encoder vs segment responder CSD averaged in each session
typi = vistrialtypes==111;
figure
for ises = 1:Nsessions
    if isempty(stCSDvisagg{ises})
        continue
    end
    for neuopt = 1:2
        switch neuopt
            case 1
                neutitle = 'IC2-encoder';
                neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.ICencoder2;
            case 2
                neutitle = 'BL+TR';
                neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indin2 | ...
                    ICsigV1agg{ises}.ICwcfg1_presentations.indin4;
            otherwise
                error('neuopt %.0f not recognized', neuopt)
        end

        lfpeleclocation = lfpelecvecvisagg{a,ises}.location;
        Nelec = numel(lfpeleclocation);
        csdelectinds = 2:Nelec-1;

        ctxelec = contains(lfpeleclocation, 'VIS');
        ctxelectop = find(ctxelec, 1, 'last');
        ctxelecbottom = find(ctxelec, 1, 'first');
        yl = [ctxelecbottom ctxelectop]+.5;

        subplot(4,6,2*(ises-1)+neuopt)
        hold all
        imagesc(sttrange, csdelectinds, squeeze(mean(stCSDvisagg{a,ises}{typi}(neuingroup,:,:),1))' )
        scatter(0, elecL23visagg(a,ises), 50, 'w*', 'linewidth', 1)
        set(gca, 'XGrid', 'on', 'YTick', csdelectinds, 'YTickLabel', lfpeleclocation(csdelectinds), 'YDir', 'normal')
        title(sprintf('Session%d %s n=%d %s CSD', ises, neutitle, nnz(neuingroup), lfpareas{a}))
        ylim(yl)
        xlim([sttrange(1) sttrange(end)])
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


%% IC-encoder vs segment responder V1L23 CSD -50~0ms relative to spike
ICsigV1all = cat(1, ICsigV1agg{:} );
ICwcfg1sigV1all = cat(1, ICsigV1all.ICwcfg1_presentations);
NneuV1all = length(cat(1,ICwcfg1sigV1all.ICencoder));

stCSDpre50V1L23 = NaN(NneuV1all,numel(vistrialtypes));
neucnt = 0;
for ises = 1:Nsessions
    if isempty(stCSDvisagg{a,ises})
        continue
    end    
    tempneuses = neucnt+1 : neucnt+size(stCSDvisagg{a,ises}{1},1);
    for typi = 1:numel(vistrialtypes)
    stCSDpre50V1L23(tempneuses, typi) = mean(stCSDvisagg{a,ises}{typi}(:, ...
        sttrange<=0 & sttrange>-50,elecL23visagg(a,ises) ), 2);
    end
    neucnt = neucnt + size(stCSDvisagg{a,ises}{typi},1);
end

if neucnt ~= NneuV1all
    error('not all neurons accounted for')
end

neugrousoi = {'ICencoder1', 'ICencoder2', 'indin1', 'indin2', 'indin3', 'indin4'};
stCSDpre50 = struct();
for g = 1:numel(neugrousoi)
    neuoi = cat(1,ICwcfg1sigV1all.(neugrousoi{g}));
    stCSDpre50.(neugrousoi{g}) =  stCSDpre50V1L23(neuoi,:);
end

IC1IC1 = stCSDpre50.ICencoder1(:,vistrialtypes==106);
IC1SR = cat(1,stCSDpre50.indin1(:,vistrialtypes==106),stCSDpre50.indin3(:,vistrialtypes==106));
IC2IC2 = stCSDpre50.ICencoder2(:,vistrialtypes==111);
IC2SR = cat(1,stCSDpre50.indin2(:,vistrialtypes==111),stCSDpre50.indin3(:,vistrialtypes==111));

ranksum(IC1IC1, IC1SR)
ranksum(IC2IC2, IC2SR)
ranksum([IC1IC1; IC2IC2], [IC1SR; IC2SR])

figure
hold all
h=histogram([IC1SR; IC2SR]);%, 'binwidth');
histogram([IC1IC1; IC2IC2], 'BinWidth', h.BinWidth)


%% IC-encoder vs segment responder V1L4 CSD -50~0ms relative to spike
ICsigV1all = cat(1, ICsigV1agg{:} );
ICwcfg1sigV1all = cat(1, ICsigV1all.ICwcfg1_presentations);
NneuV1all = length(cat(1,ICwcfg1sigV1all.ICencoder));

stCSDpre50V1L4 = NaN(NneuV1all,numel(vistrialtypes));
neucnt = 0;
for ises = 1:Nsessions
    if isempty(stCSDvisagg{a,ises})
        continue
    end    
    elecsL4 = strcmp(lfpelecvecvisagg{a,ises}.location, 'VISp4');
    tempneuses = neucnt+1 : neucnt+size(stCSDvisagg{a,ises}{1},1);
    for typi = 1:numel(vistrialtypes)
    stCSDpre50V1L4(tempneuses, typi) = mean(stCSDvisagg{a,ises}{typi}(:, sttrange<=0 & sttrange>-50,elecsL4 ), [2 3]);
    end
    neucnt = neucnt + size(stCSDvisagg{a,ises}{typi},1);
end

if neucnt ~= NneuV1all
    error('not all neurons accounted for')
end

neugrousoi = {'ICencoder1', 'ICencoder2', 'indin1', 'indin2', 'indin3', 'indin4'};
stCSDpre50 = struct();
for g = 1:numel(neugrousoi)
    neuoi = cat(1,ICwcfg1sigV1all.(neugrousoi{g}));
    stCSDpre50.(neugrousoi{g}) =  stCSDpre50V1L4(neuoi,:);
end

IC1IC1 = stCSDpre50.ICencoder1(:,vistrialtypes==106);
IC1SR = cat(1,stCSDpre50.indin1(:,vistrialtypes==106),stCSDpre50.indin3(:,vistrialtypes==106));
IC2IC2 = stCSDpre50.ICencoder2(:,vistrialtypes==111);
IC2SR = cat(1,stCSDpre50.indin2(:,vistrialtypes==111),stCSDpre50.indin3(:,vistrialtypes==111));

ranksum(IC1IC1, IC1SR)
ranksum(IC2IC2, IC2SR)
ranksum([IC1IC1; IC2IC2], [IC1SR; IC2SR])

figure
hold all
h=histogram([IC1SR; IC2SR]);%, 'binwidth');
histogram([IC1IC1; IC2IC2], 'BinWidth', h.BinWidth)

%%
fs=14;
VISlayers = {'1', '2/3', '4', '5', '6a', '6b'};
xcols = jet(numel(VISlayers));

figure
for iic = 1:3
    for neuopt = 1:3
        subplot(3,3,3*(iic-1)+neuopt)
        hold all
        for ises = 1:Nsessions
            if isempty(stCSDvisagg{a,ises})
                continue
            end
            switch iic
                case 1
                    typi = vistrialtypes==106;
                    ttdesc = 'I_C_1';
                    switch neuopt
                        case 1
                            neutitle = 'IC1-encoder';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.ICencoder1;
                        case 2
                            neutitle = 'BR+TL';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indin1 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indin3;
                        case 3
                            neutitle = 'inducerencoder1+3';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indenc1 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indenc3;
                        otherwise
                            error('neuopt %.0f not recognized', neuopt)
                    end
                case 2
                    typi = vistrialtypes==111;
                    ttdesc = 'I_C_2';
                    switch neuopt
                        case 1
                            neutitle = 'IC2-encoder';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.ICencoder2;
                        case 2
                            neutitle = 'BL+TR';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indin2 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indin4;
                        case 3
                            neutitle = 'inducerencoder2+4';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indenc2 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indenc4;
                        otherwise
                            error('neuopt %.0f not recognized', neuopt)
                    end
                case 3
                    typi = vistrialtypes==0;
                    ttdesc = 'Blank';
                    switch neuopt
                        case 1
                            neutitle = 'IC-encoder';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.ICencoder;
                        case 2
                            neutitle = 'seg. resp.';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indin1 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indin2 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indin3 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indin4;
                        case 3
                            neutitle = 'inducerencoder';
                            neuingroup = ICsigV1agg{ises}.ICwcfg1_presentations.indenc1 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indenc2 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indenc3 | ...
                                ICsigV1agg{ises}.ICwcfg1_presentations.indenc4;
                        otherwise
                            error('neuopt %.0f not recognized', neuopt)
                    end
                otherwise
                            error('iic %d not recognized', iic)
            end

            Nelec = numel(lfpelecvecvisagg{a,ises}.location);
            csdelectinds = 2:Nelec-1;
            csdeleclocation = lfpelecvecvisagg{a,ises}.location(csdelectinds);

            ctxelec = contains(csdeleclocation, 'VIS');
            ctxelecinds = find(ctxelec);
            ctxelectop = find(ctxelec, 1, 'last');
            ctxelecbottom = find(ctxelec, 1, 'first');
            [C,ya,yc]=unique(csdeleclocation(ctxelec));
            xvec = -(find(ctxelec)-ctxelecbottom)/(ctxelectop-ctxelecbottom);

            switch lfpareas{a}
                case 'V1'
            yl = 0.6*[-1 1];
                case 'LM'
            yl = 0.25*[-1 1];
                otherwise
                    yl = ylim;
            end 


            hold on
            for il = 1:numel(VISlayers)
                xoi = contains(csdeleclocation(ctxelec), VISlayers{il});
                layleg = unique(csdeleclocation(ctxelecinds(xoi)));
                plot( xvec(xoi), ...
                    squeeze(nanmean(stCSDvisagg{a,ises}{typi}(neuingroup, sttrange<=0 & sttrange>-1,ctxelecinds(xoi)), 2)), ...
                    'o', 'Color', xcols(il,:), 'MarkerFaceColor', xcols(il,:))
                text(0, yl(2)-0.06*range(yl)*(il-1), layleg, ... %VISlayers{il}, ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Color', xcols(il,:), 'FontSize', fs)
            end
        end
            set(gca, 'FontSize',fs)
            ylabel('stCSD 0ms from spike')
            xlabel('electrode')
            %legend(VISplayers)
            ylim(yl)
            title(sprintf('%s trials %s', ttdesc, neutitle))        
    end
end

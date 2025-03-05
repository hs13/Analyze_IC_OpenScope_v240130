%{
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
    clearvars -except ises nwbsessions
    sesclk = tic;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    load([pathpp, 'postprocessed.mat'])
    load([pathpp, 'psth_spontaneous.mat'])
    load([pathpp, 'qc_units.mat'])

    neuV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
    neuvis = contains(neuallloc, 'VIS');
    neuRS = unit_wfdur>0.4;
    neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);

    neuV1RS = neuV1 & neuRS;
    %neuvisRS = neuvis & neuRS;

    whichblock = 'ICwcfg1_presentations';
    trialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
    tempR = 0.4*Rall.(whichblock)(:,neuV1RS)';

    spondurs = vis.spontaneous_presentations.stop_time-vis.spontaneous_presentations.start_time;
    [mv,mi]=max(spondurs);
    whichsponind = find(vis.spontaneous_presentations.start_time==vis.ICwcfg1_presentations.stop_time(end));
    fprintf('spontaneous bock after ICwcfg1 is %.3fs long\n', spondurs(whichsponind))
    if whichsponind~=mi
        warning('spontaneous bock after ICwcfg1 is not the longest spontaneous block')
    end

    % preprocess spontaneous data
    % tempR is 400ms-window spike count formatted #neurons * #trials
    temppsth = psthspon{whichsponind}(:,neuV1RS);
    Ttot = size(temppsth,1);
    Twin = 400; % Tres must be 0.001s (1ms)
    Tslide = 25;
    Ntwins = floor( (Ttot-Twin)/Tslide );
    Tstartind = mod(Ttot-Twin, Tslide);
    Tspkinds = Tstartind+( (1:Twin)'+(0:Tslide:Ttot-Twin) );
    if Tspkinds(end,end) ~= Ttot
        error('check code')
    end
    Tctr = Tspkinds(round(Twin/2),:);
    psthsponspkcnt = zeros( nnz(neuV1RS), size(Tspkinds,2) );
    for ci = 1:nnz(neuV1RS)
        temppsthvec = temppsth(:,ci);
        psthsponspkcnt(ci,:) = sum(temppsthvec(Tspkinds), 1);
    end

    save([pathpp, 'psth_spontaneous_spikecount_V1RS.mat'], 'whichsponind', 'neuV1RS', 'Tspkinds', 'psthsponspkcnt')
end
%}

%%
mousedate = 'sub-619296';
pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
load([pathpp 'UMAP_KNN_decoding_V1RS_spontaneous.mat'])
load([pathpp 'UMAP_KNN_2fold_decoding_V1RS_spontaneous.mat'])
load([pathpp 'UMAP_KNN_decoding_V1RS_spontaneous_embeddims.mat'])
load([pathpp 'psth_spontaneous_spikecount_V1RS.mat'])
load([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'], 'trialorder')
% load([pathpp, 'psth_spontaneous.mat'])
% not needed yet. load if we want to see spontaneous reactiation at a
% higher temporal resolution

hireptt = unique(trialorder);
Nhireptt = numel(hireptt);
Nsplits = size(UMAPKNN_semisup.confusion_matrix,1);
Nneurons = nnz(neuV1RS);

trialcol = zeros(Nhireptt,3);
trialcol(hireptt==0,:) = [0 0 0];
trialcol(hireptt==101,:) = [1 0 0];
trialcol(hireptt==105,:) = [0 1 1];
trialcol(hireptt==109,:) = [0 0.5 0.5];
trialcol(hireptt==106,:) = [0 1 0];
trialcol(hireptt==111,:) = [0 0.5 0];
trialcol(hireptt==107,:) = [1 0.5 0];
trialcol(hireptt==110,:) = [0.5 0.25 0];
trialcol(hireptt==1105,:) = [0 0 1];
trialcol(hireptt==1109,:) = [0 0 0.5];
trialcol(hireptt==1201,:) = [1 0 1];
trialcol(hireptt==1299,:) = [0.5 0 0.5];

% TODO: COMPARE WITH SVM
% whichsponind=2;
% load(strcat(pathpp, 'SVMspon', num2str(whichsponind), '_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat'))

fs = 14;
%% embedding dimensions
cvaccmat_dims = zeros(Nsplits, numel(embedding_dimensions));
for idim = 1:numel(embedding_dimensions)
    cvaccmat_dims(:,idim) = UMAPKNN_dimensions{idim}.accuracy;
end

figure; hold all
plot(embedding_dimensions, cvaccmat_dims', '.')
errorbar(embedding_dimensions, mean(cvaccmat_dims,1), std(cvaccmat_dims,0,1)/sqrt(Nsplits), 'ko-', 'LineWidth', 2)
set(gca, 'FontSize', fs)
xlabel('UMAP Embedding Dimensions')
ylabel('10-fold cross-validated accuracy')
title([mousedate ' unsupervised UMAP+KNN decoding'])

[p,tbl,stats]=kruskalwallis(cvaccmat_dims);
multcompare(stats)

%% compare test accuracy between UMAP settings
cvaccmat = [UMAPKNN_semisup.accuracy' UMAPKNN_unsup.accuracy' UMAPKNN_semisup_invcv.accuracy' UMAPKNN_unsup_invcv.accuracy' ];
[p,tbl,stats]=kruskalwallis(cvaccmat);
multcompare(stats)
figure; hold all
plot(1:4, cvaccmat', '.')
errorbar(1:4, mean(cvaccmat,1), std(cvaccmat,0,1)/sqrt(Nsplits), 'ko-', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'semi-sup. 10-fold CV', 'unsup. 10-fold CV', 'semi-sup. inv.-CV', 'unsup. inv.-CV'}, 'XTickLabelRotation', 45)

%% check consistency of spontaneous reactivations between UMAP settings
threshcv = 0.9; % out of Nsplits decoders, how many agreed to the same solution?
threshinvcv = 0.6;

% visualize spontaneous reactivation
spondecode = struct();
spondecode.semisup = zeros(length(Tspkinds),Nhireptt);
spondecode.unsup = zeros(length(Tspkinds),Nhireptt);
spondecode.semisup_invcv = zeros(length(Tspkinds),Nhireptt);
spondecode.unsup_invcv = zeros(length(Tspkinds),Nhireptt);
spondecode.semisup_2fold = zeros(length(Tspkinds),Nhireptt);
spondecode.unsup_2fold = zeros(length(Tspkinds),Nhireptt);
for itt = 1:Nhireptt
    spondecode.semisup(:,itt) = mean(UMAPKNN_semisup.probe_predlabels==hireptt(itt),1);
    spondecode.unsup(:,itt) = mean(UMAPKNN_unsup.probe_predlabels==hireptt(itt),1);
    spondecode.semisup_invcv(:,itt) = mean(UMAPKNN_semisup_invcv.probe_predlabels==hireptt(itt),1);
    spondecode.unsup_invcv(:,itt) = mean(UMAPKNN_unsup_invcv.probe_predlabels==hireptt(itt),1);
    spondecode.semisup_2fold(:,itt) = mean(UMAPKNN2_semisup.probe_predlabels==hireptt(itt),1);
    spondecode.unsup_2fold(:,itt) = mean(UMAPKNN2_unsup.probe_predlabels==hireptt(itt),1);
end

decodersettings = fieldnames(spondecode);
decoderreactsurvival = zeros(numel(decodersettings), Nsplits+1);
for f = 1:numel(decodersettings)
    for ithr = 0:Nsplits
        decoderreactsurvival(f,ithr+1) = mean(any(spondecode.(decodersettings{f})>=ithr/Nsplits,2));
    end
end

decodersreactratio = zeros(numel(decodersettings),Nhireptt);
decoderreactimgind = struct();
for f = 1:numel(decodersettings)
    if contains(decodersettings{f}, 'invcv')
        tempthr = threshinvcv;
    else
        tempthr = threshcv;
    end
    decodersreactratio(f,:) = mean(spondecode.(decodersettings{f})>=tempthr,1);

    suprathreshtp = any(spondecode.(decodersettings{f})>=tempthr,2);
    [mv,mi] = max(spondecode.(decodersettings{f}),[],2);
    decoderreactimgind.(decodersettings{f}) = mi;
    decoderreactimgind.(decodersettings{f})(~suprathreshtp) = 0;
end

% unsupervised has better consistency than semisupervised
% 10-fold cross validation has more consistency between semi-supervised and unsupervised
% winner: unsupervised 10-fold cross-validation method
Nmatchtimepoints = NaN(numel(decodersettings));
Pmatchtimepoints = NaN(numel(decodersettings));
Prefmatchtimepoints = NaN(numel(decodersettings));
for f = 1:numel(decodersettings)
    for g = 1:numel(decodersettings)
        supthreshtp = decoderreactimgind.(decodersettings{f})>0 | decoderreactimgind.(decodersettings{g})>0;
        Nmatchtimepoints(f,g) = nnz(decoderreactimgind.(decodersettings{f})(supthreshtp)==decoderreactimgind.(decodersettings{g})(supthreshtp));
        Pmatchtimepoints(f,g) = mean(decoderreactimgind.(decodersettings{f})(supthreshtp)==decoderreactimgind.(decodersettings{g})(supthreshtp));
        refsupthreshtp = decoderreactimgind.(decodersettings{f})>0;
        Prefmatchtimepoints(f,g) = mean(decoderreactimgind.(decodersettings{f})(refsupthreshtp)==decoderreactimgind.(decodersettings{g})(refsupthreshtp));
    end
end

reactcoeffuk = struct();
for f = 1:numel(decodersettings)
    tempcorr = corr(psthsponspkcnt', spondecode.(decodersettings{f}));
    reactcoeffuk.(decodersettings{f}) = tempcorr;
end


%%
figure; plot((0:Nsplits)/Nsplits, decoderreactsurvival, 'linewidth', 1)
legend(decodersettings, 'FontSize', fs, 'interpreter', 'none', 'location', 'best')
set(gca, 'XTick', 0:0.1:1, 'XGrid', 'on', 'YTick', 0:0.05:1, 'YGrid', 'on', 'FontSize', fs)
xlabel('Spontaneous ReactivationThreshold')
ylabel('Portion Timepoints')

figure
for f = 1:numel(decodersettings)
    subplot(numel(decodersettings),1,f)
    imagesc(spondecode.(decodersettings{f})')
    caxis([0 1])
    title(decodersettings{f}, 'interpreter', 'none')
end
colormap redblue


figure
imagesc(decodersreactratio)
colorbar
% caxis([0 1])
set(gca, 'XTick', 1:Nhireptt, 'XTickLabel', hireptt, 'YTick', 1:4, 'YTickLabel', decodersettings)
colormap redblue

disp([hireptt; 100*decodersreactratio])

%% todo: REVISE THIS CODE to compare with SVM
% Nmatchsvm = NaN(1,numel(decodersettings));
% Pmatchsvm = NaN(1,numel(decodersettings));
% Prefmatchsvm = NaN(1,numel(decodersettings));
% for f = 1:numel(decodersettings)
%     supthreshtp = decoderreactimgind.(decodersettings{f})>0 | decoderreactimgind.(decodersettings{g})>0;
%     Nmatchsvm(f,g) = nnz(decoderreactimgind.(decodersettings{f})(supthreshtp)==decoderreactimgind.(decodersettings{g})(supthreshtp));
%     Pmatchsvm(f,g) = mean(decoderreactimgind.(decodersettings{f})(supthreshtp)==decoderreactimgind.(decodersettings{g})(supthreshtp));
%     refsupthreshtp = decoderreactimgind.(decodersettings{f})>0;
%     Prefmatchsvm(f,g) = mean(decoderreactimgind.(decodersettings{f})(refsupthreshtp)==decoderreactimgind.(decodersettings{g})(refsupthreshtp));
% end

%% visualize embedding
isplit = 1;
spinds = reshape(1:6,3,2)';
figure
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', sprintf('%s UMAP+KNN decoding: boundaries demarcate 95%% CI of test trials, dots denote spontaneous reactivations', mousedate), 'EdgeColor', 'none', 'FontSize', fs)
for f = 1:numel(decodersettings)
    if contains(decodersettings{f}, 'invcv')
        tempthr = threshinvcv;
    else
        tempthr = threshcv;
    end
    switch decodersettings{f}
        case 'semisup'
            tempUMAPKNN = UMAPKNN_semisup;
        case 'unsup'
            tempUMAPKNN = UMAPKNN_unsup;
        case 'semisup_invcv'
            tempUMAPKNN = UMAPKNN_semisup_invcv;
        case 'unsup_invcv'
            tempUMAPKNN = UMAPKNN_unsup_invcv;
        case 'semisup_2fold'
            tempUMAPKNN = UMAPKNN2_semisup;
        case 'unsup_2fold'
            tempUMAPKNN = UMAPKNN2_unsup;
        otherwise
            error('decoder setting not recognized')
    end
    subplot(2,3,spinds(f))
    hold all
    for itt = 1:Nhireptt
        trialsoi = find( tempUMAPKNN.test_truelabels(isplit,:)==hireptt(itt) );
        P = double( squeeze(tempUMAPKNN.test_embeddings(isplit,trialsoi,:)) );
        distances = mahal(P, P);
        % Define a threshold for outlier detection (e.g., 95% confidence level)
        threshold = chi2inv(0.95, size(P, 2)); % Chi-square critical value
        outliers = distances > threshold;
        % Remove outliers
        P_cleaned = P(~outliers, :);
        [k,v]=boundary(P_cleaned);

        % plot(P(:,1),P(:,2), '.', 'Color', trialcol(itt,:))
        %     plot(P_cleaned(k,1),P_cleaned(k,2), 'FaceColor', [trialcol(itt,:) 0.5])
        % fill(P_cleaned(k,1),P_cleaned(k,2), trialcol(itt,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(P_cleaned(k,1),P_cleaned(k,2), '-', 'Color', trialcol(itt,:), 'LineWidth', 1);

        sponreacttp = decoderreactimgind.(decodersettings{f})==itt;
        tempx = squeeze(tempUMAPKNN.probe_embeddings(isplit,sponreacttp,1));
        tempy = squeeze(tempUMAPKNN.probe_embeddings(isplit,sponreacttp,2));
        plot(tempx, tempy, '.', 'Color', trialcol(itt,:))
    end
    xlabel('UMAP Dim.1', 'FontSize', fs)
    ylabel('UMAP Dim.2', 'FontSize', fs)
    title( sprintf('%s reactivation threshold>=%.2f', decodersettings{f}, tempthr), 'interpreter', 'none', 'FontSize', fs)
end

%% find activity patterns that come within multi-dimensional 95% CI of test trials 
% alternative method for identifying spontaneous reactivations
%% align activity to spontaneous reactivations: 10-fold unsupervised UMAP+KNN
fprintf('spon. react. threshold %.1f\n', threshcv)
Treact = struct();
Treact.startind = cell(Nhireptt,1);
Treact.endind = cell(Nhireptt,1);
Treact.duration = cell(Nhireptt,1);
Treact.startms = cell(Nhireptt,1);
Treact.endms = cell(Nhireptt,1);

% note that we are z-scoring smoothed activity trace
psthsponreacttli = (-100:100)'; % each bin is 25ms; 2500ms before and after onset of reactivation
sponspkcntmean = mean(psthsponspkcnt,2);
sponspkcntstd = std(psthsponspkcnt,0,2);
Tsponlen = size(psthsponspkcnt,2);
psthsponzrate = (psthsponspkcnt-sponspkcntmean)./sponspkcntstd;
psthsponreactzrate = cell(Nhireptt,1);
for itt = 1:Nhireptt
    if nnz(decoderreactimgind.unsup==itt)==0
        continue
    end
    tempreacttinds=find(decoderreactimgind.unsup==itt);
    tempdtinds = find( diff(tempreacttinds)>1 ); % find timeponts that are not consecutive
    tempreactstarttind = [tempreacttinds(1); tempreacttinds(tempdtinds+1)];
    tempreactendtind = [tempreacttinds(tempdtinds); tempreacttinds(end)];

    if ~isequal(tempreactendtind-tempreactstarttind+1, diff([0; tempdtinds; length(tempreacttinds)]) )
        error('check tempdtinds')
    end
    tempreacttwin = tempreactendtind-tempreactstarttind+1;
    fprintf('spon. react. %d: all n=%d, consec.tp n=%d\n', hireptt(itt), ...
        numel(tempreacttwin), nnz(tempreacttwin>1))
    % figure; histogram(tempreacttwin)

    Treact.startind{itt} = tempreactstarttind;
    Treact.endind{itt} = tempreactendtind;
    Treact.duration{itt} = tempreacttwin;
    Treact.startms{itt} = Tspkinds(1,tempreactstarttind);
    Treact.endms{itt} = Tspkinds(end,tempreactendtind);

    tempNreacts = length(tempreactstarttind);
    temppsthtinds = reshape(tempreactstarttind,1,[])+psthsponreacttli;
    invalidtinds = temppsthtinds<=0 | temppsthtinds>Tsponlen;
    temppsthtinds(temppsthtinds<=0) = 1;
    temppsthtinds(temppsthtinds>Tsponlen) = Tsponlen;
    if ~(size(temppsthtinds,1)==length(psthsponreacttli) && size(temppsthtinds,2)==tempNreacts )
        disp(size(temppsthtinds))
        error('check temppsthtinds')
    end
    temppsthsponreactzrate = NaN(length(psthsponreacttli),tempNreacts,Nneurons);
    for ci = 1:Nneurons
        temppsthvec = psthsponzrate(ci,:);
        temppsthmat = temppsthvec(temppsthtinds);
        temppsthmat(invalidtinds) = NaN;
        temppsthsponreactzrate(:,:,ci) = temppsthmat;
    end
    psthsponreactzrate{itt} = temppsthsponreactzrate;
end

% identify spon. react. ensembles by calculating auroc between activity
% during reactivation events vs unclassified spontaneous activity (all
% classes are below 0.5)
Nboot = 0; % takes 18min for Nboot=200;
unclasstinds = find(all(spondecode.unsup<0.5,2));
neuenssponreact = struct();
neuenssponreact.unclasstinds= unclasstinds;
if Nboot==0
    neuenssponreact.auroc = NaN(Nneurons, Nhireptt);%,3);
else
    neuenssponreact.auroc = NaN(Nneurons, Nhireptt,3);
end
neuenssponreact.Pranksum = NaN(Nneurons, Nhireptt);
tic
for itt = 1:Nhireptt
    if nnz(decoderreactimgind.unsup==itt)==0
        continue
    end
    tempreacttinds=find(decoderreactimgind.unsup==itt);
    for ci = 1:Nneurons
        [p,h,stats] = ranksum(psthsponspkcnt(ci,tempreacttinds), psthsponspkcnt(ci,unclasstinds));
        if Nboot==0
        [X,Y,T,AUC] = perfcurve([ones(1,numel(tempreacttinds)), zeros(1,numel(unclasstinds))], ...
            [psthsponspkcnt(ci,tempreacttinds) psthsponspkcnt(ci,unclasstinds)], '1');
        neuenssponreact.auroc(ci,itt) = AUC;
        else
        [X,Y,T,AUC] = perfcurve([ones(1,numel(tempreacttinds)), zeros(1,numel(unclasstinds))], ...
            [psthsponspkcnt(ci,tempreacttinds) psthsponspkcnt(ci,unclasstinds)], '1', 'NBoot', Nboot);
        neuenssponreact.auroc(ci,itt,:) = AUC;
        end
        neuenssponreact.Pranksum(ci,itt) = p;
    end
    toc
end

%%
itt = find(hireptt==106);
if Nboot==0
    sponreactauroc = neuenssponreact.auroc(:,itt);
tempsigneu = neuenssponreact.Pranksum(:,itt)<0.05;
else
    % bootstrapping 200 times was *less* stringent than ranksum test!
    sponreactauroc = squeeze(neuenssponreact.auroc(:,itt,1));
    tempsigact = squeeze(all(neuenssponreact.auroc(:,itt,:)>0.5,3));
    tempsigsup = squeeze(all(neuenssponreact.auroc(:,itt,:)<0.5,3));
tempsigneu = tempsigact | tempsigsup;
end
figure; hold all
plot(sponreactauroc, reactcoeffuk.unsup(:,itt), 'o')
plot(sponreactauroc(tempsigneu), reactcoeffuk.unsup(tempsigneu,itt), 'r*')
xl=xlim; yl = ylim;
plot([0 1], [0 0], '-', 'Color', 0.5+[0 0 0])
plot([0.5 0.5], [-1 1], '-', 'Color', 0.5+[0 0 0])
set(gca, 'XGrid', 'on', 'YGrid', 'on')
axis([xl yl])
tempsigmww = neuenssponreact.Pranksum(:,itt)<0.05;


tempneuens = neuenssponreact.Pranksum(:,itt)<0.05 & ...
    neuenssponreact.auroc(:,itt)>0 & reactcoeffuk.unsup(:,itt)>0;
figure
plot(psthsponreacttli, squeeze(nanmean(psthsponreactzrate{itt}(:,:,tempneuens),2)) )

tbinms = 25;
psthsponreacttimeline = tbinms*psthsponreacttli;
t0ind = find(psthsponreacttli==0);

xt = 1:20:length(psthsponreacttli);
[sortedreactcoeff,neuord] = sort(reactcoeffuk.unsup(:,itt), 'ascend');
[~,maxvalidind] = max(sortedreactcoeff);
yt = round(1:(maxvalidind-1)/4:maxvalidind);
figure
annot = sprintf('%s UMAP+KNN decoding: trial%d spontaneous reactivations', mousedate, hireptt(itt));
annotation('textbox', [0.1 0.91 0.9 0.1], 'String', annot, 'EdgeColor', 'none');%, 'FontSize', fs)
subplot(1,2,1)
[treact, ireact] = max(Treact.duration{itt});
imagesc(squeeze(psthsponreactzrate{itt}(:,ireact,neuord))')
hold on
plot(t0ind*[1,1], 0.5+[0 Nneurons], 'k--')
caxis(2*[-1 1])
colorbar
set(gca, 'YDir', 'normal', 'YTick', yt, 'YTickLabel', sprintf('%.2f\n', sortedreactcoeff(yt)), ...
    'XTick', xt, 'XTickLabel', psthsponreacttimeline(xt))
title('longest reactivation event')
xlabel('Time (ms)')
ylabel('neuron ordered by reactivation coefficient')
subplot(1,2,2)
imagesc(squeeze(nanmean(psthsponreactzrate{itt}(:,:,neuord),2) )')
hold on
plot(t0ind*[1,1], 0.5+[0 Nneurons], 'k--')
caxis(0.5*[-1 1])
colorbar
set(gca, 'YDir', 'normal', 'YTick', yt, 'YTickLabel', sprintf('%.2f\n', sortedreactcoeff(yt)), ...
    'XTick', xt, 'XTickLabel', psthsponreacttimeline(xt))
title('average reactivation')
xlabel('Time (ms)')
colormap redblue

% neuron ordering
% order by reactivation coefficient
% order by spontaneous ensemble membership
% order by sensory responsiveness

%% order by onset of reactivation
threshzrate = 0.1;
temppsthsrzavg = squeeze(nanmean(psthsponreactzrate{itt},2));
temppsththrcross = temppsthsrzavg>threshzrate/2;
tempdtisplus = [true(1,Nneurons); diff(temppsththrcross,1,1)>0];
temppsththrfirstcross = temppsththrcross & tempdtisplus;
temppsthcum = cumsum( tempdtisplus, 1);
[tempthrcrosstime, tempthrcrossneu]= find(temppsththrfirstcross & temppsthcum==temppsthcum(t0ind,:));
thrcrosstime = NaN(Nneurons,1);
thrcrosstime(tempthrcrossneu) = tempthrcrosstime;

tempneuind = find(temppsthsrzavg(t0ind,:)>threshzrate);
[sv,si]=sort(thrcrosstime(tempneuind));
neuord = tempneuind(si);

% confirm that suprathreshold neurons have spon. react. coeff>0 and auroc>0.5

yprcts = [0 75 90 95 100];
figure; 
subplot(2,1,1); hold all
plot(reactcoeffuk.unsup(neuord,itt),'b.-')
yt = prctile(reactcoeffuk.unsup(:,itt), yprcts);
xl = [0 numel(neuord)+1];
for iy = 2:numel(yprcts)
plot(xl, yt(iy)*[1 1], '--', 'Color', yprcts(iy)/100*[1 0 0]);
text(xl(2), yt(iy), sprintf('%.0f-percentile', yprcts(iy)), 'Color', yprcts(iy)/100*[1 0 0], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
end
xlim(xl)
ylim([yt(1) yt(end)])
ylabel('react. coeff.')
annot = sprintf('%s UMAP+KNN: trial%d spon. react., neurons with z-rate>%.2f at t=0 (N=%d/%d)', mousedate, hireptt(itt), threshzrate, numel(neuord), Nneurons);
title(annot, 'fontweight', 'normal')
subplot(2,1,2); hold all
plot(sponreactauroc(neuord),'b.-'); 
tempsig = neuenssponreact.Pranksum(neuord,itt)<0.05;
plot(find(tempsig), sponreactauroc(neuord(tempsig)),'bo'); 
yt = prctile(sponreactauroc, yprcts);
xl = [0 numel(neuord)+1];
for iy = 2:numel(yprcts)
plot(xl, yt(iy)*[1 1], '--', 'Color', yprcts(iy)/100*[1 0 0]);
end
xlim(xl)
ylim([yt(1) yt(end)])
xlabel('neuron ordered by reactivation onset')
ylabel('react. auroc')

%xt = 1:20:length(psthsponreacttli);
xt = 1:4:length(psthsponreacttli);
figure
annot = sprintf('%s UMAP+KNN: trial%d spon. react., neurons with z-rate>%.2f at t=0 (N=%d/%d)', mousedate, hireptt(itt), threshzrate, numel(neuord), Nneurons);
annotation('textbox', [0 0.91 1 0.1], 'String', annot, 'EdgeColor', 'none');%, 'FontSize', fs)
imagesc(squeeze(nanmean(psthsponreactzrate{itt}(:,:,neuord),2) )')
hold on
plot(t0ind*[1,1], 0.5+[0 Nneurons], 'k--')
caxis(0.3*[-1 1])
colorbar
set(gca, 'YDir', 'reverse', 'XGrid', 'on', 'XTick', xt, 'XTickLabel', psthsponreacttimeline(xt))
title('average reactivation')
xlabel('Time (ms)')
ylabel('neurons ordered by reactivation onset')
colormap redblue
% z-rate>0.05 at -300ms could be a good cutoff for identifying pattern completion neurons... (15 neurons)
% more stringent: z-rate>0.1 at -400ms still identifies ~6 neurons

% need to check individual reactivations
figure
annot = sprintf('%s UMAP+KNN: trial%d spon. react., neurons with z-rate>%.2f at t=0 (N=%d/%d)', mousedate, hireptt(itt), threshzrate, numel(neuord), Nneurons);
annotation('textbox', [0.02 0.91 0.96 0.1], 'String', annot, 'EdgeColor', 'none');%, 'FontSize', fs)
[treactsorted, ireactsorted] = sort(Treact.duration{itt}, 'descend');
for ii = 1:8
    ireact = ireactsorted(ii);
subplot(2,4,ii)
imagesc(squeeze(psthsponreactzrate{itt}(:,ireact,neuord))')
hold on
plot(t0ind*[1,1], 0.5+[0 Nneurons], 'k--')
caxis(2*[-1 1])
colorbar
set(gca, 'YDir', 'reverse', 'XTick', xt, 'XTickLabel', psthsponreacttimeline(xt))
title(sprintf('%d-th longest reactivation event',ii))
xlabel('Time (ms)')
ylabel('neuron ordered by reactivation onset')
end
colormap redblue

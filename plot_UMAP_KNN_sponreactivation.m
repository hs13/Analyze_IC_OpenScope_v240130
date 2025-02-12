datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

testacc_dims_agg = [];
for ises = 1:Nsessions
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
load([pathpp 'UMAP_KNN_decoding_V1RS_spontaneous_embeddims.mat'])
% load([pathpp 'psth_spontaneous_spikecount_V1RS.mat'])
% load([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'])

% whichsponind=2;
% load(strcat(pathpp, 'SVMspon', num2str(whichsponind), '_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat'))

Nsplits = length(UMAPKNN_dimensions{1}.accuracy);
testacc_dims = zeros(Nsplits, numel(embedding_dimensions));
for idim = 1:numel(embedding_dimensions)
    testacc_dims(:,idim) = UMAPKNN_dimensions{idim}.accuracy;
end

testacc_dims_agg = cat(3,testacc_dims_agg,testacc_dims);
end

fs = 14;
validsesinds = 2:Nsessions;
%% embedding dimensions
testacc_dims_ses = squeeze(mean(testacc_dims_agg,1));
figure; hold all
plot(embedding_dimensions, testacc_dims_ses, '-')
errorbar(embedding_dimensions, mean(testacc_dims_ses,2), std(testacc_dims_ses,0,2)/sqrt(size(testacc_dims_ses,2)), 'ko-', 'LineWidth', 2)
set(gca, 'FontSize', fs)
xlabel('UMAP Embedding Dimensions')
ylabel('10-fold cross-validated accuracy')
title('unsupervised UMAP+KNN decoding across sessions')

[p,tbl,stats]=friedman(testacc_dims_ses(:,validsesinds)');
% [p,tbl,stats]=friedman(testacc_dims_ses');
figure; multcompare(stats)

%% RESUME EDITING HERE

hireptt = unique(trialorder);
Nhireptt = numel(hireptt);
Nsplits = size(UMAPKNN_semisup.confusion_matrix,1);

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


%% check consistency of spontaneous reactivations between UMAP settings
threshcv = 0.8; % out of Nsplits decoders, how many agreed to the same solution?
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
        % Define a threshold for outlier detection (e.g., 99% confidence level)
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

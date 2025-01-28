%% load all data
%{
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};

trialorderacc = cell(numel(nwbsessions),1);
cnt = 0;
for ises = 1:numel(nwbsessions)
    cnt = cnt+1;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    load([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'], 'trialorder')
    load([pathpp 'UMAP_V1RS_lmlvslopes.mat'])

    trialorderacc{ises} = trialorder;
    if cnt ==1
        UMAPorigagg = UMAPorig;
        UMAPorigallagg = UMAPorigall;
        UMAPorig_unsupagg = UMAPorig_unsup;
        UMAPorigall_unsupagg = UMAPorigall_unsup;
        UMAPlmlvacc = UMAPlmlv;
        UMAPlmlvallacc = UMAPlmlvall;
        UMAPlmlv_unsupacc = UMAPlmlv_unsup;
        UMAPlmlvall_unsupacc = UMAPlmlvall_unsup;
    else
        UMAPorigagg = cat(1, UMAPorigagg, UMAPorig);
        UMAPorigallagg = cat(1, UMAPorigallagg, UMAPorigall);
        UMAPorig_unsupagg = cat(1, UMAPorig_unsupagg, UMAPorig_unsup);
        UMAPorigall_unsupagg = cat(1, UMAPorigall_unsupagg, UMAPorigall_unsup);
        UMAPlmlvacc = cat(1, UMAPlmlvacc, UMAPlmlv);
        UMAPlmlvallacc = cat(1, UMAPlmlvallacc, UMAPlmlvall);
        UMAPlmlv_unsupacc = cat(1, UMAPlmlv_unsupacc, UMAPlmlv_unsup);
        UMAPlmlvall_unsupacc = cat(1, UMAPlmlvall_unsupacc, UMAPlmlvall_unsup);
    end

end

save('G:\My Drive\RESEARCH\logmean_logvar\OpenScope_UMAP_V1RS_lmlvslopes.mat', ...
    'nwbsessions', 'lmlvslope_list', 'trialorderacc', 'UMAPorigagg', 'UMAPorigallagg', 'UMAPorig_unsupagg', 'UMAPorigall_unsupagg', ...
    'UMAPlmlvacc', 'UMAPlmlvallacc', 'UMAPlmlv_unsupacc', 'UMAPlmlvall_unsupacc')
%}
%%
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])

load([drivepath 'RESEARCH/logmean_logvar/OpenScope_UMAP_V1RS_lmlvslopes.mat'])

trialorder = trialorderacc{1};
hireptt = unique(trialorder);
Nhireptt = numel(hireptt);
Ndims = 2;
Nsplits = 2;
Nsessions = numel(nwbsessions);
Nslopes = numel(lmlvslope_list);

trialorderind= trialorder;
for ii = 1:Nhireptt
    trialorderind(trialorder==hireptt(ii))=ii;
end
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

fs = 14;

isequal(UMAPlmlv_unsupacc{12, 1}.test_index, UMAPlmlv_unsupacc{12, 3}.test_index)

%% plot UMAP each session
figure;
annotation("textbox", [0.1 0.9 0.9 0.1], 'String', 'Unsupervised UMAP: all trials, original data', 'edgecolor', 'none', 'FontSize', fs)
for ises = 1:Nsessions
    subplot(3,4,ises)
    hold all
    for ii = 1:Nhireptt
        trialsoi = trialorderacc{ises}==hireptt(ii);
        tempcol = trialcol(ii,:);
        plot(UMAPorigall_unsupagg(ises).embeddings(trialsoi,1), UMAPorigall_unsupagg(ises).embeddings(trialsoi,2), '.', 'Color', tempcol)
    end
    % scatter(UMAPorigall_unsupagg(ises).embeddings(:,1), UMAPorigall_unsupagg(ises).embeddings(:,2), 2, trialorderind, 'filled')
    % colormap colorcube
    axis square
    title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
end

for slope = [0 1 2]
    islope = lmlvslope_list==slope;
    figure;
    annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Unsupervised UMAP: all trials, LMLV Slope %.1f', lmlvslope_list(islope)), 'edgecolor', 'none', 'FontSize', fs)
    for ises = 1:Nsessions
        subplot(3,4,ises)
        hold all
        for ii = 1:Nhireptt
            trialsoi = trialorderacc{ises}==hireptt(ii);
            tempcol = trialcol(ii,:);
            plot(UMAPlmlvall_unsupacc{ises,islope}.embeddings(trialsoi,1), UMAPlmlvall_unsupacc{ises,islope}.embeddings(trialsoi,2), '.', 'Color', tempcol)
        end
        axis square
        title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
    end
end

figure;
annotation("textbox", [0.1 0.9 0.9 0.1], 'String', 'Semisupervised UMAP: all trials, original data', 'edgecolor', 'none', 'FontSize', fs)
for ises = 1:Nsessions
    subplot(3,4,ises)
    hold all
    for ii = 1:Nhireptt
        trialsoi = trialorderacc{ises}==hireptt(ii);
        tempcol = trialcol(ii,:);
        plot(UMAPorigallagg(ises).embeddings(trialsoi,1), UMAPorigallagg(ises).embeddings(trialsoi,2), '.', 'Color', tempcol)
    end
    % scatter(UMAPorigallagg(ises).embeddings(:,1), UMAPorigallagg(ises).embeddings(:,2), 2, trialorderind, 'filled')
    % colormap colorcube
    axis square
    title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
end

for slope = [0 1 2]
    islope = lmlvslope_list==slope;
    figure;
    annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Semisupervised UMAP: all trials, LMLV Slope %.1f', lmlvslope_list(islope)), 'edgecolor', 'none', 'FontSize', fs)
    for ises = 1:Nsessions
        subplot(3,4,ises)
        hold all
        for ii = 1:Nhireptt
            trialsoi = trialorderacc{ises}==hireptt(ii);
            tempcol = trialcol(ii,:);
            plot(UMAPlmlvallacc{ises,islope}.embeddings(trialsoi,1), UMAPlmlvallacc{ises,islope}.embeddings(trialsoi,2), '.', 'Color', tempcol)
        end
        axis square
        title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
    end
end

%% cross-validated test trials
isplit = 1;
figure;
annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Unsupervised UMAP: %d/2-fold cross-validated test trials, original data', isplit), 'edgecolor', 'none', 'FontSize', fs)
for ises = 1:Nsessions
    subplot(3,4,ises)
    hold all
    for ii = 1:Nhireptt
        trialsoi = UMAPorig_unsupagg(ises).test_truelabels(isplit,:)==hireptt(ii);
        tempcol = trialcol(ii,:);
        plot(UMAPorig_unsupagg(ises).test_embeddings(isplit,trialsoi,1), UMAPorig_unsupagg(ises).test_embeddings(isplit,trialsoi,2), '.', 'Color', tempcol)
    end
    axis square
    title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
end

for slope = [0 1 2]
    islope = lmlvslope_list==slope;
    figure;
    annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Unsupervised UMAP: %d/2-fold cross-validated test trials, LMLV Slope %.1f', isplit, lmlvslope_list(islope)), 'edgecolor', 'none', 'FontSize', fs)
    for ises = 1:Nsessions
        subplot(3,4,ises)
        hold all
        for ii = 1:Nhireptt
            trialsoi = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(isplit,:)==hireptt(ii);
            tempcol = trialcol(ii,:);
            plot(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,trialsoi,1), UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,trialsoi,2), '.', 'Color', tempcol)
        end
        axis square
        title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
    end
end


figure;
annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Semisupervised UMAP: %d/2-fold cross-validated test trials, original data', isplit), 'edgecolor', 'none', 'FontSize', fs)
for ises = 1:Nsessions
    subplot(3,4,ises)
    hold all
    for ii = 1:Nhireptt
        trialsoi = UMAPorigagg(ises).test_truelabels(isplit,:)==hireptt(ii);
        tempcol = trialcol(ii,:);
        plot(UMAPorigagg(ises).test_embeddings(isplit,trialsoi,1), UMAPorigagg(ises).test_embeddings(isplit,trialsoi,2), '.', 'Color', tempcol)
    end
    axis square
    title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
end

for slope = [0 1 2]
    islope = lmlvslope_list==slope;
    figure;
    annotation("textbox", [0.1 0.9 0.9 0.1], 'String', sprintf('Semisupervised UMAP: %d/2-fold cross-validated test trials, LMLV Slope %.1f', isplit, lmlvslope_list(islope)), 'edgecolor', 'none', 'FontSize', fs)
    for ises = 1:Nsessions
        subplot(3,4,ises)
        hold all
        for ii = 1:Nhireptt
            trialsoi = UMAPlmlvacc{ises,islope}.test_truelabels(isplit,:)==hireptt(ii);
            tempcol = trialcol(ii,:);
            plot(UMAPlmlvacc{ises,islope}.test_embeddings(isplit,trialsoi,1), UMAPlmlvacc{ises,islope}.test_embeddings(isplit,trialsoi,2), '.', 'Color', tempcol)
        end
        axis square
        title(sprintf('%d Session %s', ises, nwbsessions{ises}), 'FontSize', fs)
    end
end

%% test accuracy across LMLV slopes
testacc_semisup_lmlv = NaN(numel(lmlvslope_list), numel(nwbsessions), Nsplits);
testacc_unsup_lmlv = NaN(numel(lmlvslope_list), numel(nwbsessions), Nsplits);
for ises = 1:numel(nwbsessions)
    for islope = 1:numel(lmlvslope_list)
        testacc_semisup_lmlv(islope,ises,:) = UMAPlmlvacc{ises,islope}.accuracy;
        testacc_unsup_lmlv(islope,ises,:) = UMAPlmlv_unsupacc{ises,islope}.accuracy;
    end
end

fs = 18;
xl = [lmlvslope_list(1) lmlvslope_list(end)];
chanceperf = 1/Nhireptt;
figure('Position', [100 100 1200 400])
subplot(1,3,1)
hold all
plot(lmlvslope_list, squeeze(mean(testacc_semisup_lmlv,3)) )
plot(lmlvslope_list, squeeze(mean(testacc_semisup_lmlv,[2,3])), 'k-', 'LineWidth', 2 )
plot(xl, chanceperf*[1 1], '--', 'Color', 0.5+[0 0 0])
ylim([0 1])
set(gca, 'FontSize', fs)
title('semi-supervised', 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
ylabel('test accuracy', 'FontSize', fs)
subplot(1,3,2)
hold all
plot(lmlvslope_list, squeeze(mean(testacc_unsup_lmlv,3)) )
plot(lmlvslope_list, squeeze(mean(testacc_unsup_lmlv,[2,3])), 'k-', 'LineWidth', 2 )
plot(xl, chanceperf*[1 1], '--', 'Color', 0.5+[0 0 0])
ylim([0 1])
set(gca, 'FontSize', fs)
title('unsupervised', 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
ylabel('test accuracy', 'FontSize', fs)
subplot(1,3,3)
hold all
errorbar(lmlvslope_list, squeeze(mean(testacc_semisup_lmlv,[2,3])), squeeze(std(mean(testacc_semisup_lmlv,3),0,2))/sqrt(numel(nwbsessions)), 'o-', 'LineWidth', 2 )
errorbar(lmlvslope_list, squeeze(mean(testacc_unsup_lmlv,[2,3])), squeeze(std(mean(testacc_unsup_lmlv,3),0,2))/sqrt(numel(nwbsessions)), 'x-', 'LineWidth', 2 )
plot(xl, chanceperf*[1 1], '--', 'Color', 0.5+[0 0 0])
legend({'semi-supervised', 'unsupervised'}, 'location', 'northeast', 'FontSize', fs)
set(gca, 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
ylabel('test accuracy', 'FontSize', fs)
title('ICwcfg1 hireptt UMAP+KNN-decoding', 'FontSize', fs)

%% centroid distance on UMAP between pairs of trial types

umapcentroid = struct();
umapcentroid.unsup_origall = NaN(Nhireptt, Ndims, numel(nwbsessions));
umapcentroid.unsup_origtest = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.unsup_origtrain = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.semisup_origall = NaN(Nhireptt, Ndims, numel(nwbsessions));
umapcentroid.semisup_origtest = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.semisup_origtrain = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);

umapcentroid.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umapcentroid.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umapcentroid.unsup_lmlvtrain = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvtrain = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umapcentroid.unsup_lmlvall{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions));
    umapcentroid.unsup_lmlvtest{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.unsup_lmlvtrain{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.semisup_lmlvall{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions));
    umapcentroid.semisup_lmlvtest{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.semisup_lmlvtrain{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
end

umapctrdist = struct();
umapctrdist.unsup_origall = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
umapctrdist.unsup_origtest = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.unsup_origtrain = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.semisup_origall = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
umapctrdist.semisup_origtest = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.semisup_origtrain = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);

umapctrdist.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umapctrdist.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umapctrdist.unsup_lmlvtrain = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvtrain = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umapctrdist.unsup_lmlvall{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
    umapctrdist.unsup_lmlvtest{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.unsup_lmlvtrain{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.semisup_lmlvall{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
    umapctrdist.semisup_lmlvtest{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.semisup_lmlvtrain{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
end

for ises = 1:numel(nwbsessions)
    tic
    for ieuc = 1:2
        switch ieuc
            case 1
                eucfield = 'semisup_origall';
                tempembedding = UMAPorigallagg(ises).embeddings;
                temptrialorder = trialorderacc{ises};
            case 2
                eucfield = 'unsup_origall';
                tempembedding = UMAPorigall_unsupagg(ises).embeddings;
                temptrialorder = trialorderacc{ises};
        end
        tempcentroid = NaN(Nhireptt, Ndims);
        for ii = 1:Nhireptt
            trialsoi = temptrialorder==hireptt(ii);
            tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
        end
        umapcentroid.(eucfield)(:,:,ises) = tempcentroid;
        
        tempctrdist = zeros(Nhireptt,Nhireptt);
        for idim = 1:Ndims
            tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
        end
        tempctrdist = sqrt(tempctrdist);
        umapctrdist.(eucfield)(:,:,ises) = tempctrdist;
        
    end
    
    for ieuc = 1:4
        for isplit = 1:Nsplits
            switch ieuc
                case 1
                    eucfield = 'semisup_origtest';
                    tempembedding = squeeze(UMAPorigagg(ises).test_embeddings(isplit,:,:));
                    temptrialorder = UMAPorigagg(ises).test_truelabels(isplit,:);
                case 2
                    eucfield = 'semisup_origtrain';
                    tempembedding = squeeze(UMAPorigagg(ises).train_embeddings(isplit,:,:));
                    temptrialorder = trialorderacc{ises}( 1+UMAPorigagg(ises).train_index(isplit,:) );
                case 3
                    eucfield = 'unsup_origtest';
                    tempembedding = squeeze(UMAPorig_unsupagg(ises).test_embeddings(isplit,:,:));
                    temptrialorder = UMAPorig_unsupagg(ises).test_truelabels(isplit,:);
                case 4
                    eucfield = 'unsup_origtrain';
                    tempembedding = squeeze(UMAPorig_unsupagg(ises).train_embeddings(isplit,:,:));
                    temptrialorder = trialorderacc{ises}( 1+UMAPorig_unsupagg(ises).train_index(isplit,:) );
            end
            tempcentroid = NaN(Nhireptt, Ndims);
            for ii = 1:Nhireptt
                trialsoi = temptrialorder==hireptt(ii);
                tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
            end
            umapcentroid.(eucfield)(:,:,ises,isplit) = tempcentroid;
            
            tempctrdist = zeros(Nhireptt,Nhireptt);
            for idim = 1:Ndims
                tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
            end
            tempctrdist = sqrt(tempctrdist);
            umapctrdist.(eucfield)(:,:,ises,isplit) = tempctrdist;
        end
    end
    
    for islope = 1:numel(lmlvslope_list)
        for ieuc = 1:2
            switch ieuc
                case 1
                    eucfield = 'semisup_lmlvall';
                    tempembedding = UMAPlmlvallacc{ises,islope}.embeddings;
                    temptrialorder = trialorderacc{ises};
                case 2
                    eucfield = 'unsup_lmlvall';
                    tempembedding = UMAPlmlvall_unsupacc{ises,islope}.embeddings;
                    temptrialorder = trialorderacc{ises};
            end
            tempcentroid = NaN(Nhireptt, Ndims);
            for ii = 1:Nhireptt
                trialsoi = temptrialorder==hireptt(ii);
                tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
            end
            umapcentroid.(eucfield){islope}(:,:,ises) = tempcentroid;
            
            tempctrdist = zeros(Nhireptt,Nhireptt);
            for idim = 1:Ndims
                tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
            end
            tempctrdist = sqrt(tempctrdist);
            umapctrdist.(eucfield){islope}(:,:,ises) = tempctrdist;
            
        end
        
        for ieuc = 1:4
            for isplit = 1:Nsplits
                switch ieuc
                    case 1
                        eucfield = 'semisup_lmlvtest';
                        tempembedding = squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(isplit,:,:));
                        temptrialorder = UMAPlmlvacc{ises,islope}.test_truelabels(isplit,:);
                    case 2
                        eucfield = 'semisup_lmlvtrain';
                        tempembedding = squeeze(UMAPlmlvacc{ises,islope}.train_embeddings(isplit,:,:));
                        temptrialorder = trialorderacc{ises}( 1+UMAPlmlvacc{ises,islope}.train_index(isplit,:) );
                    case 3
                        eucfield = 'unsup_lmlvtest';
                        tempembedding = squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,:,:));
                        temptrialorder = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(isplit,:);
                    case 4
                        eucfield = 'unsup_lmlvtrain';
                        tempembedding = squeeze(UMAPlmlv_unsupacc{ises,islope}.train_embeddings(isplit,:,:));
                        temptrialorder = trialorderacc{ises}( 1+UMAPlmlv_unsupacc{ises,islope}.train_index(isplit,:) );
                end
                tempcentroid = NaN(Nhireptt, Ndims);
                for ii = 1:Nhireptt
                    trialsoi = temptrialorder==hireptt(ii);
                    tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
                end
                umapcentroid.(eucfield){islope}(:,:,ises,isplit) = tempcentroid;
                
                tempctrdist = zeros(Nhireptt,Nhireptt);
                for idim = 1:Ndims
                    tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
                end
                tempctrdist = sqrt(tempctrdist);
                umapctrdist.(eucfield){islope}(:,:,ises,isplit) = tempctrdist;
                
            end
        end
    end
    toc
end

%%  between cross-validation sets: semi-supervised vs unsupervised
corrtype = 'Pearson';
triuttpairind = triu(true(Nhireptt),1);

% between different trials: cross-validation test sets
corrumapctrdist_unsup_origtest = NaN(numel(nwbsessions),1);
corrumapctrdist_semisup_origtest = NaN(numel(nwbsessions),1);
cossimumapctrdist_unsup_origtest = NaN(numel(nwbsessions),1);
cossimumapctrdist_semisup_origtest = NaN(numel(nwbsessions),1);
for ises = 1:numel(nwbsessions)
    vec1 = umapctrdist.unsup_origtest(:,:,ises,1);
    vec2 = umapctrdist.unsup_origtest(:,:,ises,2);
    corrumapctrdist_unsup_origtest(ises) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
    cossimumapctrdist_unsup_origtest(ises) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
    
    vec1 = umapctrdist.semisup_origtest(:,:,ises,1);
    vec2 = umapctrdist.semisup_origtest(:,:,ises,2);
    corrumapctrdist_semisup_origtest(ises) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
    cossimumapctrdist_semisup_origtest(ises) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
end

corrumapctrdist_unsup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
corrumapctrdist_semisup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_unsup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_semisup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for islope = 1:numel(lmlvslope_list)
        vec1 = umapctrdist.unsup_lmlvtest{islope}(:,:,ises,1);
        vec2 = umapctrdist.unsup_lmlvtest{islope}(:,:,ises,2);
        corrumapctrdist_unsup_lmlvtest(ises,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_unsup_lmlvtest(ises,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
        
        vec1 = umapctrdist.semisup_lmlvtest{islope}(:,:,ises,1);
        vec2 = umapctrdist.semisup_lmlvtest{islope}(:,:,ises,2);
        corrumapctrdist_semisup_lmlvtest(ises,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_semisup_lmlvtest(ises,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
    end
end

% between different sessions: cross-validated test trials
procrustesumapctr_unsup_origtest = NaN(numel(nwbsessions),1);
procrustesumapctr_semisup_origtest = NaN(numel(nwbsessions),1);
for ises = 1:numel(nwbsessions)
    vec1 = squeeze(umapcentroid.unsup_origtest(:,:,ises,1));
    vec2 = squeeze(umapcentroid.unsup_origtest(:,:,ises,2));
    procrustesumapctr_unsup_origtest(ises) = procrustes(vec1, vec2);
    
    vec1 = squeeze(umapcentroid.semisup_origtest(:,:,ises,1));
    vec2 = squeeze(umapcentroid.semisup_origtest(:,:,ises,2));
    procrustesumapctr_semisup_origtest(ises) = procrustes(vec1, vec2);
end

procrustesumapctr_unsup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
procrustesumapctr_semisup_lmlvtest = NaN(numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for islope = 1:numel(lmlvslope_list)
        vec1 = squeeze(umapcentroid.unsup_lmlvtest{islope}(:,:,ises,1));
        vec2 = squeeze(umapcentroid.unsup_lmlvtest{islope}(:,:,ises,2));
        procrustesumapctr_unsup_lmlvtest(ises,islope) = procrustes(vec1, vec2);
        
        vec1 = squeeze(umapcentroid.semisup_lmlvtest{islope}(:,:,ises,1));
        vec2 = squeeze(umapcentroid.semisup_lmlvtest{islope}(:,:,ises,2));
        procrustesumapctr_semisup_lmlvtest(ises,islope) = procrustes(vec1, vec2);
    end
end


fs = 14;
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox', [0 0.91 0.9 0.1], 'string', 'UMAP centroid per trial type: compare between test sets', 'edgecolor', 'none', 'FontSize', fs)
subplot(2,2,1)
hold all
errorbar(lmlvslope_list, squeeze(mean(corrumapctrdist_unsup_lmlvtest,1)), squeeze(std(corrumapctrdist_unsup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'bo-', 'LineWidth', 2 )
errorbar(lmlvslope_list, squeeze(mean(corrumapctrdist_semisup_lmlvtest,1)), squeeze(std(corrumapctrdist_semisup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'rx-', 'LineWidth', 2 )
plot(xl, squeeze(mean(corrumapctrdist_unsup_origtest,1))*[1 1], 'b-', 'LineWidth', 1 )
plot(xl, squeeze(mean(corrumapctrdist_semisup_origtest,1))*[1 1], 'r-', 'LineWidth', 1 )
yl = ylim; plot([1 1], yl, 'k--'); ylim(yl)
set(gca, 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
title(sprintf('%s correlation btwn ctr dist.', corrtype), 'FontSize', fs)
subplot(2,2,2)
hold all
errorbar(lmlvslope_list, squeeze(mean(cossimumapctrdist_unsup_lmlvtest,1)), squeeze(std(cossimumapctrdist_unsup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'bo-', 'LineWidth', 2 )
errorbar(lmlvslope_list, squeeze(mean(cossimumapctrdist_semisup_lmlvtest,1)), squeeze(std(cossimumapctrdist_semisup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'rx-', 'LineWidth', 2 )
plot(xl, squeeze(mean(cossimumapctrdist_unsup_origtest,1))*[1 1], 'b-', 'LineWidth', 1 )
plot(xl, squeeze(mean(cossimumapctrdist_semisup_origtest,1))*[1 1], 'r-', 'LineWidth', 1 )
yl = ylim; plot([1 1], yl, 'k--'); ylim(yl)
set(gca, 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
title('cosine similarity btwn ctr dist.', 'FontSize', fs)

subplot(2,2,3)
hold all
errorbar(lmlvslope_list, squeeze(mean(procrustesumapctr_unsup_lmlvtest,1)), squeeze(std(procrustesumapctr_unsup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'bo-', 'LineWidth', 2 )
errorbar(lmlvslope_list, squeeze(mean(procrustesumapctr_semisup_lmlvtest,1)), squeeze(std(procrustesumapctr_semisup_lmlvtest,0,1))/sqrt(numel(nwbsessions)), 'rx-', 'LineWidth', 2 )
plot(xl, squeeze(mean(procrustesumapctr_unsup_origtest,1))*[1 1], 'b-', 'LineWidth', 1 )
plot(xl, squeeze(mean(procrustesumapctr_semisup_origtest,1))*[1 1], 'r-', 'LineWidth', 1 )
legend({'unsupervised', 'semi-supervised'}, 'location', 'northeast')
yl = ylim; plot([1 1], yl, 'k--'); ylim(yl)
set(gca, 'FontSize', fs)
xlabel('LMLV slopes', 'FontSize', fs)
title('Procrustes distance', 'FontSize', fs)

%% cosine similarity between centroid distance matrices
% todo: consider adding jensen-shannon divergence between centroid distance matrices
corrtype = 'Pearson';
splitoi = 1;
triuttpairind = triu(true(Nhireptt),1);

% between different sessions: cross-validated test trials
corrumapctrdist_unsup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
corrumapctrdist_semisup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
cossimumapctrdist_unsup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
cossimumapctrdist_semisup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        vec1 = mean(umapctrdist.unsup_origtest(:,:,ises,splitoi),4);
        vec2 = mean(umapctrdist.unsup_origtest(:,:,jses,splitoi),4);
        corrumapctrdist_unsup_origsespair(ises,jses) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_unsup_origsespair(ises,jses) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
        
        vec1 = mean(umapctrdist.semisup_origtest(:,:,ises,splitoi),4);
        vec2 = mean(umapctrdist.semisup_origtest(:,:,jses,splitoi),4);
        corrumapctrdist_semisup_origsespair(ises,jses) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_semisup_origsespair(ises,jses) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
    end
end

corrumapctrdist_unsup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
corrumapctrdist_semisup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_unsup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_semisup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        for islope = 1:numel(lmlvslope_list)
            vec1 = mean(umapctrdist.unsup_lmlvtest{islope}(:,:,ises,splitoi),4);
            vec2 = mean(umapctrdist.unsup_lmlvtest{islope}(:,:,jses,splitoi),4);
            corrumapctrdist_unsup_lmlvsespair(ises,jses,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
            cossimumapctrdist_unsup_lmlvsespair(ises,jses,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
            
            vec1 = mean(umapctrdist.semisup_lmlvtest{islope}(:,:,ises,splitoi),4);
            vec2 = mean(umapctrdist.semisup_lmlvtest{islope}(:,:,jses,splitoi),4);
            corrumapctrdist_semisup_lmlvsespair(ises,jses,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
            cossimumapctrdist_semisup_lmlvsespair(ises,jses,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
        end
    end
end

% between different sessions: UMAP on all trials (not cross-validated)
corrumapctrdist_unsup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
corrumapctrdist_semisup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
cossimumapctrdist_unsup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
cossimumapctrdist_semisup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        vec1 = umapctrdist.unsup_origall(:,:,ises);
        vec2 = umapctrdist.unsup_origall(:,:,jses);
        corrumapctrdist_unsup_origallsespair(ises,jses) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_unsup_origallsespair(ises,jses) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
        
        vec1 = umapctrdist.semisup_origall(:,:,ises);
        vec2 = umapctrdist.semisup_origall(:,:,jses);
        corrumapctrdist_semisup_origallsespair(ises,jses) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
        cossimumapctrdist_semisup_origallsespair(ises,jses) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
    end
end

corrumapctrdist_unsup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
corrumapctrdist_semisup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_unsup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
cossimumapctrdist_semisup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        for islope = 1:numel(lmlvslope_list)
            vec1 = umapctrdist.unsup_lmlvall{islope}(:,:,ises);
            vec2 = umapctrdist.unsup_lmlvall{islope}(:,:,jses);
            corrumapctrdist_unsup_lmlvallsespair(ises,jses,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
            cossimumapctrdist_unsup_lmlvallsespair(ises,jses,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
            
            vec1 = umapctrdist.semisup_lmlvall{islope}(:,:,ises);
            vec2 = umapctrdist.semisup_lmlvall{islope}(:,:,jses);
            corrumapctrdist_semisup_lmlvallsespair(ises,jses,islope) = corr(vec1(triuttpairind), vec2(triuttpairind), 'type', corrtype);
            cossimumapctrdist_semisup_lmlvallsespair(ises,jses,islope) = (vec1(triuttpairind)'*vec2(triuttpairind))/(norm(vec1(triuttpairind))*norm(vec2(triuttpairind)));
        end
    end
end



triusespairind = triu(true(numel(nwbsessions)),1);
cossimumapdistmat = cat(2, cossimumapctrdist_unsup_origallsespair(triusespairind), ...
    cossimumapctrdist_semisup_origallsespair(triusespairind), ...
    cossimumapctrdist_unsup_origsespair(triusespairind), ...
    cossimumapctrdist_semisup_origsespair(triusespairind) );

figure
hold all
plot(1:4, cossimumapdistmat, '-')
errorbar(1:4, mean(cossimumapdistmat,1), std(cossimumapdistmat,0,1)/sqrt(size(cossimumapdistmat,1)), 'ko-', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'all-unsup', 'all-semi', 'cv-unsup', 'cv-semi'})

[p,tbl,stats]=friedman(cossimumapdistmat);
figure; multcompare(stats)

figure
for isp = 1:4
    switch isp
        case 1
            cossimumapdistorig = cossimumapctrdist_unsup_origallsespair;
            cossimumapdistlmlv = cossimumapctrdist_unsup_lmlvallsespair;
            sptitle = 'UMAP all-trials, unsupervised';
        case 2
            cossimumapdistorig = cossimumapctrdist_semisup_origallsespair;
            cossimumapdistlmlv = cossimumapctrdist_semisup_lmlvallsespair;
            sptitle = 'UMAP all-trials, semi-supervised';
        case 3
            cossimumapdistorig = cossimumapctrdist_unsup_origsespair;
            cossimumapdistlmlv = cossimumapctrdist_unsup_lmlvsespair;
            sptitle = 'UMAP cross-validated, unsupervised';
        case 4
            cossimumapdistorig = cossimumapctrdist_semisup_origsespair;
            cossimumapdistlmlv = cossimumapctrdist_semisup_lmlvsespair;
            sptitle = 'UMAP cross-validated, semi-supervised';
    end
    cossimumapdistmat = NaN(nnz(triusespairind), numel(lmlvslope_list));
    for islope = 1:numel(lmlvslope_list)
        tempmat = cossimumapdistlmlv(:,:,islope);
        cossimumapdistmat(:,islope) = tempmat(triusespairind);
    end
    
    subplot(2,2,isp)
    hold all
    plot([lmlvslope_list(1) lmlvslope_list(end)], mean(cossimumapdistorig(triusespairind))*[1 1], 'c-', 'LineWidth', 1)
    errorbar(lmlvslope_list, mean(cossimumapdistmat,1), std(cossimumapdistmat,0,1)/sqrt(nnz(triusespairind)), 'ko-', 'LineWidth', 1)
    title(sptitle)
    ylabel('UMAP ctr dist. cos. sim.')
end


triusespairind = triu(true(numel(nwbsessions)),1);
corrumapdistmat = cat(2, corrumapctrdist_unsup_origallsespair(triusespairind), ...
    corrumapctrdist_semisup_origallsespair(triusespairind), ...
    corrumapctrdist_unsup_origsespair(triusespairind), ...
    corrumapctrdist_semisup_origsespair(triusespairind) );

figure
hold all
plot(1:4, corrumapdistmat, '-')
errorbar(1:4, mean(corrumapdistmat,1), std(corrumapdistmat,0,1)/sqrt(size(corrumapdistmat,1)), 'ko-', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'all-unsup', 'all-semi', 'cv-unsup', 'cv-semi'})

[p,tbl,stats]=friedman(corrumapdistmat);
figure; multcompare(stats)

figure
for isp = 1:4
    switch isp
        case 1
            corrumapdistorig = corrumapctrdist_unsup_origallsespair;
            corrumapdistlmlv = corrumapctrdist_unsup_lmlvallsespair;
            sptitle = 'UMAP all-trials, unsupervised';
        case 2
            corrumapdistorig = corrumapctrdist_semisup_origallsespair;
            corrumapdistlmlv = corrumapctrdist_semisup_lmlvallsespair;
            sptitle = 'UMAP all-trials, semi-supervised';
        case 3
            corrumapdistorig = corrumapctrdist_unsup_origsespair;
            corrumapdistlmlv = corrumapctrdist_unsup_lmlvsespair;
            sptitle = 'UMAP cross-validated, unsupervised';
        case 4
            corrumapdistorig = corrumapctrdist_semisup_origsespair;
            corrumapdistlmlv = corrumapctrdist_semisup_lmlvsespair;
            sptitle = 'UMAP cross-validated, semi-supervised';
    end
    corrumapdistmat = NaN(nnz(triusespairind), numel(lmlvslope_list));
    for islope = 1:numel(lmlvslope_list)
        tempmat = corrumapdistlmlv(:,:,islope);
        corrumapdistmat(:,islope) = tempmat(triusespairind);
    end
    
    subplot(2,2,isp)
    hold all
    plot([lmlvslope_list(1) lmlvslope_list(end)], mean(corrumapdistorig(triusespairind))*[1 1], 'c-', 'LineWidth', 1)
    errorbar(lmlvslope_list, mean(corrumapdistmat,1), std(corrumapdistmat,0,1)/sqrt(nnz(triusespairind)), 'ko-', 'LineWidth', 1)
    title(sptitle)
    ylabel(sprintf('UMAP ctr dist. R%s', corrtype))
end


%% procrustes distance between UMAP centroids
isplit=1;

% between different sessions: cross-validated test trials
procrustesumapctr_unsup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
procrustesumapctr_semisup_origsespair = NaN(numel(nwbsessions),numel(nwbsessions));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        vec1 = squeeze(umapcentroid.unsup_origtest(:,:,ises,isplit));
        vec2 = squeeze(umapcentroid.unsup_origtest(:,:,jses,isplit));
        procrustesumapctr_unsup_origsespair(ises,jses) = procrustes(vec1, vec2);
        
        vec1 = squeeze(umapcentroid.semisup_origtest(:,:,ises,isplit));
        vec2 = squeeze(umapcentroid.semisup_origtest(:,:,jses,isplit));
        procrustesumapctr_semisup_origsespair(ises,jses) = procrustes(vec1, vec2);
    end
end

procrustesumapctr_unsup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
procrustesumapctr_semisup_lmlvsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        for islope = 1:numel(lmlvslope_list)
            vec1 = squeeze(umapcentroid.unsup_lmlvtest{islope}(:,:,ises,isplit));
            vec2 = squeeze(umapcentroid.unsup_lmlvtest{islope}(:,:,jses,isplit));
            procrustesumapctr_unsup_lmlvsespair(ises,jses,islope) = procrustes(vec1, vec2);
            
            vec1 = squeeze(umapcentroid.semisup_lmlvtest{islope}(:,:,ises,isplit));
            vec2 = squeeze(umapcentroid.semisup_lmlvtest{islope}(:,:,jses,isplit));
            procrustesumapctr_semisup_lmlvsespair(ises,jses,islope) = procrustes(vec1, vec2);
        end
    end
end

% between different sessions: UMAP on all trials (not cross-validated)
procrustesumapctr_unsup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
procrustesumapctr_semisup_origallsespair = NaN(numel(nwbsessions),numel(nwbsessions));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        vec1 = umapcentroid.unsup_origall(:,:,ises);
        vec2 = umapcentroid.unsup_origall(:,:,jses);
        procrustesumapctr_unsup_origallsespair(ises,jses) = procrustes(vec1, vec2);
        
        vec1 = umapcentroid.semisup_origall(:,:,ises);
        vec2 = umapcentroid.semisup_origall(:,:,jses);
        procrustesumapctr_semisup_origallsespair(ises,jses) = procrustes(vec1, vec2);
    end
end

procrustesumapctr_unsup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
procrustesumapctr_semisup_lmlvallsespair = NaN(numel(nwbsessions),numel(nwbsessions),numel(lmlvslope_list));
for ises = 1:numel(nwbsessions)
    for jses = ises+1:numel(nwbsessions)
        for islope = 1:numel(lmlvslope_list)
            vec1 = umapcentroid.unsup_lmlvall{islope}(:,:,ises);
            vec2 = umapcentroid.unsup_lmlvall{islope}(:,:,jses);
            procrustesumapctr_unsup_lmlvallsespair(ises,jses,islope) = procrustes(vec1, vec2);
            
            vec1 = umapcentroid.semisup_lmlvall{islope}(:,:,ises);
            vec2 = umapcentroid.semisup_lmlvall{islope}(:,:,jses);
            procrustesumapctr_semisup_lmlvallsespair(ises,jses,islope) = procrustes(vec1, vec2);
        end
    end
end



triusespairind = triu(true(numel(nwbsessions)),1);
procrustesumapctrmat = cat(2, procrustesumapctr_unsup_origallsespair(triusespairind), ...
    procrustesumapctr_semisup_origallsespair(triusespairind), ...
    procrustesumapctr_unsup_origsespair(triusespairind), ...
    procrustesumapctr_semisup_origsespair(triusespairind) );

figure
hold all
plot(1:4, procrustesumapctrmat, '-')
errorbar(1:4, mean(procrustesumapctrmat,1), std(procrustesumapctrmat,0,1)/sqrt(size(procrustesumapctrmat,1)), 'ko-', 'LineWidth', 2)
set(gca, 'XTick', 1:4, 'XTickLabel', {'all-unsup', 'all-semi', 'cv-unsup', 'cv-semi'})

[p,tbl,stats]=friedman(procrustesumapctrmat);
figure; multcompare(stats)

figure
for isp = 1:4
    switch isp
        case 1
            procrustesumapctrorig = procrustesumapctr_unsup_origallsespair;
            procrustesumapctrlmlv = procrustesumapctr_unsup_lmlvallsespair;
            sptitle = 'UMAP all-trials, unsupervised';
        case 2
            procrustesumapctrorig = procrustesumapctr_semisup_origallsespair;
            procrustesumapctrlmlv = procrustesumapctr_semisup_lmlvallsespair;
            sptitle = 'UMAP all-trials, semi-supervised';
        case 3
            procrustesumapctrorig = procrustesumapctr_unsup_origsespair;
            procrustesumapctrlmlv = procrustesumapctr_unsup_lmlvsespair;
            sptitle = 'UMAP cross-validated, unsupervised';
        case 4
            procrustesumapctrorig = procrustesumapctr_semisup_origsespair;
            procrustesumapctrlmlv = procrustesumapctr_semisup_lmlvsespair;
            sptitle = 'UMAP cross-validated, semi-supervised';
    end
    procrustesumapctrmat = NaN(nnz(triusespairind), numel(lmlvslope_list));
    for islope = 1:numel(lmlvslope_list)
        tempmat = procrustesumapctrlmlv(:,:,islope);
        procrustesumapctrmat(:,islope) = tempmat(triusespairind);
    end
    
    subplot(2,2,isp)
    hold all
    plot([lmlvslope_list(1) lmlvslope_list(end)], mean(procrustesumapctrorig(triusespairind))*[1 1], 'c-', 'LineWidth', 1)
    plot(lmlvslope_list, procrustesumapctrmat)
    errorbar(lmlvslope_list, mean(procrustesumapctrmat,1), std(procrustesumapctrmat,0,1)/sqrt(nnz(triusespairind)), 'ko-', 'LineWidth', 1)
    title(sptitle)
    % ylim([0.6 0.85])
    ylabel('UMAP centroid procrustes')
end

%% visualize UMAP centroid Procrustes transformation
% transform between semisupervised all trials sessions 2 and 9
ises = 2;
jses = 9;

UMAPopt = 'unsup_origall';

switch UMAPopt
    case 'semisup_origall'
        umap1 = double( UMAPorigallagg(ises).embeddings );
        umap2 = double( UMAPorigallagg(jses).embeddings );
        
        umapctr1 = umapcentroid.semisup_origall(:,:,ises);
        umapctr2 = umapcentroid.semisup_origall(:,:,jses);
        
        trialord1 = trialorderacc{ises};
        trialord2 = trialorderacc{jses};
    case 'unsup_origall'
        umap1 = double( UMAPorigall_unsupagg(ises).embeddings );
        umap2 = double( UMAPorigall_unsupagg(jses).embeddings );
        
        umapctr1 = umapcentroid.unsup_origall(:,:,ises);
        umapctr2 = umapcentroid.unsup_origall(:,:,jses);
        
        trialord1 = trialorderacc{ises};
        trialord2 = trialorderacc{jses};
end

figure
hold all
scatter(umapctr1(:,1), umapctr1(:,2), 250,1:Nhireptt, 'o', 'linewidth',1)
scatter(Z(:,1), Z(:,2), 250,1:Nhireptt, 'x', 'linewidth',1)
colormap colorcube

figure
hold all
for itt = 1:Nhireptt
    trialsoj = find( trialord2==hireptt(itt) );
    plot(umap2tf(trialsoj,1), umap2tf(trialsoj,2), '.', 'Color', trialcol(itt,:))
    scatter(Z(itt,1), Z(itt,2),100,'d', 'MarkerFaceColor', trialcol(itt,:), 'MarkerEdgeColor', [0.5 0.5 0.5], 'linewidth',2)
end

figure
for isp = 1:3
    subplot(2,2,isp)
    hold all
    for itt = 1:Nhireptt
        switch isp
            case 1
                trialsoi = find( trialorder1==hireptt(itt) );
                P = umap1(trialsoi,:);
                sptitle = sprintf('UMAP Ses%d', ises);
            case 2
                trialsoi = find( trialord2==hireptt(itt) );
                P = umap2(trialsoi,:);
                sptitle = sprintf('UMAP Ses%d', jses);
            case 3
                trialsoi = find( trialord2==hireptt(itt) );
                P = umap2tf(trialsoi,:);
                sptitle = sprintf('UMAP-transformed Ses%d', jses);
        end
        distances = mahal(P, P);
        % Define a threshold for outlier detection (e.g., 99% confidence level)
        threshold = chi2inv(0.95, size(P, 2)); % Chi-square critical value
        outliers = distances > threshold;
        % Remove outliers
        P_cleaned = P(~outliers, :);
        [k,v]=boundary(P_cleaned);
        
        plot(P(:,1),P(:,2), '.', 'Color', trialcol(itt,:))
        %     plot(P_cleaned(k,1),P_cleaned(k,2), 'FaceColor', [trialcol(itt,:) 0.5])
        fill(P_cleaned(k,1),P_cleaned(k,2), trialcol(itt,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        title(sptitle)
        
        %     [k0,v0]=boundary(P_cleaned,0);
        %     [k1,v1]=boundary(P_cleaned,1);
        %     plot(P_cleaned(k0,1),P_cleaned(k0,2), 'r-')
        %     plot(P_cleaned(k1,1),P_cleaned(k1,2), 'b-')
    end
end
subplot(2,2,4)
hold all
for ses2plt = 1:2
    for itt = 1:Nhireptt
        switch ses2plt
            case 1
                trialsoi = find( trialorderacc{ises}==hireptt(itt) );
                P = umap1(trialsoi,:);
                sptitle = sprintf('UMAP Ses%d', ises);
                ecol = [1 0 0];
            case 2
                trialsoi = find( trialord2==hireptt(itt) );
                P = umap2tf(trialsoi,:);
                sptitle = sprintf('UMAP-transformed Ses%d', jses);
                ecol = [0 0 1];
        end
        distances = mahal(P, P);
        % Define a threshold for outlier detection (e.g., 99% confidence level)
        threshold = chi2inv(0.95, size(P, 2)); % Chi-square critical value
        outliers = distances > threshold;
        % Remove outliers
        P_cleaned = P(~outliers, :);
        [k,v]=boundary(P_cleaned);
        
        %     plot(P(:,1),P(:,2), '.', 'Color', trialcol(itt,:))
        fill(P_cleaned(k,1),P_cleaned(k,2), trialcol(itt,:), 'FaceAlpha', 0.3, 'EdgeColor', ecol);
    end
end

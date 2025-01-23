%% prepare files to be loaded in UMAP_KNN_decoder.ipynb : residual-rescaled spike counts for each session
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])
%%
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_V1.mat'], 'lmlvslope', 'lmlvyintercept')


for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    disperses = -1:0.1:2;
    % disperses = 0:1:2;

    testt = [106 107 110 111];
    inferencett = [1105 1109];
    Ntt = numel(testt);
    Nhireptt = numel(hireptt);

    try
        tempspkcnt = cat(3,spkcntIChiV1agg{ises}{:});
        Nrep = size(tempspkcnt,1);
        Nneu = size(tempspkcnt,2);
    catch
        % trial repetitions were not the same across trial types
        csz = cellfun(@size, spkcntIChiV1agg{ises}, 'UniformOutput', false);
        csz = cat(1,csz{:});
        Nrep = min(csz(:,1));
        if all(csz(:,2)==csz(1,2))
            Nneu = csz(1,2);
        else
            error('check number of neurons in session %d', ises)
        end
        tempspkcnt = NaN(Nrep, Nneu, Nhireptt);
        for n = 1:Nhireptt
            tempspkcnt(:,:,n) = spkcntIChiV1agg{ises}{n}(1:Nrep,:);
        end
    end
    spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu

    spkcntlmlv = cell(size(disperses));
    for islope = 0:numel(disperses)
        if islope==0
            spkcntlmlvs = spkcntIChiV1agg{ises};
            tempR = reshape(spkcntICtt, Nrep*Nhireptt, Nneu);
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);

            spkcntorig = tempR;
        else
            % fit log(mean) vs log(var)
            spkcntres = spkcntICtt - mean(spkcntICtt,1); % Nrep * Nnstt * Nneu
            spkcntmu = mean(spkcntICtt,1); % 1XNimg X Nneurons
            spkcntvar = var(spkcntICtt,0,1); % 1XNimg X Nneurons
            temp = spkcntvar; temp(spkcntvar==0)=NaN;
            totvar = nanmean(temp,2);

            tempx = log10(spkcntmu);
            tempx(spkcntmu==0) = NaN;
            meanx = squeeze(nanmean(tempx,3)); % average across neurons: 1XNimg

            Avec = lmlvslope(ises,:);
            Bvec = lmlvyintercept(ises,:);
            Cvec = disperses(islope)*ones(1,Nhireptt);
            Dvec = (Avec-Cvec).*meanx + Bvec;
            newspkcntvar = 10.^( (Cvec./Avec).*(log10(spkcntvar)-Bvec) + Dvec);
            newspkcntres = spkcntres .* sqrt(newspkcntvar./spkcntvar);
            newspkcntICtt = mean(spkcntICtt,1)+newspkcntres;

            tempR = reshape(newspkcntICtt, Nrep*Nhireptt, Nneu);
            tempR(~isfinite(tempR)) = NaN;
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);

            spkcntlmlv{islope} = tempR;
        end
    end

    save([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'], 'disperses', 'trialorder', 'spkcntorig', 'spkcntlmlv');

end

%% run UMAP_KNN_decoder.ipynb
%% analyze output of UMAP_KNN_decoder.ipynb 
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};

testacc_origagg = [];
infperf_origagg = [];
infscore_origagg = [];
testacc_lmlvsagg = [];
infperf_lmlvsagg = [];
infscore_lmlvsagg = [];

for ises = 1:numel(nwbsessions)
mousedate = nwbsessions{ises};
fprintf('%s %d\n', mousedate, ises)
pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

load([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'])
load([pathpp 'UMAP_KNN_decoding_V1RS_lmlvslopes.mat'])
Nsplits = size(UMAPorig.probe_predlabels,1);

testt = [106 107 110 111];
inferencett = [1105 1109];
Ntt = numel(testt);

infperf = zeros( numel(inferencett), numel(testt), Nsplits);
for isplit = 1:Nsplits
    for ii = 1:numel(inferencett)
        trialsoi = trialorder==inferencett(ii);
        [v,c] = uniquecnt(UMAPorig.probe_predlabels(isplit,trialsoi));
        infperf(ii, ismember(testt,v), isplit) = c/nnz(trialsoi);
    end
end
infscore = squeeze( (( infperf(1,1,:)-infperf(1,2,:) )+( infperf(2,4,:)-infperf(2,3,:) ))/2 );
fprintf('as-is data\n')
disp(mean(infperf,3))
disp(mean(infscore))

infperf_orig = infperf;
infscore_orig = infscore;
testacc_orig = reshape(UMAPorig.accuracy,Nsplits,1);

infperf_lmlvs = zeros( numel(inferencett), numel(testt), Nsplits, numel(disperses));
infscore_lmlvs = NaN(Nsplits, numel(disperses));
testacc_lmlvs = NaN(Nsplits, numel(disperses));
for islope = 1:numel(disperses)
    infperf = zeros( numel(inferencett), numel(testt), Nsplits);
    for isplit = 1:Nsplits
        for ii = 1:numel(inferencett)
            trialsoi = trialorder==inferencett(ii);
            [v,c] = uniquecnt(UMAPlmlv{islope}.probe_predlabels(isplit,trialsoi));
            infperf(ii, ismember(testt,v), isplit) = c/nnz(trialsoi);
        end
    end
    infscore = squeeze( (( infperf(1,1,:)-infperf(1,2,:) )+( infperf(2,4,:)-infperf(2,3,:) ))/2 );
    fprintf('Slope %.2f\n', disperses(islope))
    disp(mean(infperf,3))
    disp(mean(infscore))

    infperf_lmlvs(:,:,:,islope) = infperf;
    infscore_lmlvs(:,islope) = infscore;
testacc_lmlvs(:,islope) = UMAPlmlv{islope}.accuracy;
end

testacc_origagg = cat(2, testacc_origagg, testacc_orig);
testacc_lmlvsagg = cat(3, testacc_lmlvsagg, testacc_lmlvs);

infperf_origagg = cat(4, infperf_origagg, infperf_orig);
infscore_origagg = cat(2, infscore_origagg, infscore_orig);
infperf_lmlvsagg = cat(5, infperf_lmlvsagg, infperf_lmlvs);
infscore_lmlvsagg = cat(3, infscore_lmlvsagg, infscore_lmlvs);

end

%%
slope2plt = [0 1 2];
figure
for isplit = 1:Nsplits
subplot(4,Nsplits,isplit)
scatter( squeeze(UMAPorig.probe_embeddings(isplit,:,1)), squeeze(UMAPorig.probe_embeddings(isplit,:,2)), 10, trialorder, 'filled')
axis square
title('As-Is')
for s = 1:numel(slope2plt)
    islope = disperses==slope2plt(s);
subplot(4,Nsplits,Nsplits*s+isplit)
scatter( squeeze(UMAPlmlv{islope}.probe_embeddings(isplit,:,1)), squeeze(UMAPlmlv{islope}.probe_embeddings(isplit,:,2)), 10, trialorder, 'filled')
axis square
title(sprintf('Slope %.2f', disperses(islope)))
end
end

slope2plt = [0 1 2];
figure
for isplit = 1:Nsplits
subplot(4,Nsplits,isplit)
scatter( squeeze(UMAPorig.test_embeddings(isplit,:,1)), squeeze(UMAPorig.test_embeddings(isplit,:,2)), 10, UMAPorig.test_truelabels(isplit,:), 'filled')
axis square
title('As-Is')
for s = 1:numel(slope2plt)
    islope = disperses==slope2plt(s);
subplot(4,Nsplits,Nsplits*s+isplit)
scatter( squeeze(UMAPlmlv{islope}.test_embeddings(isplit,:,1)), squeeze(UMAPlmlv{islope}.test_embeddings(isplit,:,2)), 10, UMAPlmlv{islope}.test_truelabels(isplit,:), 'filled')
axis square
title(sprintf('Slope %.2f', disperses(islope)))
end
end
%%
% figure
% subplot(1,2,1)
% plot(disperses, squeeze(mean(infscore_lmlvsagg, 1)) )
% subplot(1,2,2)
% plot(disperses, squeeze(mean(infscore_lmlvsagg, 1))' )


figure; 
subplot(2,2,1)
hold all
plot(disperses, squeeze(mean(testacc_lmlvsagg, 1)))
plot(disperses, squeeze(mean(testacc_lmlvsagg, [1, 3])), 'k-', 'LineWidth', 2)
set(gca, 'XGrid', 'on', 'YGrid', 'on')
ylabel('test accuracy')
subplot(2,2,2)
hold all
plot(disperses, squeeze(mean(infscore_lmlvsagg, 1)))
plot(disperses, squeeze(mean(infscore_lmlvsagg, [1, 3])), 'k-', 'LineWidth', 2)
set(gca, 'XGrid', 'on', 'YGrid', 'on')
ylabel('inference score')

subplot(2,2,3)
hold all
plot(disperses, squeeze(mean(infperf_lmlvsagg(1,1,:,:,:), 3)))
plot(disperses, squeeze(mean(infperf_lmlvsagg(1,1,:,:,:), [3,5])), 'k-', 'LineWidth', 2)
set(gca, 'XGrid', 'on', 'YGrid', 'on')
ylabel('P(TRE1->IC1)')
subplot(2,2,4)
hold all
plot(disperses, squeeze(mean(infperf_lmlvsagg(2,4,:,:,:), 3)))
plot(disperses, squeeze(mean(infperf_lmlvsagg(2,4,:,:,:), [3,5])), 'k-', 'LineWidth', 2)
set(gca, 'XGrid', 'on', 'YGrid', 'on')
ylabel('P(TRE2->IC2)')

%%
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
    'UMAPorigagg', 'UMAPorigallagg', 'UMAPorig_unsupagg', 'UMAPorigall_unsupagg', ...
    'UMAPlmlvacc', 'UMAPlmlvallacc', 'UMAPlmlv_unsupacc', 'UMAPlmlvall_unsupacc')

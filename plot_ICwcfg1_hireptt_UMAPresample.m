% randomly sample one trial from each trial type: exrep Ntt X Ndim matrices
% proscrutes dsitance between a pair of exrep matrices
% based on a pair of exrep matrices,
% construct Ntt X Ntt umap euclidean distance matrix,
% calculate cosine similarity between euclidean matrices
% repeat Nboot times then average cosine similarities

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

%% between pairs of sessions
[v,c] = cellfun(@uniquecnt, trialorderacc, 'uniformoutput', false);
nreps = cat(1, c{:});
Nrep = min(nreps(:));

Nboot = 100;

umappairpros = struct();
umappairpros.unsup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairpros.unsup_origtest = NaN(Nsessions, Nsessions, Nboot);
umappairpros.semisup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairpros.semisup_origtest = NaN(Nsessions, Nsessions, Nboot);

umappairpros.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umappairpros.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umappairpros.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umappairpros.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairpros.unsup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairpros.unsup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairpros.semisup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairpros.semisup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
end

umappairdistcorr = struct();
umappairdistcorr.unsup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairdistcorr.unsup_origtest = NaN(Nsessions, Nsessions, Nboot);
umappairdistcorr.semisup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairdistcorr.semisup_origtest = NaN(Nsessions, Nsessions, Nboot);

umappairdistcorr.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umappairdistcorr.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umappairdistcorr.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umappairdistcorr.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairdistcorr.unsup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcorr.unsup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcorr.semisup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcorr.semisup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
end

umappairdistcossim = struct();
umappairdistcossim.unsup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairdistcossim.unsup_origtest = NaN(Nsessions, Nsessions, Nboot);
umappairdistcossim.semisup_origall = NaN(Nsessions, Nsessions, Nboot);
umappairdistcossim.semisup_origtest = NaN(Nsessions, Nsessions, Nboot);

umappairdistcossim.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umappairdistcossim.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umappairdistcossim.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umappairdistcossim.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairdistcossim.unsup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcossim.unsup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcossim.semisup_lmlvall{islope} = NaN(Nsessions, Nsessions, Nboot);
    umappairdistcossim.semisup_lmlvtest{islope} = NaN(Nsessions, Nsessions, Nboot);
end

for supind = 1:4
    tic
    switch supind
        case 1
            UMAPoption = 'semisup_origall';
            slopes2iter = 1;
        case 2
            UMAPoption = 'unsup_origall';
            slopes2iter = 1;
        case 3
            UMAPoption = 'semisup_lmlvall';
            slopes2iter = 1:numel(lmlvslope_list);
        case 4
            UMAPoption = 'unsup_lmlvall';
            slopes2iter = 1:numel(lmlvslope_list);
    end
    disp(UMAPoption)
    for islope = slopes2iter
    for ises = 1:Nsessions
        for jses = 1:Nsessions
            
            switch UMAPoption
                case 'semisup_origall'
                    umap1 = double( UMAPorigallagg(ises).embeddings );
                    umap2 = double( UMAPorigallagg(jses).embeddings );
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
                case 'unsup_origall'
                    umap1 = double( UMAPorigall_unsupagg(ises).embeddings );
                    umap2 = double( UMAPorigall_unsupagg(jses).embeddings );
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
                case 'semisup_lmlvall'
                    umap1 = double( UMAPlmlvallacc{ises, islope}.embeddings );
                    umap2 = double( UMAPlmlvallacc{jses, islope}.embeddings );
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
                case 'unsup_lmlvall'
                    umap1 = double( UMAPlmlvall_unsupacc{ises, islope}.embeddings );
                    umap2 = double( UMAPlmlvall_unsupacc{jses, islope}.embeddings );
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
            end
            umap1rs = NaN(Nrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord1==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind)));
                umap1rs(:,itt,:) = umap1(trialsoind,:);
            end
            umap2rs = NaN(Nrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord2==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind)));
                umap2rs(:,itt,:) = umap2(trialsoind,:);
            end
            prosvec = NaN(Nboot,1);
            distcorrvec = NaN(Nboot,1);
            distcossimvec = NaN(Nboot,1);
            for iboot = 1:Nboot
                ij=randperm(Nrep,2);
                umap1ex1 = squeeze(umap1rs(ij(1),:,:));
                umap1ex2 = squeeze(umap1rs(ij(2),:,:));
                ij=randperm(Nrep,2);
                umap2ex1 = squeeze(umap2rs(ij(1),:,:));
                umap2ex2 = squeeze(umap2rs(ij(2),:,:));
                
                prosvec(iboot) = procrustes(umap1ex1, umap2ex1);
                
                eucdist1 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist1 = eucdist1 + (umap1ex1(:,idim) - umap1ex2(:,idim)').^2;
                end
                eucdist1 = sqrt(eucdist1);
                
                eucdist2 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist2 = eucdist2 + (umap2ex1(:,idim) - umap2ex2(:,idim)').^2;
                end
                eucdist2 = sqrt(eucdist2);
                
                distcorrvec(iboot) = corr(eucdist1(:), eucdist2(:));
                distcossimvec(iboot) = dot(eucdist1(:), eucdist2(:))/(norm(eucdist1(:))*norm(eucdist2(:)));
            end
            if contains(UMAPoption, 'orig')
            umappairpros.(UMAPoption)(ises,jses,:) = prosvec;
            umappairdistcorr.(UMAPoption)(ises,jses,:) = distcorrvec;
            umappairdistcossim.(UMAPoption)(ises,jses,:) = distcossimvec;
            else
            umappairpros.(UMAPoption){islope}(ises,jses,:) = prosvec;
            umappairdistcorr.(UMAPoption){islope}(ises,jses,:) = distcorrvec;
            umappairdistcossim.(UMAPoption){islope}(ises,jses,:) = distcossimvec;
            end
        end
    end
    end
    toc
end

Ntestrep = floor(Nrep/Nsplits);
for supind = 1:4
    tic
    switch supind
        case 1
            UMAPoption = 'semisup_origtest';
            slopes2iter = 1;
        case 2
            UMAPoption = 'unsup_origtest';
            slopes2iter = 1;
        case 3
            UMAPoption = 'semisup_lmlvtest';
            slopes2iter = 1:numel(lmlvslope_list);
        case 4
            UMAPoption = 'unsup_lmlvtest';
            slopes2iter = 1:numel(lmlvslope_list);
    end
    disp(UMAPoption)
    for islope = slopes2iter
    for ises = 1:Nsessions
        for jses = 1:Nsessions
            umap1rsc = cell(Nsplits,1);
            umap2rsc = cell(Nsplits,1);
            for isplit = 1:Nsplits
            switch UMAPoption
                case 'semisup_origtest'
                    umap1 = double( squeeze(UMAPorigagg(ises).test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPorigagg(jses).test_embeddings(isplit,:,:)) );
                    
                    trialord1 = UMAPorigagg(ises).test_truelabels(isplit,:);
                    trialord2 = UMAPorigagg(jses).test_truelabels(isplit,:);
                case 'unsup_origtest'
                    umap1 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPorig_unsupagg(jses).test_embeddings(isplit,:,:)) );
                    
                    trialord1 = UMAPorig_unsupagg(ises).test_truelabels(isplit,:);
                    trialord2 = UMAPorig_unsupagg(jses).test_truelabels(isplit,:);
                case 'semisup_lmlvtest'
                    umap1 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPlmlvacc{jses,islope}.test_embeddings(isplit,:,:)) );
                    
                    trialord1 = UMAPlmlvacc{ises,islope}.test_truelabels(isplit,:);
                    trialord2 = UMAPlmlvacc{jses,islope}.test_truelabels(isplit,:);
                case 'unsup_lmlvtest'
                    umap1 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPlmlv_unsupacc{jses,islope}.test_embeddings(isplit,:,:)) );
                    
                    trialord1 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(isplit,:);
                    trialord2 = UMAPlmlv_unsupacc{jses,islope}.test_truelabels(isplit,:);
            end
            umap1rsc{isplit} = NaN(Ntestrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord1==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind), Ntestrep));
                umap1rsc{isplit}(:,itt,:) = umap1(trialsoind,:);
            end
            umap2rsc{isplit} = NaN(Ntestrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord2==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind), Ntestrep));
                umap2rsc{isplit}(:,itt,:) = umap2(trialsoind,:);
            end
            end
            
            prosvec = NaN(Nboot,1);
            distcorrvec = NaN(Nboot,1);
            distcossimvec = NaN(Nboot,1);
            for iboot = 1:Nboot
                isplit = mod(iboot, Nsplits)+1;
                
                ij=randperm(Ntestrep,2);
                umap1ex1 = squeeze(umap1rsc{isplit}(ij(1),:,:));
                umap1ex2 = squeeze(umap1rsc{isplit}(ij(2),:,:));
                ij=randperm(Ntestrep,2);
                umap2ex1 = squeeze(umap2rsc{isplit}(ij(1),:,:));
                umap2ex2 = squeeze(umap2rsc{isplit}(ij(2),:,:));
                
                prosvec(iboot) = procrustes(umap1ex1, umap2ex1);
                
                eucdist1 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist1 = eucdist1 + (umap1ex1(:,idim) - umap1ex2(:,idim)').^2;
                end
                eucdist1 = sqrt(eucdist1);
                
                eucdist2 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist2 = eucdist2 + (umap2ex1(:,idim) - umap2ex2(:,idim)').^2;
                end
                eucdist2 = sqrt(eucdist2);
                
                distcorrvec(iboot) = corr(eucdist1(:), eucdist2(:));
                distcossimvec(iboot) = dot(eucdist1(:), eucdist2(:))/(norm(eucdist1(:))*norm(eucdist2(:)));
            end
            if contains(UMAPoption, 'orig')
            umappairpros.(UMAPoption)(ises,jses,:) = prosvec;
            umappairdistcorr.(UMAPoption)(ises,jses,:) = distcorrvec;
            umappairdistcossim.(UMAPoption)(ises,jses,:) = distcossimvec;
            else
            umappairpros.(UMAPoption){islope}(ises,jses,:) = prosvec;
            umappairdistcorr.(UMAPoption){islope}(ises,jses,:) = distcorrvec;
            umappairdistcossim.(UMAPoption){islope}(ises,jses,:) = distcossimvec;
            end
        end
    end
    end
    toc
end

%% between cross-validation sets
umappairpros.unsup_origcv = NaN(Nsessions, Nboot);
umappairpros.semisup_origcv = NaN(Nsessions, Nboot);

umappairpros.unsup_lmlvcv = cell(1,numel(lmlvslope_list));
umappairpros.semisup_lmlvcv = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairpros.unsup_lmlvcv{islope} = NaN(Nsessions, Nboot);
    umappairpros.semisup_lmlvcv{islope} = NaN(Nsessions, Nboot);
end

umappairdistcorr.unsup_origcv = NaN(Nsessions, Nboot);
umappairdistcorr.semisup_origcv = NaN(Nsessions, Nboot);

umappairdistcorr.unsup_lmlvcv = cell(1,numel(lmlvslope_list));
umappairdistcorr.semisup_lmlvcv = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairdistcorr.unsup_lmlvcv{islope} = NaN(Nsessions, Nboot);
    umappairdistcorr.semisup_lmlvcv{islope} = NaN(Nsessions, Nboot);
end

umappairdistcossim.unsup_origcv = NaN(Nsessions, Nboot);
umappairdistcossim.semisup_origcv = NaN(Nsessions, Nboot);

umappairdistcossim.unsup_lmlvcv = cell(1,numel(lmlvslope_list));
umappairdistcossim.semisup_lmlvcv = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umappairdistcossim.unsup_lmlvcv{islope} = NaN(Nsessions, Nboot);
    umappairdistcossim.semisup_lmlvcv{islope} = NaN(Nsessions, Nboot);
end

Ntestrep = floor(Nrep/Nsplits);
for supind = 1:4
    tic
    switch supind
        case 1
            UMAPoption = 'semisup_origcv';
            slopes2iter = 1;
        case 2
            UMAPoption = 'unsup_origcv';
            slopes2iter = 1;
        case 3
            UMAPoption = 'semisup_lmlvcv';
            slopes2iter = 1:numel(lmlvslope_list);
        case 4
            UMAPoption = 'unsup_lmlvcv';
            slopes2iter = 1:numel(lmlvslope_list);
    end
    disp(UMAPoption)
    for islope = slopes2iter
        for ises = 1:Nsessions
            switch UMAPoption
                case 'semisup_origcv'
                    umap1 = double( squeeze(UMAPorigagg(ises).test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPorigagg(ises).test_embeddings(2,:,:)) );
                    
                    trialord1 = UMAPorigagg(ises).test_truelabels(1,:);
                    trialord2 = UMAPorigagg(ises).test_truelabels(2,:);
                case 'unsup_origcv'
                    umap1 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(2,:,:)) );
                    
                    trialord1 = UMAPorig_unsupagg(ises).test_truelabels(1,:);
                    trialord2 = UMAPorig_unsupagg(ises).test_truelabels(2,:);
                case 'semisup_lmlvcv'
                    umap1 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(2,:,:)) );
                    
                    trialord1 = UMAPlmlvacc{ises,islope}.test_truelabels(1,:);
                    trialord2 = UMAPlmlvacc{ises,islope}.test_truelabels(2,:);
                case 'unsup_lmlvcv'
                    umap1 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(2,:,:)) );
                    
                    trialord1 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(1,:);
                    trialord2 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(2,:);
            end
            umap1rs = NaN(Ntestrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord1==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind), Ntestrep));
                umap1rs(:,itt,:) = umap1(trialsoind,:);
            end
            umap2rs = NaN(Ntestrep, Nhireptt, Ndims);
            for itt = 1:Nhireptt
                trialsoind = find(trialord2==hireptt(itt));
                trialsoind = trialsoind(randperm(numel(trialsoind), Ntestrep));
                umap2rs(:,itt,:) = umap2(trialsoind,:);
            end
            
            prosvec = NaN(Nboot,1);
            distcorrvec = NaN(Nboot,1);
            distcossimvec = NaN(Nboot,1);
            for iboot = 1:Nboot
                ij=randperm(Ntestrep,2);
                umap1ex1 = squeeze(umap1rs(ij(1),:,:));
                umap1ex2 = squeeze(umap1rs(ij(2),:,:));
                ij=randperm(Ntestrep,2);
                umap2ex1 = squeeze(umap2rs(ij(1),:,:));
                umap2ex2 = squeeze(umap2rs(ij(2),:,:));
                
                prosvec(iboot) = procrustes(umap1ex1, umap2ex1);
                
                eucdist1 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist1 = eucdist1 + (umap1ex1(:,idim) - umap1ex2(:,idim)').^2;
                end
                eucdist1 = sqrt(eucdist1);
                
                eucdist2 = zeros(Nhireptt);
                for idim = 1:Ndims
                    eucdist2 = eucdist2 + (umap2ex1(:,idim) - umap2ex2(:,idim)').^2;
                end
                eucdist2 = sqrt(eucdist2);
                
                distcorrvec(iboot) = corr(eucdist1(:), eucdist2(:));
                distcossimvec(iboot) = dot(eucdist1(:), eucdist2(:))/(norm(eucdist1(:))*norm(eucdist2(:)));
            end
            if contains(UMAPoption, 'orig')
                umappairpros.(UMAPoption)(ises,:) = prosvec;
                umappairdistcorr.(UMAPoption)(ises,:) = distcorrvec;
                umappairdistcossim.(UMAPoption)(ises,:) = distcossimvec;
            else
                umappairpros.(UMAPoption){islope}(ises,:) = prosvec;
                umappairdistcorr.(UMAPoption){islope}(ises,:) = distcorrvec;
                umappairdistcossim.(UMAPoption){islope}(ises,:) = distcossimvec;
            end
        end
    end
    toc
end

%% between sessions, UMAP all trials
similmetrics = {'proscrutes', 'distcorr', 'distcossim'};
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', 'UMAP all trials, similarity between sessions', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    for imet = 1:numel(similmetrics)
        switch similmetrics{imet}
            case 'proscrutes'
                tempumappairsim = umappairpros;
                ylab = 'UMAP Proscrutes distance';
            case 'distcorr'
                tempumappairsim = umappairdistcorr;
                ylab = 'UMAP distance correlation';
            case 'distcossim'
                tempumappairsim = umappairdistcossim;
                ylab = 'UMAP distance cosine similarity';
        end
        
    tempmat = mean(tempumappairsim.([UMAPopt, '_origall']),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    repsimavgorig = nanmean(tempmat);
    repsimsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    repsimavglmlvs = NaN(1,Nslopes);
    repsimsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(tempumappairsim.([UMAPopt, '_lmlvall']){islope},3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        repsimavglmlvs(islope) = nanmean(tempmat);
        repsimsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    
    subplot(2,numel(similmetrics),(supind-1)*numel(similmetrics)+imet)
    hold all
    errorbar(lmlvslope_list, repsimavglmlvs, repsimsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, repsimavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel(ylab)
    title([UMAPopt ' ' similmetrics{imet}])
    end
end

%% between sessions, cross-validated UMAP test trials
similmetrics = {'proscrutes', 'distcorr', 'distcossim'};
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', 'cross-validated UMAP test trials, similarity between sessions', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    for imet = 1:numel(similmetrics)
        switch similmetrics{imet}
            case 'proscrutes'
                tempumappairsim = umappairpros;
                ylab = 'UMAP Proscrutes distance';
            case 'distcorr'
                tempumappairsim = umappairdistcorr;
                ylab = 'UMAP distance correlation';
            case 'distcossim'
                tempumappairsim = umappairdistcossim;
                ylab = 'UMAP distance cosine similarity';
        end
        
    tempmat = mean(tempumappairsim.([UMAPopt, '_origtest']),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    repsimavgorig = nanmean(tempmat);
    repsimsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    repsimavglmlvs = NaN(1,Nslopes);
    repsimsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(tempumappairsim.([UMAPopt, '_lmlvtest']){islope},3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        repsimavglmlvs(islope) = nanmean(tempmat);
        repsimsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    
    subplot(2,numel(similmetrics),(supind-1)*numel(similmetrics)+imet)
    hold all
    errorbar(lmlvslope_list, repsimavglmlvs, repsimsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, repsimavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel(ylab)
    title([UMAPopt ' ' similmetrics{imet}])
    end
end

%% between test sets
similmetrics = {'proscrutes', 'distcorr', 'distcossim'};
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', 'cross-validated UMAP similarity between test sets', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    for imet = 1:numel(similmetrics)
        switch similmetrics{imet}
            case 'proscrutes'
                tempumappairsim = umappairpros;
                ylab = 'UMAP Proscrutes distance';
            case 'distcorr'
                tempumappairsim = umappairdistcorr;
                ylab = 'UMAP distance correlation';
            case 'distcossim'
                tempumappairsim = umappairdistcossim;
                ylab = 'UMAP distance cosine similarity';
        end
        
    tempmat = mean(tempumappairsim.([UMAPopt, '_origcv']),2);
    repsimavgorig = nanmean(tempmat);
    repsimsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    repsimavglmlvs = NaN(1,Nslopes);
    repsimsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(tempumappairsim.([UMAPopt, '_lmlvcv']){islope},2);
        repsimavglmlvs(islope) = nanmean(tempmat);
        repsimsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    
    subplot(2,numel(similmetrics),(supind-1)*numel(similmetrics)+imet)
    hold all
    errorbar(lmlvslope_list, repsimavglmlvs, repsimsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, repsimavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel(ylab)
    title([UMAPopt ' ' similmetrics{imet}])
    end
end

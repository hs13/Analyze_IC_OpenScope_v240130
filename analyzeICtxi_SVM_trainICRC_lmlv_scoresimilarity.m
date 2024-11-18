if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_V1.mat'])

% if 0, keep all neurons; if 1, exclude zero variance neurons in train trial
% types; if 2 exclude zero variance neurons in all trial types
% excludeneuvar0 = 2;
fprintf('neuron exclusion criterion %d\n', excludeneuvar0)

for ises = numel(nwbsessions):-1:1
    clearvars -except excludeneuvar0 ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    sesclk = tic;
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    fprintf('%s %d\n', mousedate, ises)
    
    pltses = true;
    preproc = 'meancenter';
    whichSVMkernel = 'Linear';
    svmdesc = 'trainICRC';
    switch excludeneuvar0
        case 0
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl_lmlvslopes.mat');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes_incl.mat');
        case 1
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes.mat_excltt');
        case 2
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes_excl.mat');
        otherwise
            error('excludeneuvar0 option not recognized')
    end
    
    load(svmfn)
    load(svmmdlfn)
    load(svmlmlvfn)
    load(svmmdllmlvfn)
    
    % parametrically change proportion silenced
    propneusilvec = 0:0.1:1;
    % randomly select neurons to silence: sample 100X
    Nsamples = 10; % estimated ~30min for 100 samples
    % propneusilvec = 0;
    % Nsamples = 1;
    Nneurons = SVMtrainICRC.Nneurons;
    Nsplits = size(SVMtrainICRC.spkcnt.testtrialinds,2);
    
    testt = [106,107,110,111];
    inferencett = [1105 1109];
    Nhireptt = numel(hireptt);
    Ntt = numel(testt);
    Nttrain = size(SVMtrainICRC.spkcnt.traintrialinds,1)/Ntt;
    Nttest = size(SVMtrainICRC.spkcnt.testtrialinds,1)/Ntt;
    
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
    
    % initialize
    rhoscorefields = {'train', 'test', 'simil', 'trainpair', 'testpair', 'similpair'};
    % rhoscorestats = {'avg', 'stdpool', 'semmean', 'medianpool', 'q1pool', 'q3pool'};
    rhoscorestats = {'avg', 'medianpool', 'prct'};
    % pool/mean indicates whether to pool or mean across k-fold splits
    meanvecscorerhoasis = struct(); % correlation between k-fold splits
    meanvecscorerholmlvs = struct();
    for r = 1:numel(rhoscorestats)
        meanvecscorerhoasis.(rhoscorestats{r}) = zeros(Nhireptt,1);
        meanvecscorerholmlvs.(rhoscorestats{r}) = zeros(Nhireptt,numel(disperses));
    end
    meanvecrankpointasis = zeros(Nhireptt,Ntt);
    meanvecrankpointlmlvs = zeros(Nhireptt,Ntt,numel(disperses));
    rhoscoreasis = struct();
    rhoscorelmlvs = struct();
    rhoxneusubasis = struct();
    rhoxneusublmlvs = struct();
    for f = 1:numel(rhoscorefields)
        if contains(rhoscorefields{f}, 'train')
            tempNt = Ntt;
            tempNtrials = Nttrain;
        elseif contains(rhoscorefields{f}, 'test')
            tempNt = Ntt;
            tempNtrials = Nttest;
        elseif contains(rhoscorefields{f}, 'simil')
            tempNt = Nhireptt;
            tempNtrials = Nrep;
        else
            error([rhoscorefields{f} ' not recognized'])
        end
        for r = 1:numel(rhoscorestats)
            rhoscoreasis.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(Nsamples, numel(propneusilvec), tempNt);
            rhoscorelmlvs.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(Nsamples, numel(propneusilvec), tempNt, numel(disperses));
            
            rhoxneusubasis.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(tempNtrials, numel(propneusilvec), tempNt);
            rhoxneusublmlvs.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(tempNtrials, numel(propneusilvec), tempNt, numel(disperses));
        end
    end
    
    similclk = tic;
    for islope = 0:numel(disperses)
        tic
        if islope==0
            tempR = reshape(spkcntICtt, Nrep*Nhireptt, Nneu)';
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);
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
            
            switch excludeneuvar0
                case 0
                    valneu = true(Nneu,1);
                case 1
                    valneu = squeeze(all(newspkcntvar(1,ismember(hireptt,testt),:)>0 & isfinite(newspkcntvar(1,ismember(hireptt,testt),:)), 2));
                case 2
                    valneu = squeeze(all(newspkcntvar(1,:,:)>0 & isfinite(newspkcntvar(1,:,:)), 2));
                otherwise
                    error('excludeneuvar0 option not recognized')
            end
            
            tempR = reshape(newspkcntICtt(:,:,valneu), Nrep*Nhireptt, nnz(valneu))';
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);
        end
        
        if islope==0
            SVMout = SVMtrainICRC;
            SVM_models = SVMtrainICRC_models.spkcnt;
            valneu = true(Nneurons,1);
        else
            SVMout = SVMtrainICRC_lmlvs(islope);
            SVM_models = SVMtrainICRC_models_lmlvs(islope).spkcnt;
            valneu = SVMout.valneu;
        end
        
        meanveclabelcv = NaN(Nhireptt, Nsplits);
        meanvecscorecv = NaN(Nhireptt, Nsplits, Ntt);
        meanvecrankptcv = zeros(Nhireptt, Nsplits, Ntt);% higher scores have higher value
        for isplit = 1:Nsplits
            traintrialinds = SVMout.spkcnt.traintrialinds(:,isplit);
            testtrialinds = SVMout.spkcnt.testtrialinds(:,isplit);
            switch preproc
                case 'none'
                    Tp = tempR';
                case 'zscore'
                    % Z-score
                    trainRmean = mean(tempR(:,traintrialinds),2);
                    trainRstd = std(tempR(:,traintrialinds),0,2);
                    
                    Tp = ( (tempR-trainRmean)./trainRstd )';
                case 'minmax'
                    trainRmin = min(tempR(:,traintrialinds),[],2);
                    trainRrange = range(tempR(:,traintrialinds),2);
                    
                    Tp = ( (tempR-trainRmin)./trainRrange )';
                case 'meancenter'
                    trainRmean = mean(tempR(:,traintrialinds),2);
                    Tp = (tempR-trainRmean)';
            end
            Tp(isnan(Tp))=0; % Ntrials * Nneurons
            
            for itt = 1:Nhireptt
                meanvec = mean(Tp(trialorder==hireptt(itt),:),1);
                [meanveclabel,meanvecscore] = predict(SVM_models{isplit}, meanvec);
                meanveclabelcv(itt, isplit) = meanveclabel;
                meanvecscorecv(itt, isplit, :) = meanvecscore;
                [sv,si]=sort(meanvecscore, 'ascend'); % ascending rank: higher scores have higher value
                meanvecrankptcv(itt,isplit,si) = 1:Ntt;
            end
        end
        % [p,tbl,stats] = friedman(squeeze(meanvecrankcv(itt,:,:)));
        % [c,m,h] = multcompare(stats);
        % [p,tbl,stats] = friedman(squeeze(-meanvecscorecv(itt,:,:)), 'display', 'off');
        % [c,m,h] = multcompare(stats);
        % m is the same as squeeze(mean(meanvecrankcv,2))
        meanvecrankpointcv = squeeze(mean(meanvecrankptcv,2)); % higher scores have higher value
        
        meanvecrho = zeros(Nhireptt, nchoosek(Nsplits,2));
        for itt = 1:Nhireptt
            rho =  corr(squeeze(meanvecscorecv(itt,:,:))', 'type', 'spearman');
            rhorank =  corr(squeeze(meanvecrankptcv(itt,:,:))', 'type', 'spearman');
            if ~isequal(rho, rhorank)
                error('sanity check failed')
            end
            meanvecrho(itt,:) = rho(triu(true(Nsplits),1));
        end
        
        if islope==0
            meanvecrankpointasis = meanvecrankpointcv;
            meanvecscorerhoasis.avg = mean(meanvecrho,2);
            meanvecscorerhoasis.medianpool = median(meanvecrho,2);
            meanvecscorerhoasis.prct = mean(meanvecrho==1,2);
        else
            meanvecrankpointlmlvs(:,:,islope) = meanvecrankpointcv;
            meanvecscorerholmlvs.avg(:,islope) = mean(meanvecrho,2);
            meanvecscorerholmlvs.medianpool(:,islope) = median(meanvecrho,2);
            meanvecscorerholmlvs.prct(:,islope) = mean(meanvecrho==1,2);
        end
        
        for n = 1:length(propneusilvec)
            propneu2sil = propneusilvec(n);
            if propneu2sil==0
                neu2silmat = false(Nneurons,1);
            elseif propneu2sil==1
                neu2silmat = true(Nneurons,1);
            else
                neu2silmat = false(Nneurons, Nsamples);
                for s = 1:Nsamples
                    neurand = randperm(Nneurons, round(propneu2sil*Nneurons));
                    neu2silmat(neurand,s) = true;
                end
            end
            
            
            trainscorepd = zeros(Nsamples, Nttrain, Ntt, Nsplits, Ntt);
            testscorepd = zeros(Nsamples, Nttest, Ntt, Nsplits, Ntt);
            similscorepd = zeros(Nsamples, Nrep, Ntt, Nsplits, Nhireptt);
            
            for s = 1:size(neu2silmat,2) % Nsamples if 0<propneu2sil<1
                trainscore = zeros(Nttrain, Ntt, Nsplits, Ntt);
                testscore = zeros(Nttest, Ntt, Nsplits, Ntt);
                similscore = zeros(Nrep, Ntt, Nsplits, Nhireptt);
                for isplit = 1:Nsplits
                    traintrialinds = SVMout.spkcnt.traintrialinds(:,isplit);
                    testtrialinds = SVMout.spkcnt.testtrialinds(:,isplit);
                    switch preproc
                        case 'none'
                            Tp = tempR';
                        case 'zscore'
                            % Z-score
                            trainRmean = mean(tempR(:,traintrialinds),2);
                            trainRstd = std(tempR(:,traintrialinds),0,2);
                            
                            Tp = ( (tempR-trainRmean)./trainRstd )';
                        case 'minmax'
                            trainRmin = min(tempR(:,traintrialinds),[],2);
                            trainRrange = range(tempR(:,traintrialinds),2);
                            
                            Tp = ( (tempR-trainRmin)./trainRrange )';
                        case 'meancenter'
                            trainRmean = mean(tempR(:,traintrialinds),2);
                            Tp = (tempR-trainRmean)';
                    end
                    Tp(isnan(Tp))=0; % Ntrials * Nneurons
                    
                    trainlabs = SVMout.trialorder(traintrialinds);
                    testlabs = SVMout.trialorder(testtrialinds);
                    inftrials = ismember(SVMout.trialorder, inferencett);
                    
                    Xsilrand = Tp;
                    Xsilrand(:,neu2silmat(valneu,s)) = 0;
                    [templabel,tempscore] = predict(SVM_models{isplit}, Xsilrand);
                    
                    % sanity check
                    if propneu2sil==0
                        if ~isequaln(templabel, SVMout.spkcnt.all.label(:,isplit))
                            error('check that Tp was calculated correctly')
                        end
                    end
                    
                    %                 if numel(unique(templabel))<Ntt
                    %                     warning('only %d/%d trial types returned by SVM, skipping...', numel(unique(templabel)), Ntt)
                    %                     continue
                    %                 end
                    
                    % train & test
                    for itt = 1:Ntt
                        trainttscore = tempscore(traintrialinds(trainlabs==testt(itt)),:);
                        trainscore(:,:,isplit,itt) = trainttscore;
                        
                        testttscore = tempscore(testtrialinds(testlabs==testt(itt)),:);
                        testscore(:,:,isplit,itt) = testttscore;
                    end
                    
                    % similarity inference
                    for itt = 1:Nhireptt
                        similttscore = tempscore(trialorder==hireptt(itt),:);
                        similscore(:,:,isplit,itt) = similttscore;
                    end
                end
                
                trainrhoscore = zeros(Nttrain, Nsplits, Ntt);
                testrhoscore = zeros(Nttest, Nsplits, Ntt);
                similrhoscore = zeros(Nrep, Nsplits, Nhireptt);
                trainrhoscorepair = zeros(Nttrain, nchoosek(Nsplits,2), Ntt);
                testrhoscorepair = zeros(Nttest, nchoosek(Nsplits,2), Ntt);
                similrhoscorepair = zeros(Nrep, nchoosek(Nsplits,2), Nhireptt);
                for itt = 1:Ntt
                    for itrial = 1:Nttrain
                        trainrhoscore(itrial,:,itt) = corr(squeeze(trainscore(itrial,:,:,itt)), meanvecrankpointcv(hireptt==testt(itt),:)', 'type', 'spearman');
                        rhomat = corr(squeeze(trainscore(itrial,:,:,itt)), 'type', 'spearman');
                        trainrhoscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                    for itrial = 1:Nttest
                        testrhoscore(itrial,:,itt) = corr(squeeze(testscore(itrial,:,:,itt)), meanvecrankpointcv(hireptt==testt(itt),:)', 'type', 'spearman');
                        rhomat = corr(squeeze(testscore(itrial,:,:,itt)), 'type', 'spearman');
                        testrhoscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                end
                for itt = 1:Nhireptt
                    for itrial = 1:Nrep
                        similrhoscore(itrial,:,itt) = corr(squeeze(similscore(itrial,:,:,itt)), meanvecrankpointcv(itt,:)', 'type', 'spearman');
                        rhomat = corr(squeeze(similscore(itrial,:,:,itt)), 'type', 'spearman');
                        similrhoscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                end
                
                % rhoscorefields = {'train', 'test', 'simil', 'tranpair', 'testpair', 'similpair'};
                % rhoscorestats = {'avg', 'stdpool', 'semmean', 'medianpool', 'q1pool', 'q3pool'};
                for f = 1:numel(rhoscorefields)
                    switch rhoscorefields{f}
                        case 'train'
                            temprhoscore = trainrhoscore;
                        case 'test'
                            temprhoscore = testrhoscore;
                        case 'simil'
                            temprhoscore = similrhoscore;
                        case 'trainpair'
                            temprhoscore = trainrhoscorepair;
                        case 'testpair'
                            temprhoscore = testrhoscorepair;
                        case 'similpair'
                            temprhoscore = similrhoscorepair;
                        otherwise
                            error([rhoscorefields{f} ' not recognized'])
                    end
                    for r = 1:numel(rhoscorestats)
                        switch rhoscorestats{r}
                            case 'avg'
                                temprhoscorestat = mean(temprhoscore,[1,2]);
                            case 'stdpool'
                                temprhoscorestat = std(reshape(temprhoscore,[],size(temprhoscore,3)),0,1);
                            case 'semmean'
                                temprhoscorestat = mean(std(temprhoscore,0,1)/sqrt(size(temprhoscore,1)), 2);
                            case 'medianpool'
                                temprhoscorestat = median(reshape(temprhoscore,[],size(temprhoscore,3)),1);
                            case 'q1pool'
                                temprhoscorestat = prctile(reshape(temprhoscore,[],size(temprhoscore,3)),25,1);
                            case 'q3pool'
                                temprhoscorestat = prctile(reshape(temprhoscore,[],size(temprhoscore,3)),75,1);
                            case 'prct'
                                temprhoscorestat = mean(temprhoscore==1,[1,2]);
                            otherwise
                                error([rhoscorestats{r} ' not recognized'])
                        end
                        if islope==0
                            rhoscoreasis.(rhoscorefields{f}).(rhoscorestats{r})(s,n,:) = squeeze(temprhoscorestat);
                        else
                            rhoscorelmlvs.(rhoscorefields{f}).(rhoscorestats{r})(s,n,:,islope) = squeeze(temprhoscorestat);
                        end
                    end
                end
                
                trainscorepd(s,:,:,:,:) = trainscore;
                testscorepd(s,:,:,:,:) = testscore;
                similscorepd(s,:,:,:,:) = similscore;
            end

            trainrhoxneusub = zeros(Nsamples, Nttrain, Nsplits, Ntt);
            testrhoxneusub = zeros(Nsamples, Nttest, Nsplits, Ntt);
            similrhoxneusub = zeros(Nsamples, Nrep, Nsplits, Nhireptt);
            trainrhoxneusubpair = zeros(nchoosek(Nsamples,2), Nttrain, Nsplits, Ntt);
            testrhoxneusubpair = zeros(nchoosek(Nsamples,2), Nttest, Nsplits, Ntt);
            similrhoxneusubpair = zeros(nchoosek(Nsamples,2), Nrep, Nsplits, Nhireptt);
            for itt = 1:Ntt
                for isplit = 1:Nsplits
                    for itrial = 1:Nttrain
                        trainrhoxneusub(:,itrial,isplit,itt) = corr(squeeze(trainscorepd(:,itrial,:,isplit,itt))', ...
                            meanvecrankpointcv(hireptt==testt(itt),:)', 'type', 'spearman');
                        rhomat = corr(squeeze(trainscorepd(:,itrial,:,isplit,itt))', 'type', 'spearman');
                        trainrhoxneusubpair(:,itrial,isplit,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                    
                    for itrial = 1:Nttest
                        testrhoxneusub(:,itrial,isplit,itt) = corr(squeeze(testscorepd(:,itrial,:,isplit,itt))', ...
                            meanvecrankpointcv(hireptt==testt(itt),:)', 'type', 'spearman');
                        rhomat = corr(squeeze(testscorepd(:,itrial,:,isplit,itt))', 'type', 'spearman');
                        testrhoxneusubpair(:,itrial,isplit,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                end
            end
            for itt = 1:Nhireptt
                for isplit = 1:Nsplits
                    for itrial = 1:Nrep
                        similrhoxneusub(:,itrial,isplit,itt) = corr(squeeze(similscorepd(:,itrial,:,isplit,itt))', ...
                            meanvecrankpointcv(itt,:)', 'type', 'spearman');
                        rhomat = corr(squeeze(similscorepd(:,itrial,:,isplit,itt))', 'type', 'spearman');
                        similrhoxneusubpair(:,itrial,isplit,itt) = rhomat(triu(true(size(rhomat)),1));
                    end
                end
            end
            
            % rhoscorefields = {'train', 'test', 'simil', 'tranpair', 'testpair', 'similpair'};
            % rhoscorestats = {'avg', 'stdpool', 'semmean', 'medianpool', 'q1pool', 'q3pool'};
            for f = 1:numel(rhoscorefields)
                switch rhoscorefields{f}
                    case 'train'
                        temprhoscore = trainrhoxneusub; % Nsamples*Nttrain*Nsplits*Ntt
                    case 'test'
                        temprhoscore = testrhoxneusub;
                    case 'simil'
                        temprhoscore = similrhoxneusub;
                    case 'trainpair'
                        temprhoscore = trainrhoxneusubpair;
                    case 'testpair'
                        temprhoscore = testrhoxneusubpair;
                    case 'similpair'
                        temprhoscore = similrhoxneusubpair;
                    otherwise
                        error([rhoscorefields{f} ' not recognized'])
                end
                temprhoscore = permute(temprhoscore,[1,3,2,4]);  % Nsamples*Nsplits*Nttrain*Ntt
                for r = 1:numel(rhoscorestats)
                    switch rhoscorestats{r}
                        case 'avg'
                            temprhoscorestat = mean(temprhoscore,[1,2]);
                        case 'stdpool'
                            temprhoscorestat = std(reshape(temprhoscore,[],size(temprhoscore,3),size(temprhoscore,4)),0,1);
                        case 'semmean'
                            temprhoscorestat = mean(std(temprhoscore,0,1)/sqrt(size(temprhoscore,1)), 2);
                        case 'medianpool'
                            temprhoscorestat = median(reshape(temprhoscore,[],size(temprhoscore,3),size(temprhoscore,4)),1);
                        case 'q1pool'
                            temprhoscorestat = prctile(reshape(temprhoscore,[],size(temprhoscore,3),size(temprhoscore,4)),25,1);
                        case 'q3pool'
                            temprhoscorestat = prctile(reshape(temprhoscore,[],size(temprhoscore,3),size(temprhoscore,4)),75,1);
                        case 'prct'
                            temprhoscorestat = mean(temprhoscore==1,[1,2]);
                        otherwise
                            error([rhoscorestats{r} ' not recognized'])
                    end
                    if islope==0
                        rhoxneusubasis.(rhoscorefields{f}).(rhoscorestats{r})(:,n,:) = squeeze(temprhoscorestat);
                    else
                        rhoxneusublmlvs.(rhoscorefields{f}).(rhoscorestats{r})(:,n,:,islope) = squeeze(temprhoscorestat);
                    end
                end
            end
            
            if islope==0
                fprintf('silence %.2f, as-is done\n', propneu2sil)
            else
                fprintf('silence %.2f, lmlv slope %.2f done\n', propneu2sil, disperses(islope))
            end
            toc
        end
    end
    
    save(similfn, 'excludeneuvar0', 'disperses', 'propneusilvec', 'meanvecrankpointasis', 'meanvecrankpointlmlvs', ...
        'meanvecscorerhoasis', 'meanvecscorerholmlvs', 'rhoscoreasis', 'rhoscorelmlvs', 'rhoxneusubasis', 'rhoxneusublmlvs')
    toc(similclk)
    
    if pltses
        figure
        subplot(2,3,1)
        hold all
        pl = plot(disperses, squeeze(meanvecscorerholmlvs.prct) );
        for itt = 1:numel(pl)
            h = squeeze(meanvecscorerhoasis.prct(itt));
            plot([disperses(1) disperses(end)], h*[1 1], 'Color', pl(itt).Color)
        end
        xlabel('LMLV slopes')
        ylabel('% rho(SVM score)==1')
        title('trial mean vector score consistency')
        for isp = [2,3,5,6]
            switch isp
                case 2
                    whichrhoscorefield = 'test';
                    ttoind = find(true(size(testt)));
                case 3
                    whichrhoscorefield = 'simil';
                    ttoind = find(~ismember(hireptt,testt));
                case 5
                    whichrhoscorefield = 'testpair';
                    ttoind = find(true(size(testt)));
                case 6
                    whichrhoscorefield = 'similpair';
                    ttoind = find(~ismember(hireptt,testt));
            end
            subplot(2,3,isp)
            hold all
            pl = plot(disperses, squeeze(rhoscorelmlvs.(whichrhoscorefield).prct(1,propneusilvec==0,ttoind,:) ) );
            yl = ylim;
            for ii = 1:numel(pl)
                if ismember(hireptt(ttoind(ii)), [1105 1109])
                    lw = 2;
                    disp(ii)
                else
                    lw = 0.2;
                end
                plot(disperses, squeeze(rhoscorelmlvs.(whichrhoscorefield).prct(1,propneusilvec==0,ttoind(ii),:) ), 'Color', pl(ii).Color, 'LineWidth', lw)
                h = squeeze(rhoscoreasis.(whichrhoscorefield).prct(1,propneusilvec==0,ttoind(ii) ));
                plot([disperses(1) disperses(end)], h*[1 1], 'Color', pl(ii).Color, 'LineWidth', lw)
                text(disperses(1), yl(1)+(numel(pl)-ii-1)*0.08*range(yl), num2str(hireptt(ttoind(ii))), 'Color', pl(ii).Color, 'FontSize', 10, 'VerticalAlignment', 'bottom')
            end
            ylim(yl)
            xlabel('LMLV slopes')
            ylabel('% rho(SVM score)==1')
            title(whichrhoscorefield)
        end
    end
    
end
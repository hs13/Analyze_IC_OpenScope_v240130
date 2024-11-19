%% score consistency across different splits of neurons (non-overlapping subsets)
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
excludeneuvar0 = 0;
fprintf('neuron exclusion criterion %d\n', excludeneuvar0)

for ises = numel(nwbsessions):-1:1
    clearvars -except excludeneuvar0 ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    sesclk = tic;
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    fprintf('%s %d\n', mousedate, ises)
    
    pltses = false;
    preproc = 'meancenter';
    whichSVMkernel = 'Linear';
    svmdesc = 'trainICRC';
    switch excludeneuvar0
        case 0
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_incl_lmlvslopes.mat');
            consistfn = strcat(pathpp, 'scoreconsistency_SVM_', svmdesc, '_lmlvslopes_incl.mat');
        case 1
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            consistfn = strcat(pathpp, 'scoreconsistency_SVM_', svmdesc, '_lmlvslopes.mat_excltt');
        case 2
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            consistfn = strcat(pathpp, 'scoreconsistency_SVM_', svmdesc, '_lmlvslopes_excl.mat');
        otherwise
            error('excludeneuvar0 option not recognized')
    end
    
    load(svmfn)
    load(svmmdlfn)
    load(svmlmlvfn)
    load(svmmdllmlvfn)
    
    % parametrically change proportion silenced
    propneu2sil = 0.5;
    % randomly select neurons to silence: sample 100X
    Nneudivs = 10; % estimated ~30min for 100 samples
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
    
    rhoxneudivasis = struct();
    rhoxneudivlmlvs = struct();
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
            rhoxneudivasis.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(Nneudivs, tempNt);
            rhoxneudivlmlvs.(rhoscorefields{f}).(rhoscorestats{r}) = NaN(Nneudivs, tempNt, numel(disperses));
        end
    end
    
    for islope = 0:numel(disperses)
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
        
        neu2silmat = false(Nneurons, Nneudivs);
        for s = 1:Nneudivs
            neurand = randperm(Nneurons, round(propneu2sil*Nneurons));
            neu2silmat(neurand,s) = true;
        end
        
        for s = 1:size(neu2silmat,2) % Nsamples if 0<propneu2sil<1
            trainscore = zeros(Nttrain, Ntt, Nsplits, Ntt);
            testscore = zeros(Nttest, Ntt, Nsplits, Ntt);
            similscore = zeros(Nrep, Ntt, Nsplits, Nhireptt);
            % complementary set of neurons
            trainscorecompl = zeros(Nttrain, Ntt, Nsplits, Ntt);
            testscorecompl = zeros(Nttest, Ntt, Nsplits, Ntt);
            similscorecompl = zeros(Nrep, Ntt, Nsplits, Nhireptt);
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

                Xsilcompl = Tp;
                Xsilcompl(:,~neu2silmat(valneu,s)) = 0;
                [compllabel,complscore] = predict(SVM_models{isplit}, Xsilcompl);                
                
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
                    
                    trainttscorecompl = complscore(traintrialinds(trainlabs==testt(itt)),:);
                    trainscorecompl(:,:,isplit,itt) = trainttscorecompl;
                    
                    testttscorecompl = complscore(testtrialinds(testlabs==testt(itt)),:);
                    testscorecompl(:,:,isplit,itt) = testttscorecompl;
                end
                
                % similarity inference
                for itt = 1:Nhireptt
                    similttscore = tempscore(trialorder==hireptt(itt),:);
                    similscore(:,:,isplit,itt) = similttscore;
                    
                    similttscorecompl = complscore(trialorder==hireptt(itt),:);
                    similscorecompl(:,:,isplit,itt) = similttscorecompl;
                end
            end
            
            trainconscorediag = zeros(Nttrain, Nsplits, Ntt);
            testconscorediag = zeros(Nttest, Nsplits, Ntt);
            similconscorediag = zeros(Nrep, Nsplits, Nhireptt);
            trainconscorepair = zeros(Nttrain, nchoosek(Nsplits,2), Ntt);
            testconscorepair = zeros(Nttest, nchoosek(Nsplits,2), Ntt);
            similconscorepair = zeros(Nrep, nchoosek(Nsplits,2), Nhireptt);
            for itt = 1:Ntt
                for itrial = 1:Nttrain
                    rhomat = corr(squeeze(trainscore(itrial,:,:,itt)), squeeze(trainscorecompl(itrial,:,:,itt)), 'type', 'spearman');
                    trainconscorediag(itrial,:,itt) = diag(rhomat);
                    trainconscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                end
                for itrial = 1:Nttest
                    rhomat = corr(squeeze(testscore(itrial,:,:,itt)), squeeze(testscorecompl(itrial,:,:,itt)), 'type', 'spearman');
                    testconscorediag(itrial,:,itt) = diag(rhomat);
                    testconscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                end
            end
            for itt = 1:Nhireptt
                for itrial = 1:Nrep
                    rhomat = corr(squeeze(similscore(itrial,:,:,itt)), squeeze(similscorecompl(itrial,:,:,itt)), 'type', 'spearman');
                    similconscorediag(itrial,:,itt) = diag(rhomat);
                    similconscorepair(itrial,:,itt) = rhomat(triu(true(size(rhomat)),1));
                end
            end
            
            % rhoscorefields = {'train', 'test', 'simil', 'tranpair', 'testpair', 'similpair'};
            % rhoscorestats = {'avg', 'stdpool', 'semmean', 'medianpool', 'q1pool', 'q3pool'};
            for f = 1:numel(rhoscorefields)
                switch rhoscorefields{f}
                    case 'train'
                        temprhoscore = trainconscorediag;
                    case 'test'
                        temprhoscore = testconscorediag;
                    case 'simil'
                        temprhoscore = similconscorediag;
                    case 'trainpair'
                        temprhoscore = trainconscorepair;
                    case 'testpair'
                        temprhoscore = testconscorepair;
                    case 'similpair'
                        temprhoscore = similconscorepair;
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
                        rhoxneudivasis.(rhoscorefields{f}).(rhoscorestats{r})(s,:) = squeeze(temprhoscorestat);
                    else
                        rhoxneudivlmlvs.(rhoscorefields{f}).(rhoscorestats{r})(s,:,islope) = squeeze(temprhoscorestat);
                    end
                end
            end
            
        end
        
        if islope==0
            fprintf('silence %.2f, as-is done\n', propneu2sil)
        else
            fprintf('silence %.2f, lmlv slope %.2f done\n', propneu2sil, disperses(islope))
        end
    end
    
    save(consistfn, 'excludeneuvar0', 'disperses', 'propneu2sil', 'rhoxneudivasis', 'rhoxneudivlmlvs')
    toc(sesclk)
    
    if pltses
        testt = [106,107,110,111];
        hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];
        figure
        for isp = 1:4
            switch isp
                case 1
                    whichrhoscorefield = 'test';
                    ttoind = find(true(size(testt)));
                case 2
                    whichrhoscorefield = 'simil';
                    ttoind = find(~ismember(hireptt,testt));
                case 3
                    whichrhoscorefield = 'testpair';
                    ttoind = find(true(size(testt)));
                case 4
                    whichrhoscorefield = 'similpair';
                    ttoind = find(~ismember(hireptt,testt));
            end
            subplot(2,2,isp)
            hold all
            pl = plot(disperses, squeeze(mean(rhoxneudivlmlvs.(whichrhoscorefield).avg(:,ttoind,:),1) ) );
            yl = ylim;
            for ii = 1:numel(pl)
                if ismember(hireptt(ttoind(ii)), [1105 1109])
                    lw = 2;
                    disp(ii)
                else
                    lw = 0.2;
                end
                plot(disperses, squeeze(mean(rhoxneudivlmlvs.(whichrhoscorefield).avg(:,ttoind(ii),:),1) ), 'Color', pl(ii).Color, 'LineWidth', lw)
                h = squeeze(mean(rhoxneudivasis.(whichrhoscorefield).avg(:,ttoind(ii) ),1));
                plot([disperses(1) disperses(end)], h*[1 1], 'Color', pl(ii).Color, 'LineWidth', lw)
                text(disperses(1), yl(1)+(numel(pl)-ii-1)*0.08*range(yl), num2str(hireptt(ttoind(ii))), 'Color', pl(ii).Color, 'FontSize', 10, 'VerticalAlignment', 'bottom')
            end
            plot(disperses, squeeze(mean(rhoxneudivlmlvs.(whichrhoscorefield).avg(:,ttoind,:),[1,2])), 'k-', 'LineWidth',2 );
            ylim(yl)
            xlabel('LMLV slopes')
            ylabel('% rho(SVM score)==1')
            title(whichrhoscorefield)
        end
    end
    
end
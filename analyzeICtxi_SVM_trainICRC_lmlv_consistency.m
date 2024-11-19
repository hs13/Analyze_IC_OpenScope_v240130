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

rng(111)
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
            optexclneu = 'incl';
        case 1
            optexclneu = 'excltt';
        case 2
            optexclneu = 'excl';
        otherwise
            error('excludeneuvar0 option not recognized')
    end
    
    svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_', optexclneu, '.mat');
    
    svmfndivs = cell(1,2);
    for d = 1:2
        svmfndivs{d} = strcat(pathpp, 'SVM', num2str(d), '_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_', optexclneu, '.mat');
        svmmdlfndivs{d} = strcat(pathpp, 'SVMmodels', num2str(d), '_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_', optexclneu, '.mat');
        svmlmlvfndivs{d} = strcat(pathpp, 'SVM', num2str(d), '_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_', optexclneu, '_lmlvslopes.mat');
        svmmdllmlvfndivs{d} = strcat(pathpp, 'SVMmodels', num2str(d), '_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_', optexclneu, '_lmlvslopes.mat');
    end
    consistfn = strcat(pathpp, 'consistency_SVM12_', svmdesc, '_lmlvslopes_', optexclneu, '.mat');
    
    load(svmfn)
    
    % parametrically change proportion silenced
    propneu2sil = 0.5;
    % randomly select neurons to silence: sample 100X
    Nneudivs = 1; % estimated ~30min for 100 samples
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
    
    SVMout1_lmlvs = struct();
    SVMmodels1_lmlvs = struct();
    SVMout2_lmlvs = struct();
    SVMmodels2_lmlvs = struct();
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
        
        neu2silvec = false(nnz(valneu), 1);
        neurand = randperm(nnz(valneu), round(propneu2sil*nnz(valneu)));
        neu2silvec(neurand) = true;

        cvtrials = struct();
        cvtrials.optimizeSVM = 2;
        cvtrials.loadcvpartition = true;
        cvtrials.traintrialinds = traintrialinds;
        cvtrials.testtrialinds = testtrialinds;
        % SVMtrainICRC_models is 14 MB
        [SVMout1, SVMmodels1] = computeICtxi_SVM(tempR, trialorder, ...
            svmdesc, 'spkcnt', preproc, whichSVMkernel, cvtrials);
        
        SVMout1.neu2silvec = neu2silvec;
        if islope>0 %&& excludeneuvar0>0
            SVMout1.valneu = valneu;
        end
        
        
        
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
                Xsilrand(:,neu2silvec) = 0;
                [templabel,tempscore] = predict(SVM_models{isplit}, Xsilrand);

                Xsilcompl = Tp;
                Xsilcompl(:,~neu2silvec) = 0;
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
           
        
        if islope==0
            save(svmfndivs{1}, 'preproc', 'whichSVMkernel', 'SVMout1', '-v7.3')
            save(svmmdlfndivs{1}, 'preproc', 'whichSVMkernel', 'SVMmodels1', '-v7.3')
            save(svmfndivs{2}, 'preproc', 'whichSVMkernel', 'SVMout2', '-v7.3')
            save(svmmdlfndivs{2}, 'preproc', 'whichSVMkernel', 'SVMmodels2', '-v7.3')
            fprintf('silence %.2f, as-is done\n', propneu2sil)
        else
            if islope==1
                SVMout1_lmlvs = SVMout1;
                SVMmodels1_lmlvs = SVMmodels1;
                SVMout2_lmlvs = SVMout2;
                SVMmodels2_lmlvs = SVMmodels2;
            else
                SVMout1_lmlvs(islope) = SVMout1;
                SVMmodels1_lmlvs(islope) = SVMmodels1;
                SVMout2_lmlvs(islope) = SVMout2;
                SVMmodels2_lmlvs(islope) = SVMmodels2;
            end
            fprintf('silence %.2f, lmlv slope %.2f done\n', propneu2sil, disperses(islope))
        end
        
    end
    if islope==numel(disperses)
        save(svmlmlvfndivs{1}, 'disperses', 'preproc', 'whichSVMkernel', 'SVMout1_lmlvs', '-v7.3')
        save(svmmdllmlvfndivs{1}, 'disperses', 'preproc', 'whichSVMkernel', 'SVMmodels1_lmlvs', '-v7.3')
        save(svmlmlvfndivs{2}, 'disperses', 'preproc', 'whichSVMkernel', 'SVMout2_lmlvs', '-v7.3')
        save(svmmdllmlvfndivs{2}, 'disperses', 'preproc', 'whichSVMkernel', 'SVMmodels2_lmlvs', '-v7.3')
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
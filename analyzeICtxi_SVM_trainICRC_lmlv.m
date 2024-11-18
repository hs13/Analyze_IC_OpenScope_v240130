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

% ises = 3;
% pathpp= [drivepath 'RESEARCH/logmean_logvar/' nwbsessions{ises} '/'];
% % pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-619293/'];
% % pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-620333/'];

for ises = 1:numel(nwbsessions)
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
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes_incl.mat');
        case 1
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_excltt_lmlvslopes.mat');
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes.mat_excltt');
        case 2
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmlmlvfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            svmmdllmlvfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes_excl.mat');
        otherwise
            error('excludeneuvar0 option not recognized')
    end

    % if exist(silrandfn, 'file')
    %     fprintf('excludeneuvar0 version %d random silencing SVM already calculated\n', excludeneuvar0)
    %     continue
    % end

    load(svmfn)
    load(svmmdlfn)
    load(svmlmlvfn)
    load(svmmdllmlvfn)

    testt = [106,107,110,111];
    inferencett = [1105 1109];
    Nhireptt = numel(hireptt);
    Ntt = numel(testt);
    Nsplits = size(SVMtrainICRC.spkcnt.testtrialinds,2);
    Nneurons = SVMtrainICRC.Nneurons;

    traincmlmlvs = zeros(Ntt, Ntt, numel(disperses));
    testcmlmlvs = zeros(Ntt, Ntt, numel(disperses));
    infcmlmlvs = zeros(numel(inferencett), Ntt, numel(disperses));
    trainperflmlvs = zeros(1, numel(disperses));
    testperflmlvs = zeros(1, numel(disperses));
    infperflmlvs = zeros(1, numel(disperses));

    for islope = 0:numel(disperses)
        if islope==0
            SVMout = SVMtrainICRC;
            SVM_models = SVMtrainICRC_models.spkcnt;
            valneu = true(Nneurons,1);
        else
            SVMout = SVMtrainICRC_lmlvs(islope);
            SVM_models = SVMtrainICRC_models_lmlvs(islope).spkcnt;
            valneu = SVMout.valneu;
        end
        for isplit = 1:Nsplits
            if ~isequal(unique(SVMout.spkcnt.all.label(:,isplit)), testt')
                warning('check trained SVM: not all trial types are respresented')
            end
        end

        % train accuracy
        trainlabs = SVMout.trialorder(SVMout.spkcnt.traintrialinds);
        trainpred = SVMout.spkcnt.train.label;
        trainacc = zeros(Ntt);
        for itt = 1:Ntt
            trialsoi = trainlabs==testt(itt);
            [v,c]=uniquecnt(trainpred(trialsoi));
            trainacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
        end

        % test accuracy
        testlabs = SVMout.trialorder(SVMout.spkcnt.testtrialinds);
        testpred = SVMout.spkcnt.test.label;
        testacc = zeros(Ntt);
        for itt = 1:Ntt
            trialsoi = testlabs==testt(itt);
            [v,c]=uniquecnt(testpred(trialsoi));
            testacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
        end

        % inference decoding
        inftrials = ismember(SVMout.trialorder, inferencett);
        infpred = SVMout.spkcnt.all.label(inftrials,:);
        infperf = zeros(numel(inferencett), Ntt);
        for itt = 1:numel(inferencett)
            trialsoi = SVMout.trialorder(inftrials)==inferencett(itt);
            [v,c]=uniquecnt(infpred(trialsoi,:));
            infperf(itt, ismember(testt,v)) = c/(Nsplits*nnz(trialsoi));
        end

        if islope==0
            traincmasis = trainacc;
            testcmasis = testacc;
            infcmasis = infperf;

            trainperfasis = mean(diag(trainacc));
            testperfasis = mean(diag(testacc));
            infperfasis = (infperf(1,1)-infperf(1,2)-infperf(2,3)+infperf(2,4))/2;
        else
            traincmlmlvs(:,:,islope) = trainacc;
            testcmlmlvs(:,:,islope) = testacc;
            infcmlmlvs(:,:,islope) = infperf;

            trainperflmlvs(islope) = mean(diag(trainacc));
            testperflmlvs(islope) = mean(diag(testacc));
            infperflmlvs(islope) = (infperf(1,1)-infperf(1,2)-infperf(2,3)+infperf(2,4))/2;
        end
    end

    if pltses
        figure;
        subplot(1,2,1)
        hold all
        plot(disperses, testperflmlvs)
        xl = [disperses(1) disperses(end)];
        plot(xl, testperfasis*[1 1], 'c-')
        xlabel('log(mean) vs log(var) slopes')
        ylabel('test accuracy')
        subplot(1,2,2)
        hold all
        plot(disperses, infperflmlvs)
        xl = [disperses(1) disperses(end)];
        plot(xl, infperfasis*[1 1], 'c-')
        xlabel('log(mean) vs log(var) slopes')
        ylabel('inference performance')
    end

    %% to test robustness silence a portion of neurons and see change in classification accuracy across slopes
    % parametrically change proportion silenced
    propneusilvec = 0:0.1:1;
    % randomly select neurons to silence: sample 100X
    Nsamples = 10; % estimated ~30min for 100 samples
    Nsplits = size(SVMtrainICRC.spkcnt.testtrialinds,2);

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
    silrandtraincmasis = NaN(Nsamples, numel(propneusilvec), Ntt, Ntt);
    silrandtestcmasis = NaN(Nsamples, numel(propneusilvec), Ntt, Ntt);
    silrandinfcmasis = NaN(Nsamples, numel(propneusilvec), numel(inferencett), Ntt);
    silrandtrainperfasis = NaN(Nsamples, numel(propneusilvec));
    silrandtestperfasis = NaN(Nsamples, numel(propneusilvec));
    silrandinfperfasis = NaN(Nsamples, numel(propneusilvec));
    silrandtraincmlmlvs = cell(size(disperses));
    silrandtestcmlmlvs = cell(size(disperses));
    silrandinfcmlmlvs = cell(size(disperses));
    silrandtrainperflmlvs = cell(size(disperses));
    silrandtestperflmlvs = cell(size(disperses));
    silrandinfperflmlvs = cell(size(disperses));
    for islope = 1:numel(disperses)
        silrandtraincmlmlvs{islope} = NaN(Nsamples, numel(propneusilvec), Ntt, Ntt);
        silrandtestcmlmlvs{islope} = NaN(Nsamples, numel(propneusilvec), Ntt, Ntt);
        silrandinfcmlmlvs{islope} = NaN(Nsamples, numel(propneusilvec), numel(inferencett), Ntt);
        silrandtrainperflmlvs{islope} = NaN(Nsamples, numel(propneusilvec));
        silrandtestperflmlvs{islope} = NaN(Nsamples, numel(propneusilvec));
        silrandinfperflmlvs{islope} = NaN(Nsamples, numel(propneusilvec));
    end
    silrandclk = tic;
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

            trainaccpd = NaN(Nsamples, Ntt,Ntt, Nsplits);
            testaccpd = NaN(Nsamples, Ntt,Ntt, Nsplits);
            infaccpd = NaN(Nsamples, numel(inferencett),Ntt, Nsplits);
            trainperfpd = NaN(Nsamples, Nsplits);
            testperfpd = NaN(Nsamples, Nsplits);
            infperfpd = NaN(Nsamples, Nsplits);
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
                Tp(isnan(Tp))=0;

                trainlabs = SVMout.trialorder(traintrialinds);
                testlabs = SVMout.trialorder(testtrialinds);
                inftrials = ismember(SVMout.trialorder, inferencett);

                for s = 1:size(neu2silmat,2) % Nsamples if 0<propneu2sil<1
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

                    % train accuracy
                    trainpred = templabel(traintrialinds);
                    trainacc = zeros(Ntt);
                    for itt = 1:Ntt
                        trialsoi = trainlabs==testt(itt);
                        [v,c]=uniquecnt(trainpred(trialsoi));
                        trainacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
                    end

                    % test accuracy
                    testpred = templabel(testtrialinds);
                    testacc = zeros(Ntt);
                    for itt = 1:Ntt
                        trialsoi = testlabs==testt(itt);
                        [v,c]=uniquecnt(testpred(trialsoi));
                        testacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
                    end

                    % inference decoding
                    infpred = templabel(inftrials);
                    infperf = zeros(numel(inferencett), Ntt);
                    for itt = 1:numel(inferencett)
                        trialsoi = SVMout.trialorder(inftrials)==inferencett(itt);
                        [v,c]=uniquecnt(infpred(trialsoi));
                        infperf(itt, ismember(testt,v)) = c/nnz(trialsoi);
                    end

                    trainaccpd(s,:,:,isplit) = trainacc;
                    testaccpd(s,:,:,isplit) = testacc;
                    infaccpd(s,:,:,isplit) = infperf;
                    trainperfpd(s,isplit) = mean(diag(trainacc));
                    testperfpd(s,isplit) = mean(diag(testacc));
                    infperfpd(s,isplit) = (infperf(1,1)-infperf(1,2)-infperf(2,3)+infperf(2,4))/2;
                end
            end

            if islope==0
                silrandtraincmasis(:,n,:,:) = mean(trainaccpd,4);
                silrandtestcmasis(:,n,:,:) = mean(testaccpd,4);
                silrandinfcmasis(:,n,:,:) = mean(infaccpd,4);

                silrandtrainperfasis(:,n) = mean(trainperfpd,2);
                silrandtestperfasis(:,n) = mean(testperfpd,2);
                silrandinfperfasis(:,n) = mean(infperfpd,2);
                fprintf('silence %.2f, as-is done\n', propneu2sil)
            else
                silrandtraincmlmlvs{islope}(:,n,:,:) = mean(trainaccpd,4);
                silrandtestcmlmlvs{islope}(:,n,:,:) = mean(testaccpd,4);
                silrandinfcmlmlvs{islope}(:,n,:,:) = mean(infaccpd,4);

                silrandtrainperflmlvs{islope}(:,n) = mean(trainperfpd,2);
                silrandtestperflmlvs{islope}(:,n) = mean(testperfpd,2);
                silrandinfperflmlvs{islope}(:,n) = mean(infperfpd,2);
                fprintf('silence %.2f, lmlv slope %.2f done\n', propneu2sil, disperses(islope))
            end
            toc
        end
    end
    toc(silrandclk)

    if pltses
        lmlvcols = cool(numel(disperses));
        lmlvlegs = cell(1,numel(disperses)+1);
        figure
        hold all
        for islope = 1:numel(disperses)
            if disperses(islope)==1
                lw = 2;
            else
                lw = 1;
            end
            errorbar( 100*(1-propneusilvec), nanmean(silrandtestperflmlvs{islope},1), ...
                nanstd(silrandtestperflmlvs{islope},0,1)./sqrt(sum(~isnan(silrandtestperflmlvs{islope}),1)), ...
                'Color', lmlvcols(islope,:), 'LineWidth', lw)
            lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
        end
        errorbar( 100*(1-propneusilvec), nanmean(silrandtestperfasis,1), ...
            nanstd(silrandtestperfasis,0,1)./sqrt(sum(~isnan(silrandtestperfasis),1)), ...
            'Color', 'k', 'LineWidth', 2)
        lmlvlegs{end} = 'as-is';
        legend(lmlvlegs)
        xlabel('%Neurons')
        ylabel('test accuracy')
        % silrandinfperflmlvs{islope}


        lmlvcols = cool(numel(disperses));
        lmlvlegs = cell(1,numel(disperses)+1);
        figure
        hold all
        for islope = 1:numel(disperses)
            if disperses(islope)==1
                lw = 2;
            else
                lw = 1;
            end
            errorbar( 100*(1-propneusilvec), nanmean(silrandinfperflmlvs{islope},1), ...
                nanstd(silrandinfperflmlvs{islope},0,1)./sqrt(sum(~isnan(silrandinfperflmlvs{islope}),1)), ...
                'Color', lmlvcols(islope,:), 'LineWidth', lw)
            lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
        end
        errorbar( 100*(1-propneusilvec), nanmean(silrandinfperfasis,1), ...
            nanstd(silrandinfperfasis,0,1)./sqrt(sum(~isnan(silrandinfperfasis),1)), ...
            'Color', 'k', 'LineWidth', 2)
        lmlvlegs{end} = 'as-is';
        legend(lmlvlegs)
        xlabel('%Neurons')
        ylabel('inference performance')
        % silrandinfperflmlvs{islope}

        lmlvcols = cool(numel(disperses));
        lmlvlegs = cell(1,numel(disperses)+1);
        figure
        for isp = 1:4
            switch isp
                case 1
                    ylab = 'P(IC1|TRE1)';
                    r=1; c=1;
                case 2
                    ylab = 'P(LC1|TRE1)';
                    r=1; c=2;
                case 3
                    ylab = 'P(LC2|TRE2)';
                    r=2; c=3;
                case 4
                    ylab = 'P(IC2|TRE2)';
                    r=2; c=4;
            end
            subplot(2,2,isp)
            hold all
            for islope = 1:numel(disperses)
                if disperses(islope)==1
                    lw = 2;
                else
                    lw = 1;
                end
                errorbar( 100*(1-propneusilvec), nanmean(silrandinfcmlmlvs{islope}(:,:,r,c),1), ...
                    nanstd(silrandinfcmlmlvs{islope}(:,:,r,c),0,1)./sqrt(sum(~isnan(silrandinfcmlmlvs{islope}(:,:,r,c)),1)), ...
                    'Color', lmlvcols(islope,:) , 'LineWidth', lw)
                lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
            end
            errorbar( 100*(1-propneusilvec), nanmean(silrandinfcmasis(:,:,r,c),1), ...
                nanstd(silrandinfcmasis(:,:,r,c),0,1)./sqrt(sum(~isnan(silrandinfcmasis(:,:,r,c)),1)), ...
                'Color', 'k', 'LineWidth', 2)
            lmlvlegs{end} = 'as-is';
            legend(lmlvlegs)
            xlabel('%Neurons')
            ylabel(ylab)
        end

    end

    %%
    save(silrandfn, 'excludeneuvar0', 'traincmasis', 'testcmasis', 'infcmasis', ...
        'trainperfasis', 'testperfasis', 'infperfasis', ...
        'traincmlmlvs', 'testcmlmlvs', 'infcmlmlvs', ...
        'trainperflmlvs', 'testperflmlvs', 'infperflmlvs', ...
        'disperses', 'propneusilvec', 'silrandtraincmasis', 'silrandtestcmasis', 'silrandinfcmasis', ...
        'silrandtrainperfasis', 'silrandtestperfasis', 'silrandinfperfasis', ...
        'silrandtraincmlmlvs', 'silrandtestcmlmlvs', 'silrandinfcmlmlvs', ...
        'silrandtrainperflmlvs', 'silrandtestperflmlvs', 'silrandinfperflmlvs', '-v7.3')

    fprintf('%d/%d %s done silencing random subsets of neurons in SVM across LMLV slopes\n', ises, numel(nwbsessions), nwbsessions{ises})
    toc(sesclk)
end

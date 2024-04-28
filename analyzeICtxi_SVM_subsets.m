% clear all; close all; clc

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

% svmdesc = 'trainICRC';
preproc = 'zscore'; % zscore or minmax or meancenter
whichSVMkernel = 'Linear';
neuopt = 'RS';

% takes ~50hrs per session...
for ises = 1:Nsessions
    clearvars -except datadir nwbsessions Nsessions ises whichSVMkernel svmdesc preproc neuopt
    sesclk = tic;

    justctx = true;
    computeSVM = true;
    optimizeSVM = true;
    allblocksareas = false;
    loadVISpCV = true;

    mousedate = nwbsessions{ises};
    fprintf(strcat('%d  ', mousedate, '\n'), ises)

    pathpp = [datadir 'postprocessed' filesep mousedate filesep];
    pathsvm = [datadir 'SVM_' svmdesc '_subsets' filesep mousedate filesep];
    if ~exist(pathsvm, 'dir')
        mkdir(pathsvm)
    end
    load(sprintf('%spostprocessed.mat', pathpp ))
    load(sprintf('%svisresponses.mat', pathpp ))
    load(sprintf('%sqc_units.mat', pathpp ))
    neuRS = unit_wfdur>0.4;
    neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
    switch neuopt
        case ''
            neucrit = true(size(neuallloc));
        case 'RS'
            neucrit = neuRS;
        case 'filtRS'
            neucrit = neufilt & neuRS;
        otherwise
            error('neuron criterion not recognized')
    end

    if allblocksareas
        areas2anal = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
        ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
    else
        areas2anal = {'VISp'};
        ICblocks = {'ICwcfg1_presentations'};
    end

    switch svmdesc
        case 'trainICRC'
            traintrialtypes = [106, 107, 110, 111];
            probetrialtypes = [1105, 1109];
        case 'trainREx'
            traintrialtypes = [1201, 1299];
            probetrialtypes = [106, 107, 110, 111];
        otherwise
            error(['set ' svmdesc])
    end
    ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
        1301 1302 1303 1304 1305 1306 1307 1308];

    whichR = 'spkcnt';
    Nsplits = 10;

    commonanaltrials = struct();
    commontraintrialinds = struct();
    commontesttrialinds = struct();
    for b = 1:numel(ICblocks)
        tempsplit = strsplit(ICblocks{b}, '_');
        whichICblock = tempsplit{1};

        trialorder = ICtrialtypes( vis.(ICblocks{b}).trialorder + 1);

        Ntt = numel(traintrialtypes);
        numtrialstt = zeros(Ntt,1);
        for typi1 = 1:Ntt
            numtrialstt(typi1) = nnz(trialorder==traintrialtypes(typi1));
        end

        Ntrialspertype = min(numtrialstt);
        Ntesttrialspertype = floor(Ntrialspertype/Nsplits);
        Ntraintrialspertype = Ntrialspertype - Ntesttrialspertype;
        Ntraintrials = Ntt*Ntraintrialspertype;
        Ntesttrials = Ntt*(Ntrialspertype-Ntraintrialspertype);

        if loadVISpCV
            whichvisarea = 'VISp';
            pathsvm0 = [datadir 'SVM_' svmdesc '_selectareas' filesep mousedate filesep];
            svmfn = strcat(pathsvm0, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
            load(svmfn)
            switch svmdesc
                case 'trainICRC'
                    SVMout = SVMtrainICRC;
                case 'trainREx'
                    SVMout = SVMtrainREx;
                otherwise
                    error(['set ' svmdesc])
            end
            commonanaltrials.(whichICblock) = SVMout.analtrials;
            if ~isequal(SVMout.analtriallabels, trialorder(commonanaltrials.(whichICblock)))
                error('check analtrials')
            end
            commontraintrialinds.(whichICblock) = SVMout.(whichR).traintrialinds;
            commontesttrialinds.(whichICblock) = SVMout.(whichR).testtrialinds;
        else

            % balance trials
            tttrials = ismember(trialorder, traintrialtypes);
            if all(SVMout.numtrials==Ntrialspertype)
                trials2anal = tttrials;
                commonanaltrials.(whichICblock) = find(trials2anal);
            else
                warning('balancing number of trials')
                trials2anal = false(numrectrials,1);
                for typi1 = 1:Ntt
                    trialsintype = find(trialorder==SVMout.trialtypes(typi1));
                    trialsintype = trialsintype(1:Ntrialspertype);
                    trials2anal(trialsintype) = true;
                end
                if all(tttrials(trials2anal)) && ~any(trials2anal(~tttrials))
                    commonanaltrials.(whichICblock) = find(trials2anal);
                else
                    error('trials to analyze was not selected correctly')
                end
            end

            trials2analind = find(trials2anal); % consider randomizing the order of this
            C = cvpartition(trialorder(trials2analind),'KFold',Nsplits, 'Stratify',true);
            if ~( all(C.TrainSize==Ntraintrials) && all(C.TestSize==Ntesttrials) )
                error('check balancing trials')
            end
            commontraintrialinds.(whichICblock) = zeros(Ntraintrials, Nsplits);
            commontesttrialinds.(whichICblock) = zeros(Ntesttrials, Nsplits);
            for isplit = 1:Nsplits

                idxTrain = training(C,isplit);
                traintrialinds = trials2analind(idxTrain);
                idxTest = test(C,isplit);
                testtrialinds = trials2analind(idxTest);

                if any(ismember(traintrialinds, testtrialinds))
                    error('train and test trials should not overlap')
                end
                if ~( all(ismember(trialorder(testtrialinds), SVMout.trialtypes)) && all(ismember(trialorder(traintrialinds), SVMout.trialtypes)) )
                    error('train and test trials of incorrect type detected')
                end

                commontraintrialinds.(whichICblock)(:,isplit) = traintrialinds;
                commontesttrialinds.(whichICblock)(:,isplit) = testtrialinds;
            end
        end
    end


    for a = 1:numel(areas2anal)
        % neu2anal = contains(neuallloc, areas2anal{a});
        whichvisarea = areas2anal{a};
        if strcmp(whichvisarea, 'VISp')
            neu2anal = neucrit & contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
        else
            neu2anal = neucrit & contains(neuallloc, whichvisarea);
        end
        if nnz(neu2anal)==0
            fprintf('%s %s has no units, skipping...\n', mousedate, areas2anal{a})
            continue
        end


        for b = 1:numel(ICblocks)
            tempsplit = strsplit(ICblocks{b}, '_');
            whichICblock = tempsplit{1};
            svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_subsets_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
            svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, neuopt, '_subsets_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
            if exist(svmfn, 'file')
                fprintf('%s already exists, skipping...\n', svmfn)
                continue
            end

            if size(Rall.(ICblocks{b}),2)~=length(neu2anal)
                error('check neu2anal')
            end
            tempR = Rall.(ICblocks{b})(:,neu2anal)';
            trialorder = ICtrialtypes( vis.(ICblocks{b}).trialorder + 1);

            Nneurons = nnz(neu2anal);
            numrectrials = size(tempR,2);

            %% determine neurons in subset
            tempneugroup = ICsigall.(ICblocks{b});
            tempneugroup.allneurons = true(size(ICsigall.(ICblocks{b}).PkwBK));
            tempneugroup.sigkwBI = ICsigall.(ICblocks{b}).PkwBI<0.05;
            tempneugroup.sigkwBK = ICsigall.(ICblocks{b}).PkwBK<0.05;
            tempneugroup.sigkwBKnotsigkwBI = tempneugroup.sigkwBK & ~tempneugroup.sigkwBI;

            tempneugroup.indin = ICsigall.(ICblocks{b}).indin1 | ICsigall.(ICblocks{b}).indin2 | ICsigall.(ICblocks{b}).indin3 | ICsigall.(ICblocks{b}).indin4;
            tempneugroup.ICRCencoder = ICsigall.(ICblocks{b}).ICencoder | ICsigall.(ICblocks{b}).RCencoder;

            tempneugroup.sigkwCRF = RFCIall.Pkw_rfclassic<0.05;
            tempneugroup.ctrCRF9 = RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
            tempneugroup.sigkwIRF = RFCIall.Pkw_rfinverse<0.05;
            tempneugroup.ctrIRF9 = RFCIall.pRFinverse<0.05 & RFCIall.RFindinverse==1;
            tempneugroup.exclctrCRF9 = RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1 & sum(RFCIall.pRrfclassic<0.05,2)==1;
            tempneugroup.ctrCRF9nonctrIRF9 = RFCIall.RFindclassic==1 & RFCIall.pRFclassic<0.05 & RFCIall.pRrfinverse(:,1)>=0.05;
            tempneugroup.faithCRF9 = false(size(tempneugroup.ICencoder));
            for irf = 1:9
                tempfaith = RFCIall.RFindclassic==irf & RFCIall.pRFclassic<0.05 & RFCIall.pRrfinverse(:,irf)>=0.05;
                tempneugroup.faithCRF9 = tempneugroup.faithCRF9 | tempfaith;
            end
            tempneugroup.exclCRF9 = false(size(tempneugroup.ICencoder));
            for irf = 1:9
                tempexcl = RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==irf & sum(RFCIall.pRrfclassic<0.05,2)==1;
                tempneugroup.exclCRF9 = tempneugroup.exclCRF9 | tempexcl;
            end
            tempneugroup.ctrCG9 = RFCIall.pRrfclassic(:,1)<0.05;

            % switch whichsubsets
            %     case 'CGIG'
            %         subsetdescs = {'sigkwCRF', 'ctrCRF9', 'exclctrCRF9', 'ctrIRF9', ...
            %             'allbutsigkwCRF', 'allbutctrCG9', 'exclCRF9'};
            %     case 'encoder'
            %         subsetdescs = {'ICRCencoder', 'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive', 'indin', ...
            %             'allbutICRCencoder'};
            %     case 'sigkw'
            %         subsetdescs = {'allneurons', 'sigkwBI', 'sigkwBK', 'allbutsigkwBI', 'allbutsigkwBK', 'sigkwBKnotsigkwBI'};
            %     otherwise
            %         error([whichsubsets ' not defined'])
            % end

            subsetdescs = {'ICRCencoder', 'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive', 'indin'};
            subsets2train = cell(size(subsetdescs));
            for isd = 1:numel(subsetdescs)
                if contains(subsetdescs{isd}, 'allbut')
                    allbutstrsplit = strsplit(subsetdescs{isd}, 'allbut');
                    if ~all(ismember(tempneugroup.(allbutstrsplit{end}), [0 1]))
                        error('cannot be converted to logical : check')
                    end
                    subsets2train{isd} = ~logical(tempneugroup.(allbutstrsplit{end})(neu2anal));
                else
                    if ~all(ismember(tempneugroup.(subsetdescs{isd}), [0 1]))
                        error(['cannot be converted to logical : check ' subsetdescs{isd}])
                    end
                    subsets2train{isd} = logical(tempneugroup.(subsetdescs{isd})(neu2anal));
                end
                if length(subsets2train{isd}) ~= Nneurons
                    error('check Nneurons')
                end
            end

            %% train SVM
            SVM_models = struct();
            SVMout = struct();

            SVMout.subsetdescs = subsetdescs;
            SVMout.subsets2train = subsets2train;

            SVMout.neu2anal = neu2anal;
            SVMout.Nneurons = Nneurons;
            SVMout.exptid = ICblocks{b};
            SVMout.ICtrialtypes = ICtrialtypes;
            SVMout.trialtypes = traintrialtypes;

            % probe trials
            probetrials = ismember(trialorder, probetrialtypes);
            SVMout.probetrials = find(probetrials);

            alltrials = true(size(trialorder));
            SVMout.alltrials = find(alltrials);

            SVMout.trialorder = trialorder;

            SVMout.analtrials = commonanaltrials.(whichICblock);
            SVMout.analtriallabels = trialorder(commonanaltrials.(whichICblock));

            SVMout.traintrialinds = commontraintrialinds.(whichICblock);
            SVMout.testtrialinds = commontesttrialinds.(whichICblock);

            for isd = 1:numel(subsetdescs)
                whichR = ['spkcnt_' subsetdescs{isd}];

                if nnz(subsets2train{isd})==0
                    fprintf('%s %s %s %s has no neurons, skipping...\n', ...
                        mousedate, whichvisarea, neuopt, subsetdescs{isd})
                    SVM_models.(whichR) = cell(1, Nsplits);
                    SVMout.(whichR) = struct();
                    continue
                end

                SVM_models.(whichR) = cell(1, Nsplits);

                for ts = 1:4
                    switch ts
                        case 1
                            svmmd = 'train';
                            tempNtrials = Ntraintrials;
                        case 2
                            svmmd = 'test';
                            tempNtrials = Ntesttrials;
                        case 3
                            svmmd = 'probe';
                            tempNtrials = nnz(probetrials);
                        case 4
                            svmmd = 'all';
                            tempNtrials = nnz(alltrials);
                    end
                    SVMout.(whichR).(svmmd).label = NaN(tempNtrials, Nsplits);
                    SVMout.(whichR).(svmmd).score = NaN(tempNtrials,Ntt, Nsplits);
                end

                for isplit = 1:Nsplits
                    ttclk = tic;

                    traintrialinds = commontraintrialinds.(whichICblock)(:,isplit);
                    testtrialinds = commontesttrialinds.(whichICblock)(:,isplit);

                    if any(ismember(traintrialinds, testtrialinds))
                        error('train and test trials should not overlap')
                    end
                    if ~( all(ismember(trialorder(testtrialinds), SVMout.trialtypes)) && all(ismember(trialorder(traintrialinds), SVMout.trialtypes)) )
                        error('train and test trials of incorrect type detected')
                    end

                    switch preproc
                        case 'none'
                            Tp = tempR';
                        case 'zscore'
                            % Z-score
                            trainRmean = mean(tempR(:,traintrialinds),2);
                            trainRstd = std(tempR(:,traintrialinds),0,2);

                            Tp = ( (tempR-trainRmean)./trainRstd )';
                            Tp(isnan(Tp))=0;
                        case 'minmax'
                            trainRmin = min(tempR(:,traintrialinds),[],2);
                            trainRrange = range(tempR(:,traintrialinds),2);

                            Tp = ( (tempR-trainRmin)./trainRrange )';
                        case 'meancenter'
                            trainRmean = mean(tempR(:,traintrialinds),2);

                            Tp = (tempR-trainRmean)';
                    end
                    Tp = Tp(:,subsets2train{isd}==1);

                    X = Tp(traintrialinds,:);
                    Y = trialorder(traintrialinds);

                    Xtest = Tp(testtrialinds,:);
                    Ytest = trialorder(testtrialinds);

                    Xprobe = Tp(probetrials,:);
                    Xall = Tp(alltrials,:);

                    % t is an SVM template. Most of its properties are empty.
                    % When the software trains the ECOC classifier, it sets the applicable properties to their default values.
                    % Train the ECOC classifier using the SVM template.
                    % Transform classification scores to class posterior probabilities
                    % (which are returned by predict or resubPredict) using the 'FitPosterior' name-value pair argument.
                    % Specify the class order using the 'ClassNames' name-value pair argument.
                    % Display diagnostic messages during training by using the 'Verbose' name-value pair argument.

                    Ylabs = unique(Y);
                    SVMout.(whichR).Ylabs = Ylabs;

                    switch whichSVMkernel
                        case 'RBF'
                            t = templateSVM('Standardize',true,'KernelFunction', 'rbf');
                        case 'Linear'
                            t = templateSVM('Standardize',true,'KernelFunction', 'linear');
                        case 'Poly2'
                            t = templateSVM('Standardize',true,'KernelFunction', 'polynomial' , 'PolynomialOrder', 2);
                    end
                    if optimizeSVM
                        SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, ...
                            'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', 'auto', ...
                            'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false));
                    else
                        SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, 'ClassNames', Ylabs, 'Verbose',0);
                    end

                    SVM_models.(whichR){isplit} = SVMModel;
                    for t = 1:4
                        switch t
                            case 1
                                Xtemp = X;
                                Ytemp = Y;
                                tempSVMmodel = SVMModel;
                                svmmd = 'train';
                            case 2
                                Xtemp = Xtest;
                                Ytemp = Ytest;
                                tempSVMmodel = SVMModel;
                                svmmd = 'test';
                            case 3
                                Xtemp = Xprobe;
                                tempSVMmodel = SVMModel;
                                svmmd = 'probe';
                            case 4
                                Xtemp = Xall;
                                tempSVMmodel = SVMModel;
                                svmmd = 'all';
                        end
                        [templabel,tempscore] = predict(tempSVMmodel,Xtemp);
                        SVMout.(whichR).(svmmd).label(:,isplit) = templabel;
                        SVMout.(whichR).(svmmd).score(:,:,isplit) = tempscore;
                    end

                    fprintf('%d/%d %s %s %s %s %s %d/%d\n', ises, Nsessions, mousedate, whichSVMkernel, preproc, whichvisarea, whichICblock, isplit, Nsplits)
                    toc(ttclk)
                end
            end

            switch svmdesc
                case 'trainICRC'
                    SVMtrainICRC = SVMout;
                    SVMtrainICRC_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models', '-v7.3')
                case 'trainREx'
                    SVMtrainREx = SVMout;
                    SVMtrainREx_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx_models', '-v7.3')
                otherwise
                    error(['set ' svmdesc])
            end

            fprintf('%d/%d %s %s %s %s %s done!\n', ises, Nsessions, mousedate, whichSVMkernel, preproc, whichvisarea, whichICblock)
            toc(sesclk)
        end
    end
end

disp(['analyzeICtx_SVM ' svmdesc ' trainsubsets FINISHED! READY TO MOVE DATA'])


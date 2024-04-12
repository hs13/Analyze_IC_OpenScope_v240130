% 0. score needs to be reordered?
% is there bias if I don't mask the order of Ylabs?
% 1. predict returns [label, NegLoss, PBScore] when there are >2 classes,
% but retnuns [label, score] when there are only two classes
% 2. before, trials2anal (balanced trials if trial repeat is different amongst train trial types) was ignored after being computed -- now train and test trials are only being drawn from trials2anal

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));

Nsessions = numel(nwbsessions);

% takes ~50hrs per session...
rng('shuffle')
for ises = 1:Nsessions
    clearvars -except datadir nwbsessions ises whichSVMkernel svmdesc preproc
    sesclk = tic;
    
    allblocks = true;
    computeSVM = true;
    optimizeSVM = true;
    computesilencesubsets = false;
    
    svmdesc = 'trainICRC';
    preproc = 'zscore'; % zscore or minmax or meancenter
    whichSVMkernel = 'Linear';
    
    Nsplits = 10;
    
    mousedate = nwbsessions{ises};
    fprintf(strcat('%d  ', mousedate, '\n'), ises)
    
    if issavio
        pathpp = [datadir mousedate filesep];
    else
        pathpp = [datadir 'postprocessed' filesep mousedate filesep];
    end
    pathsvm = [datadir 'SVM_' svmdesc '_selectareas' filesep mousedate filesep];
    if ~exist(pathsvm, 'dir')
        mkdir(pathsvm)
    end
    
    % using VISp2 instead of VISp2/3 to avoid error
    % areas2anal = {'LGd', 'LP', 'VISp2', 'VISp4', 'VISp5', 'VISp6'};
    areas2anal = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam', 'LGd', 'LP'};
    if allblocks
        ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
    else
        ICblocks = {'ICwcfg1_presentations'};
    end
    ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
        1301 1302 1303 1304 1305 1306 1307 1308];
    
    load(sprintf('%spostprocessed.mat', pathpp ))
    load(sprintf('%svisresponses.mat', pathpp ))
    load(sprintf('%sqc_units.mat', pathpp ))
    neuRS = unit_wfdur>0.4;
    neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);

    for a = 1:numel(areas2anal)
        % neu2anal = contains(neuallloc, areas2anal{a});
        whichvisarea = areas2anal{a};
        if strcmp(whichvisarea, 'VISp')
            neu2anal = neuRS & contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
        else
            neu2anal = neuRS & contains(neuallloc, whichvisarea);
        end
        if nnz(neu2anal)==0
            fprintf('%s %s has no units, skipping...\n', mousedate, areas2anal{a})
            continue
        end
        
        for b = 1:numel(ICblocks)
            tempsplit = strsplit(ICblocks{b}, '_');
            whichICblock = tempsplit{1};
            if computesilencesubsets
                svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
                svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
            else
                svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
                svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
            end
            if exist(svmfn, 'file')
                fprintf('%s already exists, skipping...\n', svmfn)
                continue
            end
            
            whichR = 'spkcnt';
            if size(Rall.(ICblocks{b}),2)~=length(neu2anal)
                error('check neu2anal')
            end
            tempR = Rall.(ICblocks{b})(:,neu2anal)';
            trialorder = ICtrialtypes( vis.(ICblocks{b}).trialorder + 1);
            
            Nneurons = nnz(neu2anal);
            numrectrials = size(tempR,2);
            
            %% determine neurons to silence
            if computesilencesubsets
                tempneugroup = ICsigall.(ICblocks{b});
                tempneugroup.sigkwBI = ICsigall.(ICblocks{b}).PkwBI<0.05;
                tempneugroup.sigkwBK = ICsigall.(ICblocks{b}).PkwBK<0.05;
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
                
                silencedescs = {'sigkwBI', 'sigkwBK', 'ICRCencoder', 'ICencoder', 'RCencoder', ...
                    'inducerencoder', 'inducerresponsive', 'indin', ...
                    'sigkwCRF', 'ctrCRF9', 'exclctrCRF9', 'ctrCG9', 'exclCRF9', 'ctrIRF9', ...
                    'allbutsigkwBI', 'allbutsigkwBK', 'allbutICRCencoder', 'allbutICencoder', 'allbutRCencoder', ...
                    'allbutinducerencoder', 'allbutinducerresponsive', 'allbutindin', ...
                    'allbutsigkwCRF', 'allbutctrCRF9', 'allbutexclctrCRF9', 'allbutctrCG9', 'allbutexclCRF9', 'allbutctrIRF9'};
                subsets2silence = cell(size(silencedescs));
                for isd = 1:numel(silencedescs)
                    if contains(silencedescs{isd}, 'allbut')
                        allbutstrsplit = strsplit(silencedescs{isd}, 'allbut');
                        subsets2silence{isd} = ~tempneugroup.(allbutstrsplit{end})(neu2anal);
                    else
                        subsets2silence{isd} = tempneugroup.(silencedescs{isd})(neu2anal);
                    end
                    if length(subsets2silence{isd}) ~= Nneurons
                        error('check Nneurons')
                    end
                end
                % disp(sum(cat(2, subsets2silence{:}), 1))
            end
            %% train SVM
            SVM_models = struct();
            SVMout = struct();
            
            if computesilencesubsets
                SVMout.silencedescs = silencedescs;
                SVMout.subsets2silence = subsets2silence;
            end
            
            alltrialtypes = unique(trialorder);
            switch svmdesc
                case 'trainICRC'
                    traintrialtypes = [106, 107, 110, 111];
                    probetrialtypes = [1105, 1109];
                case 'trainREx'
                    traintrialtypes = [1201, 1299];
                    probetrialtypes = [106, 107, 110, 111];
                case 'trainIC1RC1'
                    traintrialtypes = [106, 107];
                    probetrialtypes = [1105];
                case 'trainIC2RC2'
                    traintrialtypes = [111 110];
                    probetrialtypes = [1109];
            end
            
            Ntt = numel(traintrialtypes);
            Nprobett = numel(probetrialtypes);
            Nalltt = numel(alltrialtypes);
            
            SVMout.neu2anal = neu2anal;
            SVMout.Nneurons = Nneurons;
            SVMout.exptid = ICblocks{b};
            SVMout.ICtrialtypes = ICtrialtypes;
            SVMout.trialtypes = traintrialtypes;
            
            SVMout.numtrials = zeros(Ntt,1);
            % DI.numtrialpairs = zeros(Ntt,Ntt);
            for typi1 = 1:Ntt
                SVMout.numtrials(typi1) = nnz(trialorder==SVMout.trialtypes(typi1));
            end
            
            % balance trials
            cftttrials = ismember(trialorder, SVMout.trialtypes);
            Ntrialspertype = min(SVMout.numtrials);
            if all(SVMout.numtrials==Ntrialspertype)
                trials2anal = cftttrials;
                SVMout.analtrials = find(trials2anal);
            else
                warning('balancing number of trials')
                trials2anal = false(numrectrials,1);
                for typi1 = 1:Ntt
                    trialsintype = find(trialorder==SVMout.trialtypes(typi1));
                    trialsintype = trialsintype(1:Ntrialspertype);
                    trials2anal(trialsintype) = true;
                end
                if all(cftttrials(trials2anal)) && ~any(trials2anal(~cftttrials))
                    SVMout.analtrials = find(trials2anal);
                else
                    error('trials to analyze was not selected correctly')
                end
            end
            SVMout.analtriallabels = trialorder(trials2anal);
            
            Ntesttrialspertype = floor(Ntrialspertype/Nsplits);
            Ntraintrialspertype = Ntrialspertype - Ntesttrialspertype;
            SVMout.Ntt = Ntt;
            SVMout.Ntrialspertype = Ntrialspertype;
            SVMout.Ntraintrialspertype = Ntraintrialspertype;
            
            Ntraintrials = Ntt*Ntraintrialspertype;
            Ntesttrials = Ntt*(Ntrialspertype-Ntraintrialspertype);
            
            % probe trials
            probetrials = ismember(trialorder, probetrialtypes);
            SVMout.probetrials = find(probetrials);
            
            alltrials = true(size(trialorder));
            SVMout.alltrials = find(alltrials);
            
            SVMout.trialorder = trialorder;
            
            % randtrialorder=randperm(numrectrials);
            % SVMout.randtrialorder = randtrialorder;
            % trials2anal = randtrialorder(ismember(randtrialorder, SVMout.analtrials));
            
            SVM_models.(whichR) = cell(1, Nsplits);
            
            SVMout.(whichR).traintrialinds = zeros(Ntraintrials, Nsplits);
            SVMout.(whichR).testtrialinds = zeros(Ntesttrials, Nsplits);
            
            % SVMout.(whichR).Ylabs = cell(Ntt, Nsplits);
            
            if computesilencesubsets
                for isd = 0:numel(silencedescs)
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
                        if isd>0
                            svmmd = [svmmd '_' silencedescs{isd}];
                        end
                        SVMout.(whichR).(svmmd).label = NaN(tempNtrials, Nsplits);
                        SVMout.(whichR).(svmmd).score = NaN(tempNtrials,Ntt, Nsplits);
                    end
                end
            end

            % Nsplits-fold cross-validation
            trials2analind = find(trials2anal); % consider randomizing the order of this
            C = cvpartition(trialorder(trials2analind),'KFold',Nsplits, 'Stratify',true);
            if ~( all(C.TrainSize==Ntraintrials) && all(C.TestSize==Ntesttrials) )
                error('check balancing trials')
            end
            for isplit = 1:Nsplits
                close all
                ttclk = tic;
                
                %{
                testtrialinds = zeros(Ntesttrials,1);
                traintrialinds = zeros(Ntraintrials,1);
                for typi1 = 1:Ntt
                    tempinds = trials2anal( trialorder(trials2anal)==SVMout.trialtypes(typi1) );
                    tempinds = reshape(tempinds,[],1);
                    if size(tempinds,1) ~= Ntrialspertype
                        error('Ntrialspertype not consistent between trial types? check')
                    end
                    temptestintype = false(Ntrialspertype,1);
                    temptestintype((isplit-1)*Ntesttrialspertype+1:isplit*Ntesttrialspertype) = true;
                    temptrainintype = true(Ntrialspertype,1);
                    temptrainintype((isplit-1)*Ntesttrialspertype+1:isplit*Ntesttrialspertype) = false;
                    testtrialinds((typi1-1)*Ntesttrialspertype+1:typi1*Ntesttrialspertype) = tempinds(temptestintype);
                    traintrialinds((typi1-1)*Ntraintrialspertype+1:typi1*Ntraintrialspertype) = tempinds(temptrainintype);
                end
                testtrialinds = trials2anal(ismember(trials2anal, testtrialinds));
                traintrialinds = trials2anal(ismember(trials2anal, traintrialinds));
                %}

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
                
                SVMout.(whichR).traintrialinds(:,isplit) = traintrialinds;
                SVMout.(whichR).testtrialinds(:,isplit) = testtrialinds;
                
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
                
                X = Tp(traintrialinds,:);
                Y = reshape(trialorder(traintrialinds),[],1);
                
                %                 X = X(randomizedtraintrialorder, :);
                %                 Y = Y(randomizedtraintrialorder);
                
                Xtest = Tp(testtrialinds,:);
                Ytest = reshape(trialorder(testtrialinds),[],1);
                
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
                %                 CVMdl = crossval(SVMModel);
                
                
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
                
                if computesilencesubsets
                    for isd = 1:numel(silencedescs)
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
                            svmmd = [svmmd '_' silencedescs{isd}];
                            Xtemp(:,subsets2silence{isd})=0;
                            
                            [templabel,tempscore] = predict(tempSVMmodel,Xtemp);
                            SVMout.(whichR).(svmmd).label(:,isplit) = templabel;
                            SVMout.(whichR).(svmmd).score(:,:,isplit) = tempscore;
                        end
                    end
                end
                
                fprintf('%s %s %s %s %d/%d\n', mousedate, whichSVMkernel, whichvisarea, whichICblock, isplit, Nsplits)
                toc(ttclk)
            end
            
            
            switch svmdesc
                case 'trainICRCtestRE'
                    SVMtrainICRC = SVMout;
                    SVMtrainICRC_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models', '-v7.3')
                case 'trainRExtestICRC'
                    SVMtrainREx = SVMout;
                    SVMtrainREx_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx_models', '-v7.3')
                case 'trainIC1RC1'
                    SVMtrainIC1RC1 = SVMout;
                    SVMtrainIC1RC1_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC1RC1', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC1RC1_models', '-v7.3')
                case 'trainIC2RC2'
                    SVMtrainIC2RC2 = SVMout;
                    SVMtrainIC2RC2_models = SVM_models;
                    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC2RC2', '-v7.3')
                    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC2RC2_models', '-v7.3')
            end
            
            fprintf('%s %s %s %s done!\n', mousedate, whichSVMkernel, whichvisarea, whichICblock)
            toc(sesclk)
        end
    end
end

disp(['analyzeICtx_SVM ' svmdesc ' ' whichSVMkernel ' FINISHED! READY TO MOVE DATA'])


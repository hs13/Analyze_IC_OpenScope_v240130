function [SVMout, SVM_models] = computeICtxi_SVM(tempR, trialorder, svmdesc, whichR, preproc, whichSVMkernel, cvtrials)
sesclk = tic;

optimizeSVM = 2; % 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
if cvtrials.loadcvpartition
    Nsplits = size(cvtrials.testtrialinds,2);
else
    Nsplits = cvtrials.Nsplits;
end

% whichR = 'spkcnt';
% if size(Rall.(ICblocks{b}),2)~=length(neu2anal)
%     error('check neu2anal')
% end
% tempR = Rall.(ICblocks{b})(:,neu2anal)';
% trialorder = ICtrialtypes( vis.(ICblocks{b}).trialorder + 1);

Nneurons = size(tempR,1);
numrectrials = size(tempR,2);

%% train SVM
SVM_models = struct();
SVMout = struct();

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
    case 'trainBK'
        traintrialtypes = [0, 106, 107, 110, 111];
        probetrialtypes = [1105, 1109];
    otherwise
        error(['set ' svmdesc])
end

Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);
Nalltt = numel(alltrialtypes);

SVMout.optimizeSVM = optimizeSVM;
SVMout.preproc = preproc;
SVMout.whichSVMkernel = whichSVMkernel;
SVMout.Nneurons = Nneurons;
% SVMout.exptid = ICblocks{b};
% SVMout.ICtrialtypes = ICtrialtypes;
SVMout.trialtypes = traintrialtypes;

SVMout.numtrials = zeros(Ntt,1);
% DI.numtrialpairs = zeros(Ntt,Ntt);
for typi1 = 1:Ntt
    SVMout.numtrials(typi1) = nnz(trialorder==SVMout.trialtypes(typi1));
end

if cvtrials.loadcvpartition
    SVMout.analtrials = reshape( unique(cvtrials.traintrialinds(:)) ,1,[]);
    SVMout.analtriallabels = trialorder(SVMout.analtrials);
    Ntraintrials = size(cvtrials.traintrialinds,1);
    Ntesttrials = size(cvtrials.testtrialinds,1);
    Ntraintrialspertype = Ntraintrials/Ntt;
    Ntesttrialspertype = Ntesttrials/Ntt;
    Ntrialspertype = Ntraintrialspertype+Ntesttrialspertype;
else
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

    % Nsplits-fold cross-validation
    trials2analind = find(trials2anal); % consider randomizing the order of this

    Ntesttrialspertype = floor(Ntrialspertype/Nsplits);
    Ntraintrialspertype = Ntrialspertype - Ntesttrialspertype;

    Ntraintrials = Ntt*Ntraintrialspertype;
    Ntesttrials = Ntt*(Ntrialspertype-Ntraintrialspertype);

    C = cvpartition(trialorder(trials2analind),'KFold',Nsplits, 'Stratify',true);
    if ~( all(C.TrainSize==Ntraintrials) && all(C.TestSize==Ntesttrials) )
        error('check balancing trials')
    end
end

SVMout.Ntt = Ntt;
SVMout.Ntrialspertype = Ntrialspertype;
SVMout.Ntraintrialspertype = Ntraintrialspertype;

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
SVMout.(whichR).Xall = cell(1, Nsplits);

for isplit = 1:Nsplits
    close all
    ttclk = tic;
    
    if cvtrials.loadcvpartition
        traintrialinds = cvtrials.traintrialinds(:,isplit);
        testtrialinds = cvtrials.testtrialinds(:,isplit);
    else
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
        traintrialinds = reshape( trials2analind(idxTrain) ,[],1);
        idxTest = test(C,isplit);
        testtrialinds = reshape( trials2analind(idxTest) ,[],1);
    end

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
        case 'minmax'
            trainRmin = min(tempR(:,traintrialinds),[],2);
            trainRrange = range(tempR(:,traintrialinds),2);
            
            Tp = ( (tempR-trainRmin)./trainRrange )';
        case 'meancenter'
            trainRmean = mean(tempR(:,traintrialinds),2);
            
            Tp = (tempR-trainRmean)';
    end
    Tp(isnan(Tp))=0;
    
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
    
    % HS 241115: when [[ 'Standardize',true ]], data is z-scored. 
    % mean centering didn't actually work
    switch whichSVMkernel
        case 'RBF'
            t = templateSVM('KernelFunction', 'rbf'); % 'Standardize' is false by default
        case 'Linear'
            t = templateSVM('KernelFunction', 'linear');
        case 'Poly2'
            t = templateSVM('KernelFunction', 'polynomial' , 'PolynomialOrder', 2);
    end
    switch optimizeSVM
        case 0
        SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, 'ClassNames', Ylabs, 'Verbose',0);
        case 1
        SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, ...
            'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', 'auto', ...
            'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false, 'Verbose',0));
        case 2
        SVMModel = fitcecoc(X,Y,'Coding','onevsone', 'Learners',t,'FitPosterior',false, ...
            'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ...
            'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false, 'Verbose',0));
        case 3
        SVMModel = fitcecoc(X,Y,'Coding','onevsall', 'Learners',t,'FitPosterior',false, ...
            'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', {'BoxConstraint','KernelScale'}, ...
            'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false, 'Verbose',0));
    end
    % optimized hyperparameters: Learners = 'svm' (default) — {'BoxConstraint','KernelScale'}
    % note, hyperparameter optimization entails five-fold cross-validation
    % Find hyperparameters that minimize five-fold cross-validation loss by 
    % using automatic hyperparameter optimization. 
    % CVMdl = crossval(SVMModel);
    % isequal(CVMdl.Y, SVMModel.Y)
    
    
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
    SVMout.(whichR).Xall{isplit} = Xall;
    
    fprintf('SVM %s %s %d/%d\n', whichSVMkernel, preproc, isplit, Nsplits)
    toc(ttclk)
end


toc(sesclk)
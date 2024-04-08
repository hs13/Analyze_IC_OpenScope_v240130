% figure
% for irf = 1:numel(gaborbankresp)
%     subplot(2,2,irf)
%     imagesc( squeeze(gaborbankresp(irf).responses(1,:,:)) )
% end

gfbdir = 'C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130\visICtxiwcfg1\';
rfpixs = [26 53 79 106];
gaborrespagg = struct();
for irf = 1:length(rfpixs)
    if isempty(fieldnames(gaborrespagg))
        gaborrespagg = load([gfbdir 'ICwcfg1_gaborfiltbank_responses' num2str(rfpixs(irf)) '.mat']);
    else
        gaborresp = load([gfbdir 'ICwcfg1_gaborfiltbank_responses' num2str(rfpixs(irf)) '.mat']);
        gaborrespagg = cat(1,gaborrespagg,gaborresp);
    end
end


% gaborresp = load([gfbdir 'ICwcfg1_gaborfiltbank_responses106.mat']);
% gabornormresp = load([gfbdir 'ICwcfg1_gaborfiltbank_normresp106.mat']);
% figure; histogram(gaborresp.allresps(:))
% figure; histogram(gabornormresp.allresps(:))
% isequal(gaborresp.allresps, gabornormresp.allresps)

% image sizes (300, 480)
%%
delete(gcp('nocreate'))
c = parcluster;
c.NumWorkers = 12;
parpool(c)

gfbdir = 'C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130\visICtxiwcfg1\';
rfpixs = [26 53 79 106];

for irf = 1:length(rfpixs)
    clearvars -except gfbdir irf rfpixs
    Npixinterval = 20;
    Nsplits = 10;
    svmdesc = 'trainICRC';
    preproc = 'zscore';
    whichSVMkernel = 'Linear';
    optimizeSVM = true;

    gaborresp = load([gfbdir 'ICwcfg1_gaborfiltbank_responses' num2str(rfpixs(irf)) '.mat']);
    ICtrialtypes = str2num(gaborresp.imlabel);
    filtsz = size(gaborresp.allresps);
    Nneurons = filtsz(2)*(filtsz(3)/Npixinterval)*(filtsz(4)/Npixinterval);
    Nrep = Npixinterval*Npixinterval;
    Rgabor = NaN(numel(ICtrialtypes), Nneurons, Nrep);
    offsetrc = zeros(Nrep,2);
    cnt = 0;
    for ii = 1:Npixinterval
        for jj = 1:Npixinterval
            cnt = cnt+1;
            offsetrc(cnt,:) = [ii,jj];
            Rgabor(:,:,cnt) = reshape( gaborresp.allresps(:,:,ii:Npixinterval:end,jj:Npixinterval:end), numel(ICtrialtypes),Nneurons );
        end
    end

    trialorder = reshape( repmat(ICtrialtypes',Nrep,1) ,[],1);
    tempR = reshape( permute(Rgabor,[2,3,1]) ,Nneurons, Nrep*length(ICtrialtypes)); % Nneurons X Ntrials
    valneu = std(tempR,0,2)>0;
    tempR = tempR(valneu,:);

    svmfn = [gfbdir 'SVMtrainICRC_ICwcfg1_gaborfiltbank' num2str(rfpixs(irf)) '.mat'];
    SVMgabor = struct();

    switch svmdesc
        case 'trainICRC'
            traintrialtypes = [106, 107, 110, 111];
        case 'trainREx'
            traintrialtypes = [1201, 1299];
    end
    Ntt = numel(traintrialtypes);

    SVMgabor.npix16deg = gaborresp.npix16deg;
    SVMgabor.rfsize = gaborresp.receptive_field_size;
    SVMgabor.kernels = gaborresp.kernels;
    SVMgabor.svmdesc = svmdesc;
    SVMgabor.preproc = preproc;
    SVMgabor.whichSVMkernel = whichSVMkernel;
    SVMgabor.validneurons = valneu;

    SVMgabor.exptid = 'ICwcfg1';
    SVMgabor.ICtrialtypes = ICtrialtypes;
    SVMgabor.traintrialtypes = traintrialtypes;

    SVMgabor.numtrials = zeros(Ntt,1);
    for typi1 = 1:Ntt
        SVMgabor.numtrials(typi1) = nnz(trialorder==traintrialtypes(typi1));
    end

    % balance trials
    cftttrials = ismember(trialorder, traintrialtypes);
    Ntrialspertype = min(SVMgabor.numtrials);
    if all(SVMgabor.numtrials==Ntrialspertype)
        trials2anal = cftttrials;
        SVMgabor.analtrials = find(trials2anal);
    else
        warning('balancing number of trials')
        trials2anal = false(numtrials,1);
        for typi1 = 1:Ntt
            trialsintype = find(trialorder==traintrialtypes(typi1));
            trialsintype = trialsintype(1:Ntrialspertype);
            trials2anal(trialsintype) = true;
        end
        if all(cftttrials(trials2anal)) && ~any(trials2anal(~cftttrials))
            SVMgabor.analtrials = find(trials2anal);
        else
            error('trials to analyze was not selected correctly')
        end
    end
    SVMgabor.analtriallabels = trialorder(trials2anal);

    Ntesttrialspertype = floor(Ntrialspertype/Nsplits);
    Ntraintrialspertype = Ntrialspertype-Ntesttrialspertype;
    SVMgabor.Ntt = Ntt;
    SVMgabor.Ntrialspertype = Ntrialspertype;
    SVMgabor.Ntraintrialspertype = Ntraintrialspertype;

    Ntraintrials = Ntt*Ntraintrialspertype;
    Ntesttrials = Ntt*(Ntrialspertype-Ntraintrialspertype);

    % probe trials
    % cfprobetrials = ismember(rectrialorder, probetrialtypes);
    % SVMtrainICRC.cfprobetrials = find(cfprobetrials);

    SVMgabor.trialorder = trialorder;
    numtrials = length(trialorder);
    randtrialorder=randperm(numtrials);
    SVMgabor.randtrialorder = randtrialorder;

    SVMgabor_models = cell(1, Nsplits);
    SVMgabor.traintrialinds = zeros(Ntraintrials, Nsplits);
    SVMgabor.testtrialinds = zeros(Ntesttrials, Nsplits);
    % SVMout.Ylabs = cell(Ntt, Nsplits);
    for ts = 1:3
        switch ts
            case 1
                svmmd = 'train';
                tempNtrials = Ntraintrials;
            case 2
                svmmd = 'test';
                tempNtrials = Ntesttrials;
            case 3
                % svmmd = 'probe';
                % tempNtrials = nnz(cfprobetrials);
                svmmd = 'all';
                tempNtrials = numtrials;
        end
        SVMgabor.(svmmd).label = NaN(tempNtrials, Nsplits);
        SVMgabor.(svmmd).score = NaN(tempNtrials, Ntt, Nsplits);
    end

    % takes 20 min per trial type pair. 2000 min per session (33 hr)
    for isplit = 1:Nsplits
        close all
        ttclk = tic;

        % divide into train and test set
        % trials2anal = randtrialorder(ismember(randtrialorder, SVMtestRE.analtrials));
        testtrialinds = zeros(Ntesttrials,1);
        traintrialinds = zeros(Ntraintrials,1);
        for typi1 = 1:Ntt
            trialsintype = find( trialorder==traintrialtypes(typi1));
            tempinds = randtrialorder( trialorder(randtrialorder)==traintrialtypes(typi1) );
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
        testtrialinds = randtrialorder(ismember(randtrialorder, testtrialinds));
        traintrialinds = randtrialorder(ismember(randtrialorder, traintrialinds));

        SVMgabor.traintrialinds(:,isplit) = traintrialinds;
        SVMgabor.testtrialinds(:, isplit) = testtrialinds;

        switch preproc
            case 'none'
                Tp = tempR';
            case 'zscore'
                % Z-score
                trainRmean = mean(tempR(:,traintrialinds),2);
                trainRstd = std(tempR(:,traintrialinds),0,2);

                Tp = ( (tempR-trainRmean)./trainRstd )';
                Tp(:,trainRstd==0)=0;
            case 'minmax'
                trainRmin = min(tempR(:,traintrialinds),[],2);
                trainRrange = range(tempR(:,traintrialinds),2);

                Tp = ( (tempR-trainRmin)./trainRrange )';
        end

        X = Tp(traintrialinds,:);
        % Y = strsplit(sprintf('%d\n',rectrialorder(traintrialinds)), '\n')';
        % Y = Y(1:end-1);
        Y = reshape( trialorder(traintrialinds), [],1);

        %                 X = X(randomizedtraintrialorder, :);
        %                 Y = Y(randomizedtraintrialorder);

        Xtest = Tp(testtrialinds,:);
        % Ytest = strsplit(sprintf('%d\n',rectrialorder(testtrialinds)), '\n')';
        % Ytest = Ytest(1:end-1);
        Ytest = reshape( trialorder(testtrialinds), [],1);

        % Xprobe = Tp(cfprobetrials,:);
        Xall = Tp;

        % t is an SVM template. Most of its properties are empty.
        % When the software trains the ECOC classifier, it sets the applicable properties to their default values.
        % Train the ECOC classifier using the SVM template.
        % Transform classification scores to class posterior probabilities
        % (which are returned by predict or resubPredict) using the 'FitPosterior' name-value pair argument.
        % Specify the class order using the 'ClassNames' name-value pair argument.
        % Display diagnostic messages during training by using the 'Verbose' name-value pair argument.

        Ylabs = unique(Y);
        % Ylabs = Ylabs(randperm(numel(Ylabs)));
        % SVMtrainICRC.Ylabs(:,isplit) = Ylabs;
        SVMgabor.Ylabs = Ylabs;

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

        SVMgabor_models{isplit} = SVMModel;

        for t = 1:3
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
                    Xtemp = Xall;
                    tempSVMmodel = SVMModel;
                    svmmd = 'all';
            end
            [templabel,tempscore] = predict(tempSVMmodel,Xtemp);
            SVMgabor.(svmmd).label(:,isplit) = templabel;
            SVMgabor.(svmmd).score(:,:,isplit) = tempscore;
        end

        fprintf('%d-pixels %d/%d\n', rfpixs(irf), isplit, Nsplits)
        toc(ttclk)
    end

    save(svmfn, 'SVMgabor_models', 'SVMgabor', '-v7.3')
end

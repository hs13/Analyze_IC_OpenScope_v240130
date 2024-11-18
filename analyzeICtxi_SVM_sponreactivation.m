%% run inverted cross-validation SVM (non-overlapping training data)
addpath('C:\Users\USER\GitHub\helperfunctions')

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
    clearvars -except ises nwbsessions
    sesclk = tic;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    load([pathpp, 'postprocessed.mat'])

    % if 0, keep all neurons; if 1, exclude zero variance neurons in train trial
    % types; if 2 exclude zero variance neurons in all trial types
    twin = 0.4; % in sec. USING SPIKE COUNT RATHER THAN RATE FOR DECODING
    testt = [0 106 107 110 111];
    inferencett = [1105 1109];

    Ntt = numel(testt);
    whichblock = 'ICwcfg1_presentations';
    trialorder = vis.(whichblock).ICtrialtypes(vis.(whichblock).trialorder+1);
    tempspkcnt = 0.4*Rall.(whichblock)';

    preproc = 'zscore';
    whichSVMkernel = 'Linear';
    svmdesc = 'trainBK';
    Nsplits = 10;

    rng(111)
    % balance trials
    cftttrials = ismember(trialorder, testt);
    [v,numtrials]=uniquecnt(trialorder(cftttrials));
    Ntrialspertype = min(numtrials);
    if all(numtrials==Ntrialspertype)
        trials2anal = cftttrials;
    else
        warning('balancing number of trials')
        trials2anal = false(numrectrials,1);
        for typi1 = 1:Ntt
            trialsintype = find(trialorder==testt(typi1));
            trialsintype = trialsintype(1:Ntrialspertype);
            trials2anal(trialsintype) = true;
        end
    end

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

    traintrialinds = zeros(Ntraintrials, Nsplits);
    testtrialinds = zeros(Ntesttrials, Nsplits);
    for isplit = 1:Nsplits
        idxTrain = training(C,isplit);
        traintrialinds(:,isplit) = reshape( trials2analind(idxTrain) ,[],1);
        idxTest = test(C,isplit);
        testtrialinds(:,isplit) = reshape( trials2analind(idxTest) ,[],1);
    end


    % optimizeSVM: 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
    cvtrials = struct();
    cvtrials.optimizeSVM = 2;
    cvtrials.loadcvpartition = true;
    cvtrials.traintrialinds = testtrialinds; % NOTE THE SWITCH!
    cvtrials.testtrialinds = traintrialinds; % NOTE THE SWITCH!
    disp('inverting cross-validation: use non-overlapping training datasets')
    % SVMtrainICRC_models is 14 MB
    [SVMtrainBK, SVMtrainBK_models] = computeICtxi_SVM(tempspkcnt, trialorder, ...
        svmdesc, 'spkcnt', preproc, whichSVMkernel, cvtrials);

    % optimizeSVM: 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
    switch SVMtrainBK.optimizeSVM
        case 0
            optoptim = '_nooptim';
        case 1
            optoptim = '_alloptim';
        case 2
            optoptim = '';
        case 3
            optoptim = '_onevsall';
        otherwise
            error('optimizeSVM option %d not recognized', SVMtrainBK.optimizeSVM)
    end

    % train accuracy
    trainlabs = SVMtrainBK.trialorder(SVMtrainBK.spkcnt.traintrialinds);
    trainpred = SVMtrainBK.spkcnt.train.label;
    trainacc = zeros(Ntt);
    for itt = 1:Ntt
        trialsoi = trainlabs==testt(itt);
        [v,c]=uniquecnt(trainpred(trialsoi));
        trainacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
    end

    % test accuracy
    testlabs = SVMtrainBK.trialorder(SVMtrainBK.spkcnt.testtrialinds);
    testpred = SVMtrainBK.spkcnt.test.label;
    testacc = zeros(Ntt);
    for itt = 1:Ntt
        trialsoi = testlabs==testt(itt);
        [v,c]=uniquecnt(testpred(trialsoi));
        testacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
    end

    % inference decoding
    inftrials = ismember(SVMtrainBK.trialorder, inferencett);
    infpred = SVMtrainBK.spkcnt.all.label(inftrials,:);
    infperf = zeros(numel(inferencett), Ntt);
    for itt = 1:numel(inferencett)
        trialsoi = SVMtrainBK.trialorder(inftrials)==inferencett(itt);
        [v,c]=uniquecnt(infpred(trialsoi,:));
        infperf(itt, ismember(testt,v)) = c/(size(infpred,2)*nnz(trialsoi));
    end

    disp('SVM trainacc')
    disp(mean(trainacc,3))
    disp('SVM testacc')
    disp(mean(testacc,3))
    disp('SVM infperf')
    disp(mean(infperf,3))

    svmfn = strcat(pathpp, 'SVM_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
    svmmdlfn = strcat(pathpp, 'SVMmodels_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
    save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainBK', '-v7.3')
    save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainBK_models', '-v7.3')

    fprintf('%d/%d %s done running SVM %s\n', ises, numel(nwbsessions), nwbsessions{ises}, preproc)
    toc(sesclk)
end


%%
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};
Nsessions = numel(nwbsessions);
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_V1.mat'])
load([drivepath 'RESEARCH/ICexpts_revision23/openscope_psthavgall.mat'])

optimizeSVM = 2;

for ises = 1:numel(nwbsessions)
    clearvars -except optimizeSVM ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    preproc = 'meancenter';
    whichSVMkernel = 'Linear';
    svmdesc = 'trainICRC';
    % optimizeSVM: 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
    switch optimizeSVM
        case 0
            optoptim = '_nooptim';
        case 1
            optoptim = '_alloptim';
        case 2
            optoptim = '';
        case 3
            optoptim = '_onevsall';
        otherwise
            error('optimizeSVM option %d not recognized', optimizeSVM)
    end

    svmfn = strcat(pathpp, 'SVM_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
    svmmdlfn = strcat(pathpp, 'SVMmodels_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');

    load(svmfn)
    load(svmmdlfn)

    %% decode spontaneous activity
    load([pathpp, 'psth_spontaneous.mat'])
    load([pathpp, 'qc_units.mat'])

    neuV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
    neuvis = contains(neuallloc, 'VIS');
    neuRS = unit_wfdur>0.4;
    neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);

    neuV1RS = neuV1 & neuRS;
    %neuvisRS = neuvis & neuRS;

    spondurs = vis.spontaneous_presentations.stop_time-vis.spontaneous_presentations.start_time;
    [mv,mi]=max(spondurs);
    whichsponind = find(vis.spontaneous_presentations.start_time==vis.ICwcfg1_presentations.stop_time(end));
    fprintf('spontaneous bock after ICwcfg1 is %.3fs long\n', spondurs(whichsponind))
    if whichsponind~=mi
        warning('spontaneous bock after ICwcfg1 is not the longest spontaneous block')
    end

    tempspkcnt = cat(3,spkcntIChiV1agg{ises}{:});
    Nrep = size(tempspkcnt,1);
    Nneu = size(tempspkcnt,2);
    Nhireptt = numel(hireptt);
    spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu
    tempR = reshape(spkcntICtt, Nrep*Nhireptt, Nneu)';


    % preprocess spontaneous data
    % tempR is 400ms-window spike count formatted #neurons * #trials
    temppsth = psthspon{whichsponind}(:,neuV1RS);
    Ttot = size(temppsth,1);
    Twin = 400; % Tres must be 0.001s (1ms)
    Tslide = 25;
    Ntwins = floor( (Ttot-Twin)/Tslide );
    Tstartind = mod(Ttot-Twin, Tslide);
    Tspkinds = Tstartind+( (1:Twin)'+(0:Tslide:Ttot-Twin) );
    if Tspkinds(end,end) ~= Ttot
        error('check code')
    end
    Tctr = Tspkinds(round(Twin/2),:);
    tempspon = zeros( nnz(neuV1RS), size(Tspkinds,2) );
    for ci = 1:nnz(neuV1RS)
        temppsthvec = temppsth(:,ci);
        tempspon(ci,:) = sum(temppsthvec(Tspkinds), 1);
    end

    testt = SVMtrainICRC.trialtypes;
    Ntt = numel(testt);
    Nsplits = size(SVMtrainICRC.spkcnt.testtrialinds,2);
    SVMtrainICRC.spkcnt.spon.activity = tempspon';
    SVMtrainICRC.spkcnt.spon.label = NaN(size(Tspkinds,2), Nsplits);
    SVMtrainICRC.spkcnt.spon.score = NaN(size(Tspkinds,2), Ntt, Nsplits);
    for isplit = 1:Nsplits
        traintrialinds = SVMtrainICRC.spkcnt.traintrialinds(:,isplit);
        % testtrialinds = SVMtrainICRC.spkcnt.testtrialinds(:,isplit);
        switch preproc
            case 'none'
                Tp = tempspon';
            case 'zscore'
                % Z-score
                trainRmean = mean(tempR(:,traintrialinds),2);
                trainRstd = std(tempR(:,traintrialinds),0,2);

                Tp = ( (tempspon-trainRmean)./trainRstd )';
            case 'minmax'
                trainRmin = min(tempR(:,traintrialinds),[],2);
                trainRrange = range(tempR(:,traintrialinds),2);

                Tp = ( (tempspon-trainRmin)./trainRrange )';
            case 'meancenter'
                trainRmean = mean(tempR(:,traintrialinds),2);
                Tp = (tempspon-trainRmean)';
        end
        Tp(isnan(Tp))=0; % Ntrials * Nneurons
        [templabel,tempscore] = predict(SVMtrainICRC_models.spkcnt{isplit}, Tp);

        SVMtrainICRC.spkcnt.spon.label(:,isplit) = templabel;
        SVMtrainICRC.spkcnt.spon.score(:,:,isplit) = tempscore;
    end

    % calculate reactivation coefficient
    reactivcoeff = zeros(nnz(neuV1RS), Ntt);
    for itt = 1:Ntt
        sponpredtt = mean(SVMtrainICRC.spkcnt.spon.label==testt(itt),2);
        reactivcoeff(:,itt) = corr(SVMtrainICRC.spkcnt.spon.activity, sponpredtt);
    end


    svmsponfn = strcat(pathpp, 'SVMspon', num2str(whichsponind), '_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
    save(svmsponfn, 'neuV1RS', 'Twin', 'Tslide', 'Tstartind', 'Tctr', 'SVMtrainICRC', 'reactivcoeff')

    %
    pltses = false;
    if pltses
        load([pathpp, 'visresponses.mat'])
        figure;
        for itt = 1:4
            switch itt
                case 1
                    neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RS);
                    neuoj = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RS);
                case 2
                    neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RS);
                    neuoj = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RS);
                case 3
                    neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RS);
                    neuoj = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RS);
                case 4
                    neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RS);
                    neuoj = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RS);
                otherwise
                    error('check itt')
            end
            subplot(2,2,itt)
            hold all
            h=histogram(reactivcoeff(:,itt));
            histogram(reactivcoeff(neuoi==1,itt), 'BinEdges', h.BinEdges)
            histogram(reactivcoeff(neuoj==1,itt), 'BinEdges', h.BinEdges)
            title(testt(itt))
        end

        kerwinhalf = 12; kersigma = 5;
        kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
        kergauss = (kergauss/sum(kergauss));
        sespsthavg = psthavgall.ICwcfg1_presentations(:,:,sesneuall==ises);
        temppsthavg = convn(sespsthavg(:,:,neuV1RS), kergauss, 'same');
        neu2plt = find(neuV1RS);

        tt2p = testt;
        ttcol = [0 .4 0; .5 0.25 0; 1 0.5 0; 0 1 0];
        Ntop = 7;
        figure
        for itt = 1:4
            [sv,si]=sort(reactivcoeff(:,itt), 'descend', 'MissingPlacement', 'last');
            for ii = 1:Ntop
                subplot(Ntt,Ntop,Ntop*(itt-1)+ii)
                hold all
                for typi = 1:numel(tt2p)
                    plot(psthtli, temppsthavg(:,ICtrialtypes==tt2p(typi),si(ii)), 'Color', ttcol(typi,:), 'LineWidth', 1)
                end
                xlim([-200 600])
                title(sprintf('Trial%d Cell%d reactivation coeff. %.4f', testt(itt), neu2plt(si(ii)), sv(ii)))
            end
        end
    end


end
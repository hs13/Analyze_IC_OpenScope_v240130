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

% for ises = 1:numel(nwbsessions)
pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-619293/'];
% pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-620333/'];

load([pathpp 'SVM_trainICRC_V1_Linear_meancenter.mat'])
load([pathpp 'SVMmodels_trainICRC_V1_Linear_meancenter.mat'])
load([pathpp 'SVM_trainICRC_V1_Linear_meancenter_lmlvslopes.mat'])
load([pathpp 'SVMmodels_trainICRC_V1_Linear_meancenter_lmlvslopes.mat'])

testt = [106,107,110,111];
inferencett = [1105 1109];
Nhireptt = numel(hireptt);
Ntt = numel(testt);
Nsplits = size(SVMtrainICRC.spkcnt.testtrialinds,2);

traincmlmlvs = zeros(Ntt, Ntt, numel(disperses));
testcmlmlvs = zeros(Ntt, Ntt, numel(disperses));
infcmlmlvs = zeros(numel(inferencett), Ntt, numel(disperses));

trainperflmlvs = zeros(1, numel(disperses));
testperflmlvs = zeros(1, numel(disperses));
infperflmlvs = zeros(1, numel(disperses));

for islope = 0:numel(disperses)
    if islope==0
        SVMout = SVMtrainICRC;
    else
        SVMout = SVMtrainICRC_lmlvs(islope);
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

%% to test robustness silence a portion of neurons and see change in classification accuracy across slopes
% parametrically change proportion silenced
propneusilvec = 0:0.1:1;
% randomly select neurons to silence: sample 1000X
Nsamples = 1000;
Nneurons = SVMtrainICRC.Nneurons;
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
            neu2silmat(neurand,:) = true;
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
            
            excludeneuvar0 = 2;
            switch excludeneuvar0
                case 0
                    valneu = true(Nneu,1);
                case 1
                    valneu = squeeze(all(newspkcntvar(1,ismember(hireptt,testt),:)>0, 2));
                case 2
                    valneu = squeeze(all(newspkcntvar(1,:,:)>0, 2));
                otherwise
                    error('excludeneuvar0 option not recognized')
            end
            
            tempR = reshape(newspkcntICtt, Nrep*Nhireptt, nnz(valneu))';
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
        
        trainaccpd = zeros(Nsamples, Ntt,Ntt, Nsplits);
        testaccpd = zeros(Nsamples, Ntt,Ntt, Nsplits);
        infaccpd = zeros(Nsamples, numel(inferencett),Ntt, Nsplits);
        trainperfpd = zeros(Nsamples, Nsplits);
        testperfpd = zeros(Nsamples, Nsplits);
        infperfpd = zeros(Nsamples, Nsplits);
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
        else
            silrandtraincmlmlvs{islope}(:,n,:,:) = mean(trainaccpd,4);
            silrandtestcmlmlvs{islope}(:,n,:,:) = mean(testaccpd,4);
            silrandinfcmlmlvs{islope}(:,n,:,:) = mean(infaccpd,4);
            
            silrandtrainperflmlvs{islope}(:,n) = mean(trainperfpd,2);
            silrandtestperflmlvs{islope}(:,n) = mean(testperfpd,2);
            silrandinfperflmlvs{islope}(:,n) = mean(infperfpd,2);
        end
        
    end
end

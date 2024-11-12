if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])

pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-619293/'];
% pathpp= [drivepath 'RESEARCH/logmean_logvar/sub-620333/'];

load([pathpp 'SVM_trainICRC_V1_Linear_meancenter.mat'])
load([pathpp 'SVMmodels_trainICRC_V1_Linear_meancenter.mat'])
load([pathpp 'SVM_trainICRC_V1_Linear_meancenter_lmlvslopes.mat'])
load([pathpp 'SVMmodels_trainICRC_V1_Linear_meancenter_lmlvslopes.mat'])

testt = [106,107,110,111];
Ntt = numel(testt);
inferencett = [1105 1109];
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
    trainlabs = SVMtrainICRC.trialorder(SVMout.spkcnt.traintrialinds);
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

figure;hold all
plot(disperses, testperflmlvs)
xl = [disperses(1) disperses(end)];
plot(xl, testperfasis*[1 1], 'c-')

%% to test robustness silence a portion of neurons and see change in classification accuracy across slopes
% parametrically change proportion silenced 
propneusilvec = 0:0.1:1;
%for n = 1:length(propneusilvec)
n = 6;
propneu2sil = propneusilvec(n);
% randomly select neurons to silence: sample 1000X
Nsample = 1000;

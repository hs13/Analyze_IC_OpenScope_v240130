if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])


ises = 3;
nwbsessions{ises}

%% comparison point: linear SVM
load(['S:\OpenScopeData\00248_v240130\SVM_trainICRC_selectareas\' nwbsessions{ises} '\SVM_trainICRC_VISpRS_Linear_zscore_ICwcfg1.mat'])
Ntt = numel(SVMtrainICRC.trialtypes);

testtrialinds = SVMtrainICRC.spkcnt.testtrialinds;

% test accuracy
testlabs = SVMtrainICRC.trialorder(SVMtrainICRC.spkcnt.testtrialinds);
testpred = SVMtrainICRC.spkcnt.test.label;
testacc = zeros(Ntt);
for itt = 1:Ntt
    trialsoi = testlabs==SVMtrainICRC.trialtypes(itt);
    [v,c]=uniquecnt(testpred(trialsoi));
    testacc(itt, ismember(SVMtrainICRC.trialtypes,v)) = c/nnz(trialsoi);
end

% inference decoding
inferencett = [1105 1109];
inftrials = ismember(SVMtrainICRC.trialorder, inferencett);
infpred = SVMtrainICRC.spkcnt.all.label(inftrials,:);
infperf = zeros(numel(inferencett), Ntt);
for itt = 1:numel(inferencett)
    trialsoi = SVMtrainICRC.trialorder(inftrials)==inferencett(itt);
    [v,c]=uniquecnt(infpred(trialsoi));
    infperf(itt, ismember(SVMtrainICRC.trialtypes,v)) = c/nnz(trialsoi);
end

%% Bayesian image decoding (inspired by position decoding in hippocampus literature, e.g., Buzsaki lab, Fenton lab)



%% Naive Bayes Decoder (independent neurons)

%% AODE (Sugden...Andermann 2020)

%% UMAP + Bayesian decoding
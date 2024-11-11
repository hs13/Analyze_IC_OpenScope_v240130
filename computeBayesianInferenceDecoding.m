if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/Analyze_IC_OpenScope_v240130/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])

twin = 0.4; % in sec
testt = [106 107 110 111];
Ntt = numel(testt);
inferencett = [1105 1109];

ises = 3;
spkcntses = spkcntIChiV1agg{ises};

%% comparison point: linear SVM
load(['S:\OpenScopeData\00248_v240130\SVM_trainICRC_selectareas\' nwbsessions{ises} '\SVM_trainICRC_VISpRS_Linear_zscore_ICwcfg1.mat'])

if ~isequal(SVMtrainICRC.trialtypes, testt)
    error('mismatch in test trial types: check that you loaded trainICRC')
end

% test accuracy
testlabs = SVMtrainICRC.trialorder(SVMtrainICRC.spkcnt.testtrialinds);
testpred = SVMtrainICRC.spkcnt.test.label;
testacc = zeros(Ntt);
for itt = 1:Ntt
    trialsoi = testlabs==testt(itt);
    [v,c]=uniquecnt(testpred(trialsoi));
    testacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
end

% inference decoding
inftrials = ismember(SVMtrainICRC.trialorder, inferencett);
infpred = SVMtrainICRC.spkcnt.all.label(inftrials,:);
infperf = zeros(numel(inferencett), Ntt);
for itt = 1:numel(inferencett)
    trialsoi = SVMtrainICRC.trialorder(inftrials)==inferencett(itt);
    [v,c]=uniquecnt(infpred(trialsoi));
    infperf(itt, ismember(testt,v)) = c/nnz(trialsoi);
end

%% Bayesian image decoding (inspired by position decoding in hippocampus literature, e.g., Buzsaki lab, Fenton lab)

filtervalneu = false; % whether to filter out neurons that hardly ever fire
% NO NEED TO WORRY about this because 0^0=1

% train vs test trial divide
if exist('SVMtrainICRC', 'var')
    kfold = size(SVMtrainICRC.spkcnt.testtrialinds,2);
    traintrialinds = cell(size(hireptt));
    testtrials = cell(size(hireptt));
    for itt = 1:Ntt
        trialsoind = find(SVMtrainICRC.trialorder==inferencett(itt));
        typi = hireptt==testt(itt);
        ntrials = size(spkcntses{typi},1);
        testtrials{typi} = false(ntrials,kfold);
        for k = 1:kfold
            testtrials{typi}(:,k) = ismember(trialsoind, SVMtrainICRC.spkcnt.testtrialinds(:,k));
        end
    end
else
    kfold = 10;
    testtrials = cell(size(hireptt));
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        ntrials = size(spkcntses{typi},1);
        testtrials{typi} = false(ntrials,kfold);
        c = cvpartition(ntrials,"KFold",kfold);
        for k = 1:kfold
            testtrials{typi}(:,k) = test(c,k);
        end
    end
end

Nneurons = size(spkcntses{1},2);

% filter out neurons whose FRtrainavg could be zero
neuval = true(Nneurons,1);
if filtervalneu
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        tempvalneu = mean(spkcntses{typi}>0,1) >= 1/kfold;
        neuval(~tempvalneu) = false;
    end
end

bayesimage = struct();
bayesimage.postprob = cell(numel(hireptt),kfold);
bayesimage.postprobnorm = cell(numel(hireptt),kfold);
bayesimage.trainacc = zeros(Ntt,Ntt,kfold);
bayesimage.testacc = zeros(Ntt,Ntt,kfold);
bayesimage.infperf = zeros(numel(inferencett),Ntt,kfold);
for k = 1:kfold
    
    FRtrainavg = NaN(nnz(neuval), Ntt);
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        FRtrainavg(:,itt) = (1/twin) * mean(spkcntses{typi}(~testtrials{typi}(:,k), neuval),1);
    end
    
    for typi = 1:numel(hireptt)
        ntrials = size(spkcntses{typi},1);
        bayesimage.postprob{typi,k} = NaN( ntrials, Ntt ); % intialize
        for itrial = 1:ntrials
            spkcntvec = spkcntses{typi}(itrial,neuval)';
            for jtt = 1:Ntt
                FRimvec = FRtrainavg(:,jtt);
                bayesimage.postprob{typi,k}(itrial,jtt) = prod(FRimvec.^spkcntvec) * exp(-twin*sum(FRimvec));
            end
        end
        bayesimage.postprobnorm{typi,k} = bayesimage.postprob{typi,k}./sum(bayesimage.postprob{typi,k},2);
    end
        
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        temptesttrials = testtrials{typi}(:,k);
        [mv,mi] = max( bayesimage.postprobnorm{typi,k}(~temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        bayesimage.trainacc(itt, ismember([1,2,3,4], v), k) = c/nnz(~temptesttrials);
        
        [mv,mi] = max( bayesimage.postprobnorm{typi,k}(temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        bayesimage.testacc(itt, ismember([1,2,3,4], v), k) = c/nnz(temptesttrials);
    end    
    
    for itt = 1:numel(inferencett)
        typi = hireptt==inferencett(itt);
        ntrials = size(bayesimage.postprobnorm{typi,k},1);
        [mv,mi] = max( bayesimage.postprobnorm{typi,k},[],2);
        [v,c]=uniquecnt(mi);
        bayesimage.infperf(itt, ismember([1,2,3,4], v), k) = c/ntrials;
    end
end

disp('bayesimage.trainacc')
disp(mean(bayesimage.trainacc,3))
disp('bayesimage.testacc')
disp(mean(bayesimage.testacc,3))
disp('bayesimage.infperf')
disp(mean(bayesimage.infperf,3))

%% Naive Bayes Gaussian Decoder: fit gaussian to log spike count
% contruct P(r|s) assuming that spike counts are log normal (train vs test split)
% sum P(r|s) across neurons
% normalize across stimuli

% adding spike count offset of 1 *increases* train, test and inference
% performance!!!
spkcntoffset = 1;

%{
% check gaussian fit of log(spike counts)
typi = hireptt==111;
[sv,si] = sort(mean(spkcntses{typi},1), 'descend');
neuex = si(100);

figure
for offset = 0:1

    tempspk = offset+spkcntses{typi}(:,neuex);

xval = log10(tempspk(tempspk>0));
[muHat,sigmaHat] = normfit(xval);
mu = mean(xval);
sigma = std(xval);
subplot(2,2,1+offset*2); hold all
h = histogram(xval, 'normalization', 'pdf');
xt = h.BinEdges(1):0.001:h.BinEdges(end);
y = normpdf(xt,mu,sigma);
yHat = normpdf(xt,muHat,sigmaHat);
y1 = normpdf(xt,mean(log10(tempspk+10^-1)), std(log10(tempspk+10^-1)) );
y2 = normpdf(xt,mean(log10(tempspk+10^0)), std(log10(tempspk+10^0)) );
plot(xt, y, 'k-', 'linewidth',2)
plot(xt, yHat, 'r--', 'linewidth',1.5)
plot(xt, y1, 'b--', 'linewidth',2)
plot(xt, y2, 'c--', 'linewidth',2)

subplot(2,2,2+offset*2); hold all
h = histogram(tempspk, 'normalization', 'pdf');
xt = log10(h.BinEdges(2:end));
y = normpdf(xt,mu,sigma);
yHat = normpdf(xt,muHat,sigmaHat);
y1 = normpdf(xt,mean(log10(tempspk+10^-1)), std(log10(tempspk+10^-1)) );
y2 = normpdf(xt,mean(log10(tempspk+10^0 )), std(log10(tempspk+10^0)) );
plot(10.^xt, 10.^y*sum(h.Values)/sum(10.^y), 'k-', 'linewidth',2)
plot(10.^xt, 10.^yHat*sum(h.Values)/sum(10.^yHat), 'r--', 'linewidth',2)
plot(10.^xt, 10.^y1*sum(h.Values)/sum(10.^y1), 'b--', 'linewidth',2)
plot(10.^xt, 10.^y2*sum(h.Values)/sum(10.^y2), 'c--', 'linewidth',2)
end
%}

naivegauss = struct();
naivegauss.postprob = cell(numel(hireptt),kfold);
naivegauss.postprobnorm = cell(numel(hireptt),kfold);
naivegauss.trainacc = zeros(Ntt,Ntt,kfold);
naivegauss.testacc = zeros(Ntt,Ntt,kfold);
naivegauss.infperf = zeros(numel(inferencett),Ntt,kfold);
for k = 1:kfold
    for typi = 1:numel(hireptt)
        ntrials = size(spkcntses{typi},1);
        ratelikelihood = zeros(Nneurons, ntrials, Ntt);
        for jtt = 1:Ntt
            typj = hireptt==testt(jtt);
            temptesttrials = testtrials{typj}(:,k);
            for ci = 1:Nneurons
                tempspk = spkcntoffset+spkcntses{typj}(~temptesttrials, ci);
                xval = log10(tempspk(tempspk>0));
                mu = mean(xval);
                sigma = std(xval);
                xtt = log10(spkcntoffset+spkcntses{typi}(:,ci));
                ratelikelihood(ci,:,jtt) = normpdf(xtt,mu,sigma);
            end
        end
        naivegauss.postprob{typi,k} = squeeze(nansum(ratelikelihood,1));
        
%         % sum across neurons then normalize across trialtypes
%         normsumratelikelihood = squeeze( nansum(ratelikelihood,1)./sum(nansum(ratelikelihood,1),3) );
%         naivegauss.postprobnorm{typi,k} = normsumratelikelihood;
        
        % normalize each neuron, then sum across neurons
        % this leads to higher train and test accuracy
        sumnormratelikelihood = squeeze(nansum( ratelikelihood./sum(ratelikelihood,3), 1));
        naivegauss.postprobnorm{typi,k} = sumnormratelikelihood;
    end

    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        temptesttrials = testtrials{typi}(:,k);
        [mv,mi] = max( naivegauss.postprobnorm{typi,k}(~temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        naivegauss.trainacc(itt, ismember([1,2,3,4], v), k) = c/nnz(~temptesttrials);
        
        [mv,mi] = max( naivegauss.postprobnorm{typi,k}(temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        naivegauss.testacc(itt, ismember([1,2,3,4], v), k) = c/nnz(temptesttrials);
    end    
    
    for itt = 1:numel(inferencett)
        typi = hireptt==inferencett(itt);
        ntrials = size(naivegauss.postprobnorm{typi,k},1);
        [mv,mi] = max( naivegauss.postprobnorm{typi,k},[],2);
        [v,c]=uniquecnt(mi);
        naivegauss.infperf(itt, ismember([1,2,3,4], v), k) = c/ntrials;
    end
end

disp('naivegauss.trainacc')
disp(mean(naivegauss.trainacc,3))
disp('naivegauss.testacc')
disp(mean(naivegauss.testacc,3))
disp('naivegauss.infperf')
disp(mean(naivegauss.infperf,3))

%% Bayesian multivariate gaussian decoder
% for N neurons, fit N-dimensional multivarite gaussian
% unlike Naive Bayes, covariance is factored in
% however, fitting multivariate gaussian performs much more poorly than naive bayes, probably due to poor fit

filtervalneu = false; % whether to filter out neurons that hardly ever fire
% filter out neurons whose FRtrainavg could be zero
neuval = true(Nneurons,1);
if filtervalneu
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        tempvalneu = mean(spkcntses{typi}>0,1) >= 1/kfold;
        neuval(~tempvalneu) = false;
    end
end

%{
% check multivariate gaussian fit of log(spike counts) for a pair of neurons
typi = hireptt==106;
[sv,si] = sort(mean(spkcntses{typi},1), 'descend');
neuexi = si(1);
neuexj = si(2);

tempspki = spkcntoffset+spkcntses{typi}(:,neuexi);
tempspkj = spkcntoffset+spkcntses{typi}(:,neuexj);

tempspkn0 = [tempspki tempspkj];
tempspkn0(tempspkn0==0) = NaN;
X = log10(tempspkn0);
mu = nanmean(X,1);
Sigma = cov(X, 'partialrows');

% xxt = log10(1):0.001:log10(max(tempspki));
% yyt = log10(1):0.001:log10(max(tempspkj));
xxt = log10(1:1:max(tempspki));
yyt = log10(1:1:max(tempspkj));
gridxyticks = zeros(length(yyt), length(xxt), 2);
gridxyticks(:,:,1) = repmat(xxt, length(yyt), 1);
gridxyticks(:,:,2) = repmat(yyt', 1, length(xxt) );
gridxy = reshape(gridxyticks, length(yyt)*length(xxt), 2);
gridpdf = mvnpdf(gridxy,mu,Sigma);
gridpdf = reshape(gridpdf, length(yyt), length(xxt));

redcm = 1-gray;
redcm(:,1)=1;
figure
hold all
imagesc(10.^xxt, 10.^yyt, gridpdf);%, 'AlphaData',0.1)
scatter(tempspki, tempspkj, 'o', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k')
colormap(redcm)
%}

mvnbayes = struct();
mvnbayes.postprob = cell(numel(hireptt),kfold);
mvnbayes.postprobnorm = cell(numel(hireptt),kfold);
mvnbayes.trainacc = zeros(Ntt,Ntt,kfold);
mvnbayes.testacc = zeros(Ntt,Ntt,kfold);
mvnbayes.infperf = zeros(numel(inferencett),Ntt,kfold);
for k = 1:kfold
    for typi = 1:numel(hireptt)
        ntrials = size(spkcntses{typi},1);
        mvnbayes.postprob{typi,k} = NaN( ntrials, Ntt );
        for jtt = 1:Ntt
            typj = hireptt==testt(jtt);
            temptesttrials = testtrials{typj}(:,k);
            
            tempspkn0 = spkcntoffset+spkcntses{typi}(~temptesttrials,neuval);
            tempspkn0(tempspkn0==0) = NaN;
            X = log10(tempspkn0);
            mu = nanmean(X,1);
            Sigma = cov(X, 'partialrows');
            xtt = log10( spkcntoffset+spkcntses{typi}(:,neuval) );
            try
                mvnbayes.postprob{typi,k}(:,jtt) = mvnpdf(xtt,mu,Sigma);
            catch
                % SIGMA must be a square, symmetric, positive definite matrix.
                % Positive definite matrix: all eigenvalues are positive;
                % the matrix is guaranteed to be invertible.
                % Step 1: Symmetrize A
                Sigma = (Sigma + Sigma') / 2;
                % Step 2: Adjust eigenvalues to ensure positive definiteness
                [eigVec, eigVal] = eig(Sigma);
                eigVal(eye(size(eigVal)) & eigVal<=1e-20) = 1e-20; % Shift non-positive eigenvalues
                Sigma = eigVec * eigVal * eigVec';
                mv = max(abs(cov(X, 'partialrows')-eigVec * eigVal * eigVec'),[],'all');
                if mv>10^-10
                warning('difference after making Sigma square, symmetric, positive definite')
                end
                %figure; plot(cov(X, 'partialrows'), eigVec * eigVal * eigVec', '.')
                
                mvnbayes.postprob{typi,k}(:,jtt) = mvnpdf(xtt,mu,Sigma);
            end
            
        end
        mvnbayes.postprobnorm{typi,k} = mvnbayes.postprob{typi,k}./sum(mvnbayes.postprob{typi,k},2);
    end
    
    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        temptesttrials = testtrials{typi}(:,k);
        [mv,mi] = max( mvnbayes.postprobnorm{typi,k}(~temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        mvnbayes.trainacc(itt, ismember([1,2,3,4], v), k) = c/nnz(~temptesttrials);
        
        [mv,mi] = max( mvnbayes.postprobnorm{typi,k}(temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        mvnbayes.testacc(itt, ismember([1,2,3,4], v), k) = c/nnz(temptesttrials);
    end
    
    for itt = 1:numel(inferencett)
        typi = hireptt==inferencett(itt);
        ntrials = size(mvnbayes.postprobnorm{typi,k},1);
        [mv,mi] = max( mvnbayes.postprobnorm{typi,k},[],2);
        [v,c]=uniquecnt(mi);
        mvnbayes.infperf(itt, ismember([1,2,3,4], v), k) = c/ntrials;
    end
end

disp('mvnbayes.trainacc')
disp(mean(mvnbayes.trainacc,3))
disp('mvnbayes.testacc')
disp(mean(mvnbayes.testacc,3))
disp('mvnbayes.infperf')
disp(mean(mvnbayes.infperf,3))

%% Naive Bayes Gaussian PCA Decoder: fit gaussian to PCA
% independence assumtion becomes valid after PCA
logspkopt = false; % logspkopt==true returns a much poorer performance!!!

pcanaivegauss = struct();
pcanaivegauss.postprob = cell(numel(hireptt),kfold);
pcanaivegauss.postprobnorm = cell(numel(hireptt),kfold);
pcanaivegauss.trainacc = zeros(Ntt,Ntt,kfold);
pcanaivegauss.testacc = zeros(Ntt,Ntt,kfold);
pcanaivegauss.infperf = zeros(numel(inferencett),Ntt,kfold);
for k = 1:kfold
    for typi = 1:numel(hireptt)
        ntrials = size(spkcntses{typi},1);
        

        ratelikelihood = zeros(Nneurons, ntrials, Ntt);
        for jtt = 1:Ntt
            typj = hireptt==testt(jtt);
            temptesttrials = testtrials{typj}(:,k);
            if logspkopt
                tempspk = spkcntoffset+spkcntses{typj}(~temptesttrials, :);
                Xval = log10(tempspk(tempspk>0));
                Xtt = log10(spkcntoffset+spkcntses{typi});
            else
                Xval = spkcntses{typj}(~temptesttrials, :);
                Xtt = spkcntses{typi};
            end
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(Xval);
            mcspkcnt = Xtt-MU;
            SCOREtt = mcspkcnt*COEFF;
            for ci = 1:size(COEFF,2)
                mu = mean(SCORE(:,ci));
                sigma = std(SCORE(:,ci));
                ratelikelihood(ci,:,jtt) = normpdf(SCOREtt(:,ci),mu,sigma);
                
                %{
                % the following takes much longer and does not improve performance
                try
                    GMModel = fitgmdist(SCORE(:,ci),3);
                    ratelikelihood(ci,:,jtt) = pdf(GMModel, SCOREtt(:,ci));
                catch
                    GMModel = fitgmdist(SCORE(:,ci),1);
                    ratelikelihood(ci,:,jtt) = pdf(GMModel, SCOREtt(:,ci));
                end
                %}
            end
        end
        pcanaivegauss.postprob{typi,k} = squeeze(nansum(ratelikelihood,1));
        
%         % sum across neurons then normalize across trialtypes
%         normsumratelikelihood = squeeze( nansum(ratelikelihood,1)./sum(nansum(ratelikelihood,1),3) );
%         pcanaivegauss.postprobnorm{typi,k} = normsumratelikelihood;
        
        % normalize each neuron, then sum across neurons
        % this leads to higher train and test accuracy
        sumnormratelikelihood = squeeze(nansum( ratelikelihood./sum(ratelikelihood,3), 1));
        pcanaivegauss.postprobnorm{typi,k} = sumnormratelikelihood;
    end

    for itt = 1:Ntt
        typi = hireptt==testt(itt);
        temptesttrials = testtrials{typi}(:,k);
        [mv,mi] = max( pcanaivegauss.postprobnorm{typi,k}(~temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        pcanaivegauss.trainacc(itt, ismember([1,2,3,4], v), k) = c/nnz(~temptesttrials);
        
        [mv,mi] = max( pcanaivegauss.postprobnorm{typi,k}(temptesttrials,:),[],2);
        [v,c]=uniquecnt(mi);
        pcanaivegauss.testacc(itt, ismember([1,2,3,4], v), k) = c/nnz(temptesttrials);
    end    
    
    for itt = 1:numel(inferencett)
        typi = hireptt==inferencett(itt);
        ntrials = size(pcanaivegauss.postprobnorm{typi,k},1);
        [mv,mi] = max( pcanaivegauss.postprobnorm{typi,k},[],2);
        [v,c]=uniquecnt(mi);
        pcanaivegauss.infperf(itt, ismember([1,2,3,4], v), k) = c/ntrials;
    end
end

disp('pcanaivegauss.trainacc')
disp(mean(pcanaivegauss.trainacc,3))
disp('pcanaivegauss.testacc')
disp(mean(pcanaivegauss.testacc,3))
disp('pcanaivegauss.infperf')
disp(mean(pcanaivegauss.infperf,3))
mean(diag(mean(pcanaivegauss.testacc,3)))

%% Naive Bayes PCA Decoder: 


%% Naive Bayes Decoder (independent neurons)
% contruct P(r|s) based on observed spike counts (train vs test split)

%% AODE (Sugden...Andermann 2020)
% acrivity was binarized based on a threshold


%% UMAP + Bayesian decoding
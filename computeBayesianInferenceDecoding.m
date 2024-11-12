if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/Analyze_IC_OpenScope_v240130/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_V1.mat'])

warning('off')
for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    sesclk = tic;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    computeSVM = false;
    if computeSVM
        preproc = 'meancenter';
        whichSVMkernel = 'Linear';
        svmdesc = 'trainICRC';
    end
    testt = [106 107 110 111];
    inferencett = [1105 1109];

    twin = 0.4; % in sec
    disperses = -1:0.1:2;
    Ntt = numel(testt);
    Nhireptt = numel(hireptt);

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

    if computeSVM
        SVMtrainICRC_lmlvs = struct();
        SVMtrainICRC_models_lmlvs = struct();
    end
    validneurons_lmlvs = true(Nneu,length(disperses));
    newspkcnt_lmlvs = cell(Nhireptt,length(disperses));
    bayesimage_lmlvs = struct();
    naivegauss_lmlvs = struct();
    mvnbayes_lmlvs = struct();
    pcanaivegauss_lmlvs = struct();
    for islope = 0:numel(disperses)
        if islope==0
            spkcntlmlvs = spkcntIChiV1agg{ises};
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

            valneu = squeeze(all(newspkcntvar(1,:,:)>0, 2));
            validneurons_lmlvs(:,islope) = valneu;

            spkcntlmlvs = cell(Nhireptt,1);
            for n = 1:Nhireptt
                spkcntlmlvs{n} = squeeze(newspkcntICtt(:,n,valneu));
            end
            newspkcnt_lmlvs(:,islope) = spkcntlmlvs;

            tempR = reshape(newspkcntICtt(:,:,valneu), Nrep*Nhireptt, nnz(valneu))';
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);
        end
        % tempinvalneu = squeeze(any(newspkcntvar(1,:,:)<=0, 2));
        % typi = 4;
        % figure; hold all
        % plot(log(mean(spkcntIChiV1agg{ises}{typi},1)), log(var(spkcntIChiV1agg{ises}{typi},0,1)), 'k.')
        % plot(log(mean(spkcntlmlvs{typi},1)), log(var(spkcntlmlvs{typi},0,1)), 'r.')
        % plot(log(mean(spkcntlmlvs{typi}(:,tempinvalneu),1)), log(var(spkcntlmlvs{typi}(:,tempinvalneu),0,1)), 'ro')

        %% comparison point: linear SVM ~7min per slope
        % load(['S:\OpenScopeData\00248_v240130\SVM_trainICRC_selectareas\' nwbsessions{ises} '\SVM_trainICRC_VISpRS_Linear_zscore_ICwcfg1.mat'])
        % SVMtrainICRC_models is 14 MB
        if computeSVM % && (islope==0 || mod(disperses(islope)*10,5)==0)
            kfold = size(SVMtrainICRC.spkcnt.testtrialinds,2);
            traintrials = false(length(trialorder),kfold);
            testtrials = false(length(trialorder),kfold);
            for itt = 1:Ntt
                trialsoind = find(SVMtrainICRC.trialorder==testt(itt));
                typi = hireptt==testt(itt);
                ntrials = size(spkcntlmlvs{typi},1);
                temptesttrials = false(ntrials,kfold);
                for k = 1:kfold
                    temptesttrials(:,k) = ismember(trialsoind, SVMtrainICRC.spkcnt.testtrialinds(:,k));
                end
                traintrials(trialorder==testt(itt),:) = ~temptesttrials;
                testtrials(trialorder==testt(itt),:) = temptesttrials;
            end
            traintrialinds = zeros(max(sum(traintrials,1)), kfold);
            testtrialinds = zeros(max(sum(testtrials,1)), kfold);
            for k = 1:kfold
                traintrialinds(:,k) = find(traintrials(:,k));
                testtrialinds(:,k) = find(testtrials(:,k));
            end

            cvtrials = struct();
            cvtrials.loadcvpartition = true;
            cvtrials.traintrialinds = traintrialinds;
            cvtrials.testtrialinds = testtrialinds;

            [SVMtrainICRC, SVMtrainICRC_models] = computeICtxi_SVM(tempR, trialorder, ...
                svmdesc, 'spkcnt', preproc, whichSVMkernel, cvtrials);

            if ~isequal(SVMtrainICRC.trialtypes, testt)
                error('mismatch in test trial types: check that you loaded trainICRC')
            end

            % train accuracy
            trainlabs = SVMtrainICRC.trialorder(SVMtrainICRC.spkcnt.traintrialinds);
            trainpred = SVMtrainICRC.spkcnt.train.label;
            trainacc = zeros(Ntt);
            for itt = 1:Ntt
                trialsoi = trainlabs==testt(itt);
                [v,c]=uniquecnt(trainpred(trialsoi));
                trainacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
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

            disp('SVM trainacc')
            disp(mean(trainacc,3))
            disp('SVM testacc')
            disp(mean(testacc,3))
            disp('SVM infperf')
            disp(mean(infperf,3))
        else
            SVMtrainICRC = struct();
            SVMtrainICRC_models = struct();
        end

        %% Bayesian image decoding (inspired by position decoding in hippocampus literature, e.g., Buzsaki lab, Fenton lab)

        filtervalneu = false; % whether to filter out neurons that hardly ever fire
        % NO NEED TO WORRY about this because 0^0=1

        % train vs test trial divide
        load(['S:\OpenScopeData\00248_v240130\SVM_trainICRC_selectareas\' nwbsessions{ises} '\SVM_trainICRC_VISpRS_Linear_zscore_ICwcfg1.mat'])
        if exist('SVMtrainICRC', 'var')
            kfold = size(SVMtrainICRC.spkcnt.testtrialinds,2);
            testtrials = cell(size(hireptt));
            for itt = 1:Ntt
                trialsoind = find(SVMtrainICRC.trialorder==testt(itt));
                typi = hireptt==testt(itt);
                ntrials = size(spkcntlmlvs{typi},1);
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
                ntrials = size(spkcntlmlvs{typi},1);
                testtrials{typi} = false(ntrials,kfold);
                c = cvpartition(ntrials,"KFold",kfold);
                for k = 1:kfold
                    testtrials{typi}(:,k) = test(c,k);
                end
            end
        end

        Nneurons = size(spkcntlmlvs{1},2);

        % filter out neurons whose FRtrainavg could be zero
        neuval = true(Nneurons,1);
        if filtervalneu
            for itt = 1:Ntt
                typi = hireptt==testt(itt);
                tempvalneu = mean(spkcntlmlvs{typi}>0,1) >= 1/kfold;
                neuval(~tempvalneu) = false;
            end
        end

        bayesimage = struct();
        bayesimage.neuval = neuval;
        bayesimage.postprob = cell(numel(hireptt),kfold);
        bayesimage.postprobnorm = cell(numel(hireptt),kfold);
        bayesimage.trainacc = zeros(Ntt,Ntt,kfold);
        bayesimage.testacc = zeros(Ntt,Ntt,kfold);
        bayesimage.infperf = zeros(numel(inferencett),Ntt,kfold);
        for k = 1:kfold

            FRtrainavg = NaN(nnz(neuval), Ntt);
            for itt = 1:Ntt
                typi = hireptt==testt(itt);
                FRtrainavg(:,itt) = (1/twin) * mean(spkcntlmlvs{typi}(~testtrials{typi}(:,k), neuval),1);
            end

            for typi = 1:numel(hireptt)
                ntrials = size(spkcntlmlvs{typi},1);
                bayesimage.postprob{typi,k} = NaN( ntrials, Ntt ); % intialize
                for itrial = 1:ntrials
                    spkcntvec = spkcntlmlvs{typi}(itrial,neuval)';
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
        naivegauss.spkcntoffset = spkcntoffset;
        naivegauss.postprob = cell(numel(hireptt),kfold);
        naivegauss.postprobnorm = cell(numel(hireptt),kfold);
        naivegauss.trainacc = zeros(Ntt,Ntt,kfold);
        naivegauss.testacc = zeros(Ntt,Ntt,kfold);
        naivegauss.infperf = zeros(numel(inferencett),Ntt,kfold);
        for k = 1:kfold
            for typi = 1:numel(hireptt)
                ntrials = size(spkcntlmlvs{typi},1);
                ratelikelihood = zeros(Nneurons, ntrials, Ntt);
                for jtt = 1:Ntt
                    typj = hireptt==testt(jtt);
                    temptesttrials = testtrials{typj}(:,k);
                    for ci = 1:Nneurons
                        tempspk = spkcntoffset+spkcntlmlvs{typj}(~temptesttrials, ci);
                        xval = log10(tempspk(tempspk>0));
                        mu = mean(xval);
                        sigma = std(xval);
                        xtt = log10(spkcntoffset+spkcntlmlvs{typi}(:,ci));
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
                tempvalneu = mean(spkcntlmlvs{typi}>0,1) >= 1/kfold;
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
        mvnbayes.spkcntoffset = spkcntoffset;
        mvnbayes.neuval = neuval;
        mvnbayes.postprob = cell(numel(hireptt),kfold);
        mvnbayes.postprobnorm = cell(numel(hireptt),kfold);
        mvnbayes.trainacc = zeros(Ntt,Ntt,kfold);
        mvnbayes.testacc = zeros(Ntt,Ntt,kfold);
        mvnbayes.infperf = zeros(numel(inferencett),Ntt,kfold);
        try
            for k = 1:kfold
                for typi = 1:numel(hireptt)
                    ntrials = size(spkcntlmlvs{typi},1);
                    mvnbayes.postprob{typi,k} = NaN( ntrials, Ntt );
                    for jtt = 1:Ntt
                        typj = hireptt==testt(jtt);
                        temptesttrials = testtrials{typj}(:,k);

                        tempspkn0 = spkcntoffset+spkcntlmlvs{typi}(~temptesttrials,neuval);
                        tempspkn0(tempspkn0<=0) = NaN;
                        X = log10(tempspkn0);
                        mu = nanmean(X,1);
                        Sigma = cov(X, 'partialrows');
                        xtt = log10( spkcntoffset+spkcntlmlvs{typi}(:,neuval) );
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
                            eigVal(eye(size(eigVal)) & eigVal<=1e-10) = 1e-10; % Shift non-positive eigenvalues
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
        catch
            mvnbayes.trainacc = NaN(Ntt,Ntt,kfold);
            mvnbayes.testacc = NaN(Ntt,Ntt,kfold);
            mvnbayes.infperf = NaN(numel(inferencett),Ntt,kfold);
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
        pcanaivegauss.logspkopt = logspkopt;
        pcanaivegauss.postprob = cell(numel(hireptt),kfold);
        pcanaivegauss.postprobnorm = cell(numel(hireptt),kfold);
        pcanaivegauss.trainacc = zeros(Ntt,Ntt,kfold);
        pcanaivegauss.testacc = zeros(Ntt,Ntt,kfold);
        pcanaivegauss.infperf = zeros(numel(inferencett),Ntt,kfold);
        for k = 1:kfold
            for typi = 1:numel(hireptt)
                ntrials = size(spkcntlmlvs{typi},1);

                ratelikelihood = zeros(Nneurons, ntrials, Ntt);
                for jtt = 1:Ntt
                    typj = hireptt==testt(jtt);
                    temptesttrials = testtrials{typj}(:,k);
                    if logspkopt
                        tempspk = spkcntoffset+spkcntlmlvs{typj}(~temptesttrials, :);
                        Xval = log10(tempspk(tempspk>0));
                        Xtt = log10(spkcntoffset+spkcntlmlvs{typi});
                    else
                        Xval = spkcntlmlvs{typj}(~temptesttrials, :);
                        Xtt = spkcntlmlvs{typi};
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
        fprintf('pcanaivegauss test accuracy %.4f\n', mean(diag(mean(pcanaivegauss.testacc,3))))

        %% Naive Bayes PCA Decoder:

        %% Naive Bayes Decoder (independent neurons)
        % contruct P(r|s) based on observed spike counts (train vs test split)

        %% AODE (Sugden...Andermann 2020)
        % acrivity was binarized based on a threshold

        %% UMAP + Bayesian decoding

        if islope==0
            if computeSVM
                svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
                svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
                save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC', '-v7.3')
                save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models', '-v7.3')
            end

            bayesfn = strcat(pathpp, 'bayesinferencedecoding_V1.mat');
            save(bayesfn, 'bayesimage', 'naivegauss', 'mvnbayes', 'pcanaivegauss')
        else
            if islope==1
                if computeSVM
                    SVMtrainICRC_lmlvs = SVMtrainICRC;
                    SVMtrainICRC_models_lmlvs = SVMtrainICRC_models;
                end
                bayesimage_lmlvs = bayesimage;
                naivegauss_lmlvs = naivegauss;
                mvnbayes_lmlvs = mvnbayes;
                pcanaivegauss_lmlvs = pcanaivegauss;
            else
                if computeSVM
                    SVMtrainICRC_lmlvs(islope) = SVMtrainICRC;
                    SVMtrainICRC_models_lmlvs(islope) = SVMtrainICRC_models;
                end
                bayesimage_lmlvs(islope) = bayesimage;
                naivegauss_lmlvs(islope) = naivegauss;
                mvnbayes_lmlvs(islope) = mvnbayes;
                pcanaivegauss_lmlvs(islope) = pcanaivegauss;
            end
        end

    end
    if computeSVM
        svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
        svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
        save(svmfn, 'disperses', 'preproc', 'whichSVMkernel', 'SVMtrainICRC_lmlvs', '-v7.3')
        save(svmmdlfn, 'disperses', 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models_lmlvs', '-v7.3')
    end

    bayesfn = strcat(pathpp, 'bayesinferencedecoding_V1_lmlvslopes.mat');
    save(bayesfn, 'disperses', 'validneurons_lmlvs', 'newspkcnt_lmlvs', ...
        'bayesimage_lmlvs', 'naivegauss_lmlvs', 'mvnbayes_lmlvs', 'pcanaivegauss_lmlvs', '-v7.3')

    toc(sesclk)
end

warning('on')

%% run just SVM acros LMLV slopes
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/Analyze_IC_OpenScope_v240130/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end

addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_V1.mat'])

for ises = 1:numel(nwbsessions)
    clearvars -except ises nwbsessions spkcntIChiV1agg hireptt lmlvslope lmlvyintercept
    sesclk = tic;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    excludeneuvar0 = 2;
    % if 0, keep all neurons; if 1, exclude 0 variance neurons in train trial
    % types; if 2 exclude 0 variance neurons in all trial types
    twin = 0.4; % in sec
    testt = [106 107 110 111];
    inferencett = [1105 1109];

    disperses = -1:0.5:2;
    Ntt = numel(testt);
    Nhireptt = numel(hireptt);

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

    SVMtrainICRC_lmlvs = struct();
    SVMtrainICRC_models_lmlvs = struct();
    for islope = 0:numel(disperses)
        if islope==0
            spkcntlmlvs = spkcntIChiV1agg{ises};
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

            % note, fitcsvm removes entire rows of data corresponding to a missing response. 
            % rows correspond to observations
            % When computing total weights (see the next bullets), 
            % fitcsvm ignores any weight corresponding to an observation 
            % with at least one missing predictor. 
            % This action can lead to unbalanced prior probabilities in 
            % balanced-class problems. Consequently, observation box constraints might not equal BoxConstraint.
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

            spkcntlmlvs = cell(Nhireptt,1);
            for n = 1:Nhireptt
                spkcntlmlvs{n} = squeeze(newspkcntICtt(:,n,valneu));
            end

            tempR = reshape(newspkcntICtt(:,:,valneu), Nrep*Nhireptt, nnz(valneu))';
            trialorder = reshape( repmat(hireptt,Nrep,1), 1,[]);
        end
        % typi = 4;
        % figure; hold all
        % plot(log(mean(spkcntIChiV1agg{ises}{typi},1)), log(var(spkcntIChiV1agg{ises}{typi},0,1)), 'k.')
        % plot(log(mean(spkcntlmlvs{typi},1)), log(var(spkcntlmlvs{typi},0,1)), 'r.')

        %% comparison point: linear SVM ~7min per slope
        load(['S:\OpenScopeData\00248_v240130\SVM_trainICRC_selectareas\' nwbsessions{ises} '\SVM_trainICRC_VISpRS_Linear_zscore_ICwcfg1.mat'], 'SVMtrainICRC')
        kfold = size(SVMtrainICRC.spkcnt.testtrialinds,2);
        traintrials = false(length(trialorder),kfold);
        testtrials = false(length(trialorder),kfold);
        for itt = 1:Ntt
            trialsoind = find(SVMtrainICRC.trialorder==testt(itt));
            typi = hireptt==testt(itt);
            ntrials = size(spkcntlmlvs{typi},1);
            temptesttrials = false(ntrials,kfold);
            for k = 1:kfold
                temptesttrials(:,k) = ismember(trialsoind, SVMtrainICRC.spkcnt.testtrialinds(:,k));
            end
            traintrials(trialorder==testt(itt),:) = ~temptesttrials;
            testtrials(trialorder==testt(itt),:) = temptesttrials;
        end
        traintrialinds = zeros(max(sum(traintrials,1)), kfold);
        testtrialinds = zeros(max(sum(testtrials,1)), kfold);
        for k = 1:kfold
            traintrialinds(:,k) = find(traintrials(:,k));
            testtrialinds(:,k) = find(testtrials(:,k));
        end

        preproc = 'meancenter';
        whichSVMkernel = 'Linear';
        svmdesc = 'trainICRC';

        cvtrials = struct();
        cvtrials.loadcvpartition = true;
        cvtrials.traintrialinds = traintrialinds;
        cvtrials.testtrialinds = testtrialinds;
        % SVMtrainICRC_models is 14 MB
        [SVMtrainICRC, SVMtrainICRC_models] = computeICtxi_SVM(tempR, trialorder, ...
            svmdesc, 'spkcnt', preproc, whichSVMkernel, cvtrials);
        if islope>0 %&& excludeneuvar0>0
            SVMtrainICRC.valneu = valneu;
        end

        if ~isequal(SVMtrainICRC.trialtypes, testt)
            error('mismatch in test trial types: check that you loaded trainICRC')
        end

        % train accuracy
        trainlabs = SVMtrainICRC.trialorder(SVMtrainICRC.spkcnt.traintrialinds);
        trainpred = SVMtrainICRC.spkcnt.train.label;
        trainacc = zeros(Ntt);
        for itt = 1:Ntt
            trialsoi = trainlabs==testt(itt);
            [v,c]=uniquecnt(trainpred(trialsoi));
            trainacc(itt, ismember(testt,v)) = c/nnz(trialsoi);
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
            [v,c]=uniquecnt(infpred(trialsoi,:));
            infperf(itt, ismember(testt,v)) = c/(size(infpred,2)*nnz(trialsoi));
        end

        if islope==0
        fprintf('%d/%d %s as-is\n', ises, numel(nwbsessions), nwbsessions{ises})
        else
        fprintf('%d/%d %s LMLV slope %.2f\n', ises, numel(nwbsessions), nwbsessions{ises}, disperses(islope))
        end
        disp('SVM trainacc')
        disp(mean(trainacc,3))
        disp('SVM testacc')
        disp(mean(testacc,3))
        disp('SVM infperf')
        disp(mean(infperf,3))

        if islope==0
            svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '.mat');
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models', '-v7.3')
        else
            if islope==1
                SVMtrainICRC_lmlvs = SVMtrainICRC;
                SVMtrainICRC_models_lmlvs = SVMtrainICRC_models;
            else
                SVMtrainICRC_lmlvs(islope) = SVMtrainICRC;
                SVMtrainICRC_models_lmlvs(islope) = SVMtrainICRC_models;
            end
        end

    end
    svmfn = strcat(pathpp, 'SVM_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
    svmmdlfn = strcat(pathpp, 'SVMmodels_', svmdesc, '_V1_', whichSVMkernel, '_', preproc, '_lmlvslopes.mat');
    save(svmfn, 'disperses', 'preproc', 'whichSVMkernel', 'SVMtrainICRC_lmlvs', '-v7.3')
    save(svmmdlfn, 'disperses', 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models_lmlvs', '-v7.3')

    toc(sesclk)
end

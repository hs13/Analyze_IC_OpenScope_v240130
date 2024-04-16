%% RESUME HERE: EDIT CODE FROM trainICRC to trainREx
discardbelowNneurons = 50;
whichICblock = 'ICwcfg1';
SVMall = SVMtrainRExagg.(whichICblock);
Nsessions = numel(SVMall);

visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
traintrialtypes = SVMall(1).VISp.trialtypes;

%% test trials: check that there is correlation between area decoders
% this correlation must exceed null distribution, where null distribution
% is obtained by shuffling within trial types

% perhaps train and test trial divide should have been the same between all
% areas...
Nshuf = 1000;

testaccuracy = NaN(numel(visareas), Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    testtrialstt = SVMall(ises).(whichvisarea).trialorder(SVMall(ises).(whichvisarea).spkcnt.testtrialinds);
    testaccuracy(a,ises) = mean(SVMall(ises).(whichvisarea).spkcnt.test.label == testtrialstt, 'all');
end
end
% disp(testaccuracy)

testpredmatchchance = NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
testpredmatchchance(:,:,ises) = testaccuracy(:,ises) * testaccuracy(:,ises)';
end

tic
testpredmatch = NaN(numel(visareas), numel(visareas), Nsessions);
testpredmatchnull = NaN(numel(visareas), numel(visareas), Nshuf, Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    [ttiA,testtrialorderA]=sort(SVMall(ises).(whichvisareaA).spkcnt.testtrialinds(:));
    testlabelA = SVMall(ises).(whichvisareaA).spkcnt.test.label(testtrialorderA);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        if Nneuronsperarea(ises,b)<discardbelowNneurons
            continue
        end
        [ttiB,testtrialorderB]=sort(SVMall(ises).(whichvisareaB).spkcnt.testtrialinds(:));
        if ~isequal(ttiA, ttiB)
            error('set of test trial inds do not match between areas')
        end
        testlabelB = SVMall(ises).(whichvisareaB).spkcnt.test.label(testtrialorderB);
        testtrialsttB = SVMall(ises).(whichvisareaB).trialorder(ttiB);
        testpredmatch(a,b,ises) = mean(testlabelA==testlabelB);
        
        testaccB = mean(SVMall(ises).(whichvisareaB).spkcnt.test.label==SVMall(ises).(whichvisareaB).trialorder(SVMall(ises).(whichvisareaB).spkcnt.testtrialinds), 'all');
        if mean(testlabelB(:)==testtrialsttB(:)) ~= testaccB
            error('check code')
        end
        
        for ishuf = 1:Nshuf
            testlabelBshuf = testlabelB(:);
            for itt = 1:numel(traintrialtypes)
                trialsoind = find(testtrialsttB(:)==traintrialtypes(itt));
                trialsoindperm = trialsoind(randperm(numel(trialsoind)));
                testlabelBshuf(trialsoind) = testlabelB(trialsoindperm);
            end
            % isequal(mean(testlabelBshuf==testtrialsttB(:)), mean(testlabelB(:)==testtrialsttB(:)))
            testpredmatchnull(a,b,ishuf,ises) = mean(testlabelA==testlabelBshuf);
        end
    end
end
end
toc

testpredmatchprctile =  NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
    for a = 1:numel(visareas)
        for b = a+1:numel(visareas)
            if ~isnan(testpredmatch(a,b,ises))
                testpredmatchprctile(a,b,ises) = 100*mean(testpredmatchnull(a,b,:,ises)<testpredmatch(a,b,ises),3);
            end
        end
    end
end


ICtrialtypes = [106, 111];
probehc2 = struct();
probeNtrials = struct();
probehc2val = struct(); % invalidate trials where there was a tie for mode predictions across K-fold cross-validation
probeNvaltrials = struct();
for ab = 1:2
switch ab
    case 1
whichvisareaA = 'VISp';
whichvisareaB = 'VISl';
    case 2
whichvisareaA = 'VISp';
whichvisareaB = 'VISal';
    otherwise
        error('specify whichvisareaA and whichvisareaB')
end
a = find(strcmp(visareas, whichvisareaA));
b = find(strcmp(visareas, whichvisareaB));
if ~( numel(a)==1 && numel(b)==1 )
    error('check whichvisareaA and whichvisareaB')
end
ABfield = [whichvisareaA '_' whichvisareaB];

ttbe = 0.5*([traintrialtypes(1)-1 traintrialtypes] + [traintrialtypes traintrialtypes(end)+1]);
probehc2.(ABfield) = NaN(length(traintrialtypes), length(traintrialtypes), length(ICtrialtypes), Nsessions);
probeNtrials.(ABfield) = zeros(length(ICtrialtypes), Nsessions);
probehc2val.(ABfield) = NaN(length(traintrialtypes), length(traintrialtypes), length(ICtrialtypes), Nsessions);
probeNvaltrials.(ABfield) = zeros(length(ICtrialtypes), Nsessions);
for ises = 1:Nsessions
    if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
        continue
    end
for iprobe = 1:length(ICtrialtypes)
    probetrialsA = SVMall(ises).(whichvisareaA).trialorder==ICtrialtypes(iprobe);
    [probepredA,FA,CA] = mode( SVMall(ises).(whichvisareaA).spkcnt.all.label(probetrialsA,:),2 );    
    validtrialsA = cellfun(@numel, CA)==1;
    % figure; hold all; histogram(FA,-0.5:10.5); histogram(FA(~validtrialsA),-0.5:10.5)

    probetrialsB = SVMall(ises).(whichvisareaB).trialorder==ICtrialtypes(iprobe);
    [probepredB,FB,CB] = mode( SVMall(ises).(whichvisareaB).spkcnt.all.label(probetrialsB,:),2 );
    validtrialsB = cellfun(@numel, CB)==1;
    
    % rows are probepredA
    hc = histcounts2(probepredA, probepredB, 'XBinEdges', ttbe, 'YBinEdges', ttbe, 'normalization', 'probability');
    probehc2.(ABfield)(:,:,iprobe,ises) = hc;
probeNtrials.(ABfield)(iprobe,ises) = length(probepredA);

    validtrials = validtrialsA & validtrialsB;
    hcval = histcounts2(probepredA(validtrials), probepredB(validtrials), 'XBinEdges', ttbe, 'YBinEdges', ttbe, 'normalization', 'probability');
    probehc2val.(ABfield)(:,:,iprobe,ises) = hcval;
probeNvaltrials.(ABfield)(iprobe,ises) = nnz(validtrials);
end
end
end

% match null distribution for inference trials:
% calculate for each trial type
% use mode prediction across K-fold cross-validation (i.e., Nsplits)
ICpredmatch = cell(size(ICtrialtypes));
ICpredmatchnull = cell(size(ICtrialtypes));
ICpredmatchprctile = cell(size(ICtrialtypes));
tic
for ii = 1:numel(ICtrialtypes)
ICpredmatch{ii} = NaN(numel(visareas), numel(visareas), Nsessions);
ICpredmatchnull{ii} = NaN(numel(visareas), numel(visareas), Nshuf, Nsessions);
ICpredmatchprctile{ii} =  NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    ICtrials = SVMall(ises).(whichvisareaA).trialorder==ICtrialtypes(ii);
    IClabelA = SVMall(ises).(whichvisareaA).spkcnt.all.label(ICtrials);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        if Nneuronsperarea(ises,b)<discardbelowNneurons
            continue
        end
        if ~isequal( ICtrials, SVMall(ises).(whichvisareaB).trialorder==ICtrialtypes(ii))
            error('trialorder should be same across all decoders within same session')
        end
        IClabelB = SVMall(ises).(whichvisareaB).spkcnt.all.label(ICtrials);
        ICpredmatch{ii}(a,b,ises) = mean(IClabelA==IClabelB);
        
        
        for ishuf = 1:Nshuf
            IClabelBshuf = IClabelB(randperm(numel(IClabelB)));
            ICpredmatchnull{ii}(a,b,ishuf,ises) = mean(IClabelA==IClabelBshuf);
        end

        ICpredmatchprctile{ii}(a,b,ises) = 100*mean(ICpredmatchnull{ii}(a,b,:,ises)<ICpredmatch{ii}(a,b,ises),3);
        
    end
end
end
toc
end

XREtrialtypes = [1201 1299];
ICasXREmatch = cell(size(ICtrialtypes));
ICasXREmatchnull = cell(size(ICtrialtypes));
ICasXREmatchprctile = cell(size(ICtrialtypes));
for ii = 1:numel(ICtrialtypes)
tic
ICasXREmatch{ii} = NaN(numel(visareas), numel(visareas), Nsessions);
ICasXREmatchnull{ii} = NaN(numel(visareas), numel(visareas), Nshuf, Nsessions);
ICasXREmatchprctile{ii} =  NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    ICtrials = SVMall(ises).(whichvisareaA).trialorder==ICtrialtypes(ii);
    IClabelA = SVMall(ises).(whichvisareaA).spkcnt.all.label(ICtrials)==XREtrialtypes(ii);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        if Nneuronsperarea(ises,b)<discardbelowNneurons
            continue
        end
        if ~isequal( ICtrials, SVMall(ises).(whichvisareaB).trialorder==ICtrialtypes(ii))
            error('trialorder should be same across all decoders within same session')
        end
        IClabelB = SVMall(ises).(whichvisareaB).spkcnt.all.label(ICtrials)==XREtrialtypes(ii);
        ICasXREmatch{ii}(a,b,ises) = mean(IClabelA & IClabelB);
        
        for ishuf = 1:Nshuf
            IClabelBshuf = IClabelB(randperm(numel(IClabelB)));
            ICasXREmatchnull{ii}(a,b,ishuf,ises) = mean(IClabelA & IClabelBshuf);
        end
        ICasXREmatchprctile{ii}(a,b,ises) = 100*mean(ICasXREmatchnull{ii}(a,b,:,ises)<ICasXREmatch{ii}(a,b,ises),3);
    end
end
end
toc
end
isequaln(ICpredmatch, ICasXREmatch)
figure; plot(ICpredmatch{ii}(:), ICasXREmatch{ii}(:), 'o')

%% plot test trial prediction match between pairs of areas
figure
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        subplot(numel(visareas), numel(visareas), (a-1)*numel(visareas)+b)
        hold all
        plot(testaccuracy(a,:), testaccuracy(b,:), 'o')
        xl = xlim;
        plot(xl,xl,'r-')
        % xlabel([whichvisareaA ' test accuracy'])
        % ylabel([whichvisareaB ' test accuracy'])
        xlabel(whichvisareaA)
        ylabel(whichvisareaB)
        p = signrank(testaccuracy(a,:), testaccuracy(b,:));
        title(sprintf('test accuracy p=%.4f',p))
    end
end

figure; hold all
for ises = 1:Nsessions
plot(reshape(testpredmatchchance(:,:,ises),[],1), reshape(nanmean(testpredmatchnull(:,:,:,ises),3),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('chance (if independent)')
ylabel('null distribution mean')

figure; hold all
for ises = 1:Nsessions
plot(reshape(nanmean(testpredmatchnull(:,:,:,ises),3),[],1), reshape(testpredmatch(:,:,ises),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')

figure; 
subplot(2,2,1)
imagesc(squeeze(nanmean(testpredmatch,3))); 
set(gca, 'XTick', 1:size(testpredmatch,2), 'XTickLabel', visarealabels, ...
    'YTick', 1:size(testpredmatch,1), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title('test trial prediction match')
subplot(2,2,2)
hold all
tempmat = squeeze(nanmean(testpredmatchprctile,3));
imagesc(tempmat); 
for a = 1:numel(visareas)
    for b = a+1:numel(visareas)
        text(b,a, sprintf('%.0f',tempmat(a,b)), 'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    end
end
axis([0.5 numel(visareas)+0.5 0.5 numel(visareas)+0.5])
set(gca, 'XTick', 1:size(testpredmatchprctile,2), 'XTickLabel', visarealabels, ...
    'YTick', 1:size(testpredmatchprctile,1), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title({'test trial prediction match', 'percentile w.r.t. null distribution'})
subplot(2,2,3)
imagesc(squeeze(nanmean(testpredmatch-squeeze(mean(testpredmatchnull,3)),3))); 
set(gca, 'XTick', 1:size(testpredmatch,2), 'XTickLabel', visarealabels, ...
    'YTick', 1:size(testpredmatch,1), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title('test trial prediction match - shuf. dist.')
subplot(2,2,4)
imagesc(squeeze(nanmean(testpredmatch-testpredmatchchance,3))); 
set(gca, 'XTick', 1:size(testpredmatch,2), 'XTickLabel', visarealabels, ...
    'YTick', 1:size(testpredmatch,1), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title('test trial prediction match - chance')

disp('V1 & LM testpredmatchprctile')
disp(squeeze(testpredmatchprctile(1,2,:)))

disp('V1 & AL testpredmatchprctile')
disp(squeeze(testpredmatchprctile(1,4,:)))

for a=1:numel(visareas)
    for b=a+1:numel(visareas)
ptestpred=signrank(reshape(nanmean(testpredmatchnull(a,b,:,:),3),[],1), reshape(testpredmatch(a,b,:),[],1));
fprintf('%s & %s test match N=%d p=%.4f\n', visarealabels{a}, visarealabels{b}, ...
    nnz(~isnan(testpredmatch(a,b,:))), ptestpred)
    end
end

figure; 
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};

        subplot(numel(visareas), numel(visareas), (a-1)*numel(visareas)+b)
        hold all
plot(reshape(nanmean(testpredmatchnull(a,b,:,:),3),[],1), reshape(testpredmatch(a,b,:),[],1), 'o')
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
ptestpred=signrank(reshape(nanmean(testpredmatchnull(a,b,:,:),3),[],1), reshape(testpredmatch(a,b,:),[],1));
title(sprintf('%s vs %s test match N=%d p=%.4f\n', visarealabels{a}, visarealabels{b}, ...
    nnz(~isnan(testpredmatch(a,b,:))), ptestpred))
    end
end

% V1 & LM test match N=8 p=0.0391
% V1 & RL test match N=11 p=0.0830
% V1 & AL test match N=10 p=0.0645
% V1 & PM test match N=12 p=0.0923
% V1 & AM test match N=11 p=0.0049
% LM & RL test match N=7 p=0.0156
% LM & AL test match N=6 p=0.0312
% LM & PM test match N=8 p=0.3125
% LM & AM test match N=7 p=0.0781
% RL & AL test match N=9 p=0.0742
% RL & PM test match N=11 p=0.0186
% RL & AM test match N=10 p=0.0137
% AL & PM test match N=10 p=0.0371
% AL & AM test match N=9 p=0.0195
% PM & AM test match N=11 p=0.0029

%% inference trials:
figure
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    tempprobe = squeeze(mean(HR_SVMtrainREx.(whichICblock).(whichvisareaA).probe,3));
    infscoreA = squeeze( (tempprobe(1,1,:)+tempprobe(4,2,:))/2 - 0.5 );
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
    tempprobe = squeeze(mean(HR_SVMtrainREx.(whichICblock).(whichvisareaB).probe,3));
    infscoreB = squeeze( (tempprobe(1,1,:)+tempprobe(4,2,:))/2 - 0.5 );

        subplot(numel(visareas), numel(visareas), (a-1)*numel(visareas)+b)
        hold all
        plot(infscoreA, infscoreB, 'o')
        xl = xlim;
        plot(xl,xl,'r-')
        % xlabel([whichvisareaA ' test accuracy'])
        % ylabel([whichvisareaB ' test accuracy'])
        xlabel(whichvisareaA)
        ylabel(whichvisareaB)
        p = signrank(testaccuracy(a,:), testaccuracy(b,:));
        title(sprintf('inference score p=%.4f',p))
    end
end

figure;
for ab = 1:2
    switch ab
        case 1
            whichvisareaA = 'VISp';
            whichvisareaB = 'VISl';
        case 2
            whichvisareaA = 'VISp';
            whichvisareaB = 'VISal';
        otherwise
            error('specify whichvisareaA and whichvisareaB')
    end
    ABfield = [whichvisareaA '_' whichvisareaB];
    for ii = 1:numel(ICtrialtypes)
        subplot(2,2,2*(ab-1)+ii)
        imagesc(squeeze(nanmean(probehc2val.(ABfield)(:,:,ii,:),4)))
set(gca, 'XTick', 1:numel(traintrialtypes), 'XTickLabel', traintrialtypes, ...
    'YTick', 1:numel(traintrialtypes), 'YTickLabel', traintrialtypes)
        xlabel(whichvisareaB)
        ylabel(whichvisareaA)
        title(sprintf('trial %d', ICtrialtypes(ii) ))
caxis([0 0.25])
    end
end

% ICpredmatch
figure; 
for ii = 1:numel(ICtrialtypes)
subplot(2,3,3*(ii-1)+1)
imagesc(squeeze(nanmean(ICpredmatch{ii},3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('I_C_%d trial prediction match', ii))

subplot(2,3,3*(ii-1)+2)
imagesc(squeeze(nanmean(ICpredmatch{ii}-squeeze(mean(ICpredmatchnull{ii},3)),3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('I_C_%d trial prediction match - shuf. dist.', ii))
subplot(2,3,3*(ii-1)+3)
hold all
tempmat = squeeze(nanmean(ICpredmatchprctile{ii},3));
imagesc(tempmat); 
for a = 1:numel(visareas)
    for b = a+1:numel(visareas)
        text(b,a, sprintf('%.0f',tempmat(a,b)), 'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    end
end
axis([0.5 numel(visareas)+0.5 0.5 numel(visareas)+0.5])
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
caxis([50 100])
colorbar
title({sprintf('I_C_%d trial prediction match', ii), 'percentile w.r.t. null distribution'})
end

% not significant -- inference predictions does not match across areas,
% which argues against global attractive dynamics...
disp('V1 & LM ICpredmatchprctile')
disp([squeeze(ICpredmatchprctile{1}(1,2,:)) squeeze(ICpredmatchprctile{2}(1,2,:))])

disp('V1 & AL ICpredmatchprctile')
disp([squeeze(ICpredmatchprctile{1}(1,4,:)) squeeze(ICpredmatchprctile{2}(1,4,:))])

figure; 
for ii = 1:numel(ICtrialtypes)
subplot(2,3,3*(ii-1)+1)
imagesc(squeeze(nanmean(ICasXREmatch{ii},3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('I_C_%d as X_R_E_%d prediction match', ii, ii))

subplot(2,3,3*(ii-1)+2)
imagesc(squeeze(nanmean(ICasXREmatch{ii}-squeeze(mean(ICasXREmatchnull{ii},3)),3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('I_C_%d as X_R_E_%d prediction match - shuf. dist.', ii,ii))
subplot(2,3,3*(ii-1)+3)
hold all
tempmat = squeeze(nanmean(ICasXREmatchprctile{ii},3));
imagesc(tempmat); 
for a = 1:numel(visareas)
    for b = a+1:numel(visareas)
        text(b,a, sprintf('%.0f',tempmat(a,b)), 'color', 'k', 'horizontalalignment', 'center', 'verticalalignment', 'middle')
    end
end
axis([0.5 numel(visareas)+0.5 0.5 numel(visareas)+0.5])
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
caxis([50 100])
colorbar
title({sprintf('I_C_%d as X_R_E_%d trial prediction match', ii,ii), 'percentile w.r.t. null distribution'})
end

% not significant
disp('V1 & LM ICasXREmatchprctile')
disp([squeeze(ICasXREmatchprctile{1}(1,2,:)) squeeze(ICasXREmatchprctile{2}(1,2,:))])

disp('V1 & AL ICasXREmatchprctile')
disp([squeeze(ICasXREmatchprctile{1}(1,4,:)) squeeze(ICasXREmatchprctile{2}(1,4,:))])

figure; 
for ii = 1:numel(ICtrialtypes)
subplot(2,2,ii)
hold all
for ises = 1:Nsessions
plot(reshape(nanmean(ICpredmatchnull{ii}(:,:,:,ises),3),[],1), reshape(ICpredmatch{ii}(:,:,ises),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
xvec = reshape(nanmean(ICpredmatchnull{ii},3),[],1);
yvec = reshape(ICpredmatch{ii},[],1);
p=signrank(xvec(~isnan(xvec) & ~isnan(yvec)), yvec(~isnan(xvec) & ~isnan(yvec)));
% figure; plot(xvec,yvec,'o');hold on; xl-xlim;plot(xl,xl,'r-')
p=signrank(reshape(nanmean(ICpredmatchnull{ii},3),[],1), reshape(ICpredmatch{ii},[],1));
title(sprintf('I_C_%d pred p=%.4f', ii,p))

subplot(2,2,2+ii)
hold all
for ises = 1:Nsessions
plot(reshape(nanmean(ICasXREmatchnull{ii}(:,:,:,ises),3),[],1), reshape(ICasXREmatch{ii}(:,:,ises),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
p = signrank(reshape(nanmean(ICasXREmatchnull{ii},3),[],1), reshape(ICasXREmatch{ii},[],1));
title(sprintf('I_C_%d as X_R_E_%d p=%.4f', ii,ii,p))
end

for ab = 1:2
    switch ab
        case 1
a=1;b=2;
        case 2
a=1;b=4;
    end
pTREpred=signrank(reshape(nanmean(ICpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICpredmatch{ii}(a,b,:),[],1));
pTREasIC=signrank(reshape(nanmean(ICasXREmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICasXREmatch{ii}(a,b,:),[],1));
fprintf('%s & %s match IC p=%.4f IC as XRE p=%.4f\n', visarealabels{a}, visarealabels{b}, pTREpred, pTREasIC)
end

for a=1:numel(visareas)
    for b=a+1:numel(visareas)        
pTREpred=signrank(reshape(nanmean(ICpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICpredmatch{ii}(a,b,:),[],1));
pTREasIC=signrank(reshape(nanmean(ICasXREmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICasXREmatch{ii}(a,b,:),[],1));
fprintf('%s & %s match N=%d IC p=%.4f IC as XRE p=%.4f\n', visarealabels{a}, visarealabels{b}, ...
    nnz(~isnan(ICpredmatch{ii}(a,b,:))), pTREpred, pTREasIC)
    end
end

figure
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};

        subplot(numel(visareas), numel(visareas), (a-1)*numel(visareas)+b)
hold all
plot(reshape(nanmean(ICpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICpredmatch{ii}(a,b,:),[],1), 'o')
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
p=signrank(reshape(nanmean(ICpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICpredmatch{ii}(a,b,:),[],1));
title(sprintf('%s vs %s I_C_%d pred p=%.4f', whichvisareaA, whichvisareaB, ii,p))

        subplot(numel(visareas), numel(visareas), (b-1)*numel(visareas)+a)
hold all
plot(reshape(nanmean(ICasXREmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICasXREmatch{ii}(a,b,:),[],1), 'o')
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
p = signrank(reshape(nanmean(ICasXREmatchnull{ii}(a,b,:,:),3),[],1), reshape(ICasXREmatch{ii}(a,b,:),[],1));
title(sprintf('%s vs %s I_C_%d as X_R_E_%d p=%.4f', whichvisareaA, whichvisareaB, ii,ii,p))        

    end
end


% *V1 & LM match N=8 IC p=0.0078 IC as XRE p=0.0078
% V1 & RL match N=11 IC p=0.4648 IC as XRE p=0.5195
% #V1 & AL match N=10 IC p=0.0840 IC as XRE p=0.0840
% V1 & PM match N=12 IC p=0.9697 IC as XRE p=0.9697
% V1 & AM match N=11 IC p=0.2061 IC as XRE p=0.2061
% LM & RL match N=7 IC p=0.1562 IC as XRE p=0.2188
% LM & AL match N=6 IC p=0.5625 IC as XRE p=0.4375
% LM & PM match N=8 IC p=0.6406 IC as XRE p=0.6406
% LM & AM match N=7 IC p=0.2188 IC as XRE p=0.2969
% RL & AL match N=9 IC p=0.2031 IC as XRE p=0.2031
% *RL & PM match N=11 IC p=0.0186 IC as XRE p=0.0186
% *RL & AM match N=10 IC p=0.0371 IC as XRE p=0.0371
% AL & PM match N=10 IC p=0.6250 IC as XRE p=0.5566
% AL & AM match N=9 IC p=0.2500 IC as XRE p=0.2500
% PM & AM match N=11 IC p=0.4648 IC as XRE p=0.4648

% CONCLUSION: ASSESSED WHETHER IC->XRE INFERENCE IS CORRELATED BETWEEN AREAS (focus on V1 & LM, V1 & AL)
% within-session statistics: match percent is not significantly different shuffled distribution
% across-session statistics: match is significantly higher than chance across sessions
% despite being significant across sessions, the difference between actual
% match and chance match is only very slightly different


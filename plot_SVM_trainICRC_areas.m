discardbelowNneurons = 50;
whichICblock = 'ICwcfg1';
SVMall = SVMtrainICRCagg.(whichICblock);
Nsessions = numel(SVMall);

visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};

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


TREtrialtypes = [1105, 1109];
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
probehc2.(ABfield) = NaN(length(traintrialtypes), length(traintrialtypes), length(TREtrialtypes), Nsessions);
probeNtrials.(ABfield) = zeros(length(TREtrialtypes), Nsessions);
probehc2val.(ABfield) = NaN(length(traintrialtypes), length(traintrialtypes), length(TREtrialtypes), Nsessions);
probeNvaltrials.(ABfield) = zeros(length(TREtrialtypes), Nsessions);
for ises = 1:Nsessions
    if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
        continue
    end
for iprobe = 1:length(TREtrialtypes)
    probetrialsA = SVMall(ises).(whichvisareaA).trialorder==TREtrialtypes(iprobe);
    [probepredA,FA,CA] = mode( SVMall(ises).(whichvisareaA).spkcnt.all.label(probetrialsA,:),2 );    
    validtrialsA = cellfun(@numel, CA)==1;
    % figure; hold all; histogram(FA,-0.5:10.5); histogram(FA(~validtrialsA),-0.5:10.5)

    probetrialsB = SVMall(ises).(whichvisareaB).trialorder==TREtrialtypes(iprobe);
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
TREtrialtypes = [1105, 1109];
TREpredmatch = cell(size(TREtrialtypes));
TREpredmatchnull = cell(size(TREtrialtypes));
TREpredmatchprctile = cell(size(TREtrialtypes));
tic
for ii = 1:numel(TREtrialtypes)
TREpredmatch{ii} = NaN(numel(visareas), numel(visareas), Nsessions);
TREpredmatchnull{ii} = NaN(numel(visareas), numel(visareas), Nshuf, Nsessions);
TREpredmatchprctile{ii} =  NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    TREtrials = SVMall(ises).(whichvisareaA).trialorder==TREtrialtypes(ii);
    TRElabelA = SVMall(ises).(whichvisareaA).spkcnt.all.label(TREtrials);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        if Nneuronsperarea(ises,b)<discardbelowNneurons
            continue
        end
        if ~isequal( TREtrials, SVMall(ises).(whichvisareaB).trialorder==TREtrialtypes(ii))
            error('trialorder should be same across all decoders within same session')
        end
        TRElabelB = SVMall(ises).(whichvisareaB).spkcnt.all.label(TREtrials);
        TREpredmatch{ii}(a,b,ises) = mean(TRElabelA==TRElabelB);
        
        
        for ishuf = 1:Nshuf
            TRElabelBshuf = TRElabelB(randperm(numel(TRElabelB)));
            TREpredmatchnull{ii}(a,b,ishuf,ises) = mean(TRElabelA==TRElabelBshuf);
        end

        TREpredmatchprctile{ii}(a,b,ises) = 100*mean(TREpredmatchnull{ii}(a,b,:,ises)<TREpredmatch{ii}(a,b,ises),3);
        
    end
end
end
toc
end

TREtrialtypes = [1105, 1109];
I_C_trialtypes = [106 111];
TREasICmatch = cell(size(TREtrialtypes));
TREasICmatchnull = cell(size(TREtrialtypes));
TREasICmatchprctile = cell(size(TREtrialtypes));
for ii = 1:numel(TREtrialtypes)
tic
TREasICmatch{ii} = NaN(numel(visareas), numel(visareas), Nsessions);
TREasICmatchnull{ii} = NaN(numel(visareas), numel(visareas), Nshuf, Nsessions);
TREasICmatchprctile{ii} =  NaN(numel(visareas), numel(visareas), Nsessions);
for ises = 1:Nsessions
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    if Nneuronsperarea(ises,a)<discardbelowNneurons
        continue
    end
    TREtrials = SVMall(ises).(whichvisareaA).trialorder==TREtrialtypes(ii);
    TRElabelA = SVMall(ises).(whichvisareaA).spkcnt.all.label(TREtrials)==I_C_trialtypes(ii);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        if Nneuronsperarea(ises,b)<discardbelowNneurons
            continue
        end
        if ~isequal( TREtrials, SVMall(ises).(whichvisareaB).trialorder==TREtrialtypes(ii))
            error('trialorder should be same across all decoders within same session')
        end
        TRElabelB = SVMall(ises).(whichvisareaB).spkcnt.all.label(TREtrials)==I_C_trialtypes(ii);
        TREasICmatch{ii}(a,b,ises) = mean(TRElabelA & TRElabelB);
        
        for ishuf = 1:Nshuf
            TRElabelBshuf = TRElabelB(randperm(numel(TRElabelB)));
            TREasICmatchnull{ii}(a,b,ishuf,ises) = mean(TRElabelA & TRElabelBshuf);
        end
        TREasICmatchprctile{ii}(a,b,ises) = 100*mean(TREasICmatchnull{ii}(a,b,:,ises)<TREasICmatch{ii}(a,b,ises),3);
    end
end
end
toc
end

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

for a=1:numel(visareas)
    for b=a+1:numel(visareas)
ptestpred=signrank(reshape(nanmean(testpredmatchnull(a,b,:,:),3),[],1), reshape(testpredmatch(a,b,:),[],1));
fprintf('%s & %s test match N=%d p=%.4f\n', visarealabels{a}, visarealabels{b}, ...
    nnz(~isnan(testpredmatch(a,b,:))), ptestpred)
    end
end
% V1 & LM test match N=8 p=0.0078
% V1 & RL test match N=11 p=0.0029
% V1 & AL test match N=10 p=0.0020
% V1 & PM test match N=12 p=0.0005
% V1 & AM test match N=11 p=0.0020
% LM & RL test match N=7 p=0.0312
% LM & AL test match N=6 p=0.0312
% LM & PM test match N=8 p=0.0156
% LM & AM test match N=7 p=0.1094
% RL & AL test match N=9 p=0.0039
% RL & PM test match N=11 p=0.0020
% RL & AM test match N=10 p=0.0098
% AL & PM test match N=10 p=0.0039
% AL & AM test match N=9 p=0.0039
% PM & AM test match N=11 p=0.0049

%% inference trials:
figure
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    tempREt = squeeze(mean(HR_SVMtrainICRC.(whichICblock).(whichvisareaA).REt,3));
    infscoreA = squeeze( (tempREt(1,1,:)+tempREt(2,4,:))/2 - (tempREt(1,2,:)+tempREt(2,3,:))/2 );
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
    tempREt = squeeze(mean(HR_SVMtrainICRC.(whichICblock).(whichvisareaB).REt,3));
    infscoreB = squeeze( (tempREt(1,1,:)+tempREt(2,4,:))/2 - (tempREt(1,2,:)+tempREt(2,3,:))/2 );

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
    for ii = 1:numel(TREtrialtypes)
        subplot(2,2,2*(ab-1)+ii)
        imagesc(squeeze(nanmean(probehc2val.(ABfield)(:,:,ii,:),4)))
set(gca, 'XTick', 1:numel(traintrialtypes), 'XTickLabel', traintrialtypes, ...
    'YTick', 1:numel(traintrialtypes), 'YTickLabel', traintrialtypes)
        xlabel(whichvisareaB)
        ylabel(whichvisareaA)
        title(sprintf('trial %d', TREtrialtypes(ii) ))
caxis([0 0.25])
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
    for ii = 1:numel(TREtrialtypes)
        switch TREtrialtypes(ii)
            case 1105
                AIC_BLC = squeeze(probehc2val.(ABfield)(1,2,1,:));
                ALC_BIC = squeeze(probehc2val.(ABfield)(2,1,1,:));
            case 1109
                AIC_BLC = squeeze(probehc2val.(ABfield)(4,3,2,:));
                ALC_BIC = squeeze(probehc2val.(ABfield)(3,4,2,:));
        end
        subplot(2,2,2*(ab-1)+ii)
        hold all
        plot(AIC_BLC, ALC_BIC, 'o')
        xl = xlim;
        plot(xl,xl, 'r-')
        xlabel(sprintf('%s as IC, %s as LC', whichvisareaA, whichvisareaB))
        ylabel(sprintf('%s as LC, %s as IC', whichvisareaA, whichvisareaB))
        p = signrank(AIC_BLC, ALC_BIC);
        title(sprintf('T_R_E_%d p=%.4f', ii, p))
    end
end

% TREpredmatch
figure; 
for ii = 1:numel(TREtrialtypes)
subplot(2,3,3*(ii-1)+1)
imagesc(squeeze(nanmean(TREpredmatch{ii},3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('T_R_E_%d trial prediction match', ii))

subplot(2,3,3*(ii-1)+2)
imagesc(squeeze(nanmean(TREpredmatch{ii}-squeeze(mean(TREpredmatchnull{ii},3)),3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('T_R_E_%d trial prediction match - shuf. dist.', ii))
subplot(2,3,3*(ii-1)+3)
hold all
tempmat = squeeze(nanmean(TREpredmatchprctile{ii},3));
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
title({sprintf('T_R_E_%d trial prediction match', ii), 'percentile w.r.t. null distribution'})
end

figure; 
for ii = 1:numel(TREtrialtypes)
subplot(2,3,3*(ii-1)+1)
imagesc(squeeze(nanmean(TREasICmatch{ii},3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('T_R_E_%d as I_C_%d prediction match', ii, ii))

subplot(2,3,3*(ii-1)+2)
imagesc(squeeze(nanmean(TREasICmatch{ii}-squeeze(mean(TREasICmatchnull{ii},3)),3))); 
set(gca, 'XTick', 1:numel(visareas), 'XTickLabel', visarealabels, ...
    'YTick', 1:numel(visareas), 'YTickLabel', visarealabels, 'YDir', 'reverse')
colormap redblue
colorbar
title(sprintf('T_R_E_%d as I_C_%d prediction match - shuf. dist.', ii,ii))
subplot(2,3,3*(ii-1)+3)
hold all
tempmat = squeeze(nanmean(TREasICmatchprctile{ii},3));
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
title({sprintf('T_R_E_%d as I_C_%d trial prediction match', ii,ii), 'percentile w.r.t. null distribution'})
end

figure; 
for ii = 1:numel(TREtrialtypes)
subplot(2,2,ii)
hold all
for ises = 1:Nsessions
plot(reshape(nanmean(TREpredmatchnull{ii}(:,:,:,ises),3),[],1), reshape(TREpredmatch{ii}(:,:,ises),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
xvec = reshape(nanmean(TREpredmatchnull{ii},3),[],1);
yvec = reshape(TREpredmatch{ii},[],1);
p=signrank(xvec(~isnan(xvec) & ~isnan(yvec)), yvec(~isnan(xvec) & ~isnan(yvec)));
figure; plot(xvec,yvec,'o');hold on; xl-xlim;plot(xl,xl,'r-')
p=signrank(reshape(nanmean(TREpredmatchnull{ii},3),[],1), reshape(TREpredmatch{ii},[],1));
title(sprintf('T_R_E_%d pred p=%.4f', ii,p))

subplot(2,2,2+ii)
hold all
for ises = 1:Nsessions
plot(reshape(nanmean(TREasICmatchnull{ii}(:,:,:,ises),3),[],1), reshape(TREasICmatch{ii}(:,:,ises),[],1), 'o')
end
xl = xlim;
plot(xl, xl, 'r-')
plot(xl, [1 1]/numel(traintrialtypes), 'k--')
xlabel('null distribution mean')
ylabel('actual match')
p = signrank(reshape(nanmean(TREasICmatchnull{ii},3),[],1), reshape(TREasICmatch{ii},[],1));
title(sprintf('T_R_E_%d as I_C_%d p=%.4f', ii,ii,p))
end

for ab = 1:2
    switch ab
        case 1
a=1;b=2;
        case 2
a=1;b=4;
    end
pTREpred=signrank(reshape(nanmean(TREpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(TREpredmatch{ii}(a,b,:),[],1));
pTREasIC=signrank(reshape(nanmean(TREasICmatchnull{ii}(a,b,:,:),3),[],1), reshape(TREasICmatch{ii}(a,b,:),[],1));
fprintf('%s & %s match TRE p=%.4f TRE as IC p=%.4f\n', visarealabels{a}, visarealabels{b}, pTREpred, pTREasIC)
end

for a=1:numel(visareas)
    for b=a+1:numel(visareas)        
pTREpred=signrank(reshape(nanmean(TREpredmatchnull{ii}(a,b,:,:),3),[],1), reshape(TREpredmatch{ii}(a,b,:),[],1));
pTREasIC=signrank(reshape(nanmean(TREasICmatchnull{ii}(a,b,:,:),3),[],1), reshape(TREasICmatch{ii}(a,b,:),[],1));
fprintf('%s & %s match N=%d TRE p=%.4f TRE as IC p=%.4f\n', visarealabels{a}, visarealabels{b}, ...
    nnz(~isnan(TREpredmatch{ii}(a,b,:))), pTREpred, pTREasIC)
    end
end
%********** V1 & RL match N=11 TRE p=0.0068 TRE as IC p=0.0186
% V1 & LM match N=8 TRE p=0.1094 TRE as IC p=0.1094
% V1 & RL match N=11 TRE p=0.0068 TRE as IC p=0.0186
% V1 & AL match N=10 TRE p=0.1055 TRE as IC p=0.1934
% V1 & PM match N=12 TRE p=0.1294 TRE as IC p=0.1514
% V1 & AM match N=11 TRE p=0.7646 TRE as IC p=0.4131
% LM & RL match N=7 TRE p=1.0000 TRE as IC p=0.3750
% LM & AL match N=6 TRE p=1.0000 TRE as IC p=0.8438
% LM & PM match N=8 TRE p=0.1484 TRE as IC p=0.4609
% LM & AM match N=7 TRE p=0.1562 TRE as IC p=0.2188
% RL & AL match N=9 TRE p=0.5703 TRE as IC p=0.2031
% RL & PM match N=11 TRE p=0.1230 TRE as IC p=0.2402
% RL & AM match N=10 TRE p=0.5566 TRE as IC p=0.5566
% AL & PM match N=10 TRE p=0.0840 TRE as IC p=0.2324
% AL & AM match N=9 TRE p=0.0195 TRE as IC p=0.5703
% PM & AM match N=11 TRE p=0.7002 TRE as IC p=0.8984

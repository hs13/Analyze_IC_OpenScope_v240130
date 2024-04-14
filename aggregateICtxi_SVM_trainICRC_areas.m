pathpp = '/Users/hyeyoung/Documents/DATA/OpenScopeData/00248_v240130/postprocessed/sub-619296/';
pathsvm = '/Users/hyeyoung/Documents/DATA/OpenScopeData/00248_v240130/SVM_trainICRC_selectareas/sub-619296/';

%%
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

svmdesc = 'trainICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
whichblock = [whichICblock '_presentations'];
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};

switch svmdesc
    case 'trainICRC'
        traintrialtypes = [106, 107, 110, 111];
        probetrialtypes = [1105, 1109];
    case 'trainREx'
        traintrialtypes = [1201, 1299];
        probetrialtypes = [106, 107, 110, 111];
end

%%
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
pathsv = [datadir 'SVM_' svmdesc '_selectareas' filesep];
pathsvm = [pathsv nwbsessions{ises} filesep];

load([pathpp 'postprocessed.mat'])

SVMall = struct();
for a = 1:numel(visareas)
    clearvars SVMtrainICRC SVMtrainREx
    whichvisarea = visareas{a};
    svmfn = [pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'];
    if ~exist(svmfn, 'file')
        SVMall.(whichvisarea) = [];
    else
        load(svmfn)
        switch svmdesc
            case 'trainICRC'
                SVMall.(whichvisarea) = SVMtrainICRC;
            case 'trainREx'
                SVMall.(whichvisarea) = SVMtrainREx;
            otherwise
                error([svmdesc ' not recognized'])
        end
    end
end

%%
trainaccuracy = NaN(numel(visareas), 1);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    if ~isempty(SVMall.(whichvisarea))
    traintrialstt = SVMall.(whichvisarea).trialorder(SVMall.(whichvisarea).spkcnt.traintrialinds);
    trainaccuracy(a) = mean(SVMall.(whichvisarea).spkcnt.train.label == traintrialstt, 'all');
    end
end
disp(trainaccuracy)

%% test trials: check that there is correlation between area decoders
% this correlation must exceed null distribution, where null distribution
% is obtained by shuffling within trial types

% perhaps train and test trial divide should have been the same between all
% areas...
Nshuf = 1000;

testaccuracy = NaN(numel(visareas), 1);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    if ~isempty(SVMall.(whichvisarea))
    testtrialstt = SVMall.(whichvisarea).trialorder(SVMall.(whichvisarea).spkcnt.testtrialinds);
    testaccuracy(a) = mean(SVMall.(whichvisarea).spkcnt.test.label == testtrialstt, 'all');
    end
end
disp(testaccuracy)

testpredmatchchance = testaccuracy * testaccuracy';

tic
testpredmatch = NaN(numel(visareas), numel(visareas));
testpredmatchnull = NaN(numel(visareas), numel(visareas), Nshuf);
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    [ttiA,testtrialorderA]=sort(SVMall.(whichvisareaA).spkcnt.testtrialinds(:));
    testlabelA = SVMall.(whichvisareaA).spkcnt.test.label(testtrialorderA);
    for b = a+1:numel(visareas)
        whichvisareaB = visareas{b};
        [ttiB,testtrialorderB]=sort(SVMall.(whichvisareaB).spkcnt.testtrialinds(:));
        if ~isequal(ttiA, ttiB)
            error('set of test trial inds do not match between areas')
        end
        testlabelB = SVMall.(whichvisareaB).spkcnt.test.label(testtrialorderB);
        testtrialsttB = SVMall.(whichvisareaB).trialorder(ttiB);
        testpredmatch(a,b) = mean(testlabelA==testlabelB);
        
        testaccB = mean(SVMall.(whichvisareaB).spkcnt.test.label==SVMall.(whichvisareaB).trialorder(SVMall.(whichvisareaB).spkcnt.testtrialinds), 'all');
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
            testpredmatchnull(a,b,ishuf) = mean(testlabelA==testlabelBshuf);
        end
    end
end
toc

figure; hold all
plot(testpredmatchchance, squeeze(nanmean(testpredmatchnull,3)), 'o')
plot(testpredmatchchance, testpredmatchchance, '-')
xlabel('chance (if independent)')
ylabel('null distribution mean')

figure; hold all
plot(squeeze(nanmean(testpredmatchnull,3)), testpredmatch, 'o');
plot(testpredmatchchance, testpredmatchchance, '-')
xlabel('null distribution mean')
ylabel('actual match')

figure; imagesc(testpredmatch); colorbar

testpredmatchprctile =  NaN(numel(visareas), numel(visareas));
for a = 1:numel(visareas)
    for b = a+1:numel(visareas)
        testpredmatchprctile(a,b) = 100*mean(testpredmatchnull(a,b,:)<testpredmatch(a,b));
    end
end

%% inference trials: 
whichvisareaA = 'VISp';
whichvisareaB = 'VISal';

ttbe = 0.5*([traintrialtypes(1) traintrialtypes] + [traintrialtypes traintrialtypes(end)]);
probehc2 = NaN(length(traintrialtypes), length(traintrialtypes), length(probetrialtypes));
for iprobe = 1:length(probetrialtypes)

probetrialsA = SVMall.(whichvisareaA).trialorder==probetrialtypes(iprobe);
probepredA = mode( SVMall.(whichvisareaA).spkcnt.all.label(probetrialsA,:),2 );

probetrialsB = SVMall.(whichvisareaB).trialorder==probetrialtypes(iprobe);
probepredB = mode( SVMall.(whichvisareaB).spkcnt.all.label(probetrialsB,:),2 );

hc = histcounts2(probepredA, probepredB, 'XBinEdges', ttbe, 'YBinEdges', ttbe, 'normalization', 'probability');
probehc2(:,:,iprobe) = hc;
end

figure; 
for iprobe = 1:length(probetrialtypes)
    hc = squeeze(probehc2(:,:,iprobe));
    subplot(1,2,iprobe)
imagesc(hc')
set(gca, 'XTick', 1:numel(traintrialtypes), 'XTickLabel', traintrialtypes, ...
    'YTick', 1:numel(traintrialtypes), 'YTickLabel', traintrialtypes)
xlabel(whichvisareaA)
ylabel(whichvisareaB)
caxis([0 0.25])
end

figure; histogram2(probepredA, probepredB, 'displaystyle', 'tile', 'XBinEdges', ttbe, 'YBinEdges', ttbe, 'normalization', 'probability');
caxis([0 0.25])
xlabel(whichvisareaA)
ylabel(whichvisareaB)

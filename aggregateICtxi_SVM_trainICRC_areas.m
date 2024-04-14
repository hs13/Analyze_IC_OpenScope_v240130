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
whichvisarea = 'VISp';
traintrialtypes = SVMall.(whichvisarea).trialorder(SVMall.(whichvisarea).spkcnt.traintrialinds);

trainaccuracy = NaN(numel(visareas), 1);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    if ~isempty(SVMall.(whichvisarea))
    trainaccuracy(a) = mean(SVMall.(whichvisarea).spkcnt.train.label == traintrialtypes, 'all');
    end
end

%% test trials: check that there is correlation between area decoders
% this correlation must exceed null distribution, where null distribution
% is obtained by shuffling within trial types
Nshuf = 1000;
testtrialtypes = SVMall.(whichvisarea).trialorder(SVMall.(whichvisarea).spkcnt.testtrialinds);

testaccuracy = NaN(numel(visareas), 1);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    testaccuracy(a) = mean(SVMall.(whichvisarea).spkcnt.test.label == testtrialtypes, 'all');
end
testpredmatchchance

testpredmatch = NaN(numel(visareas), numel(visareas));
testpredmatchnull = NaN(numel(visareas), numel(visareas), Nshuf);
for a = 1:numel(visareas)
    whichvisareaA = visareas{a};
    testlabelA = SVMall.(whichvisareaA).spkcnt.test.label;
for b = a+1:numel(visareas)
    whichvisareaB = visareas{b};
    testlabelB = SVMall.(whichvisareaB).spkcnt.test.label;    
    testpredmatch(a,b) = mean(testlabelA==testlabelB, 'all');
    for ishuf = 1:Nshuf
        for itt = 1:numel(traintrialtypes)
        end
    end
end
end


%% inference trials: 

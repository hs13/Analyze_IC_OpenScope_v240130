datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

svmdesc = 'trainICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';
Twin = 5;

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
    otherwise
        error([svmdesc ' not recognized'])
end
%%
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
pathsvm = [pathsv nwbsessions{ises} filesep];

load([pathpp 'postprocessed.mat'])

SVMall = struct();
for a = 1:numel(visareas)
    clearvars SVMtrainICRC SVMtrainREx
    whichvisarea = visareas{a};
    load([pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])    
    switch svmdesc
        case 'trainICRC'
            SVMall.(whichvisarea) = SVMtrainICRC;
        case 'trainREx'
            SVMall.(whichvisarea) = SVMtrainREx;
        otherwise
            error([svmdesc ' not recognized'])
    end
end

SVMcumpsthall = struct();
for a = 1:numel(visareas)
    clearvars SVMcumpsth
    whichvisarea = visareas{a};
    load([pathsvm, 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
    SVMcumpsthall.(whichvisarea) = SVMcumpsth;
end

whichvisarea = 'VISp';
figure; plot(SVMall.(whichvisarea).spkcnt.all.score(:), reshape( SVMcumpsthall.(whichvisarea).score(end,:,:,:),[],1), '.')

[sv,si]=max(SVMall.(whichvisarea).spkcnt.all.score,[],2);
temp = squeeze(traintrialtypes(si));
temp1 = SVMall.(whichvisarea).spkcnt.all.label;
isequal( squeeze(traintrialtypes(si)), SVMall.(whichvisarea).spkcnt.all.label)

isequal(SVMall.(whichvisarea).spkcnt.traintrialinds, SVMcumpsthall.(whichvisarea).traintrialinds)
isequal(SVMall.(whichvisarea).spkcnt.testtrialinds, SVMcumpsthall.(whichvisarea).testtrialinds)

isplit = 1;
testtrialinds = SVMcumpsthall.(whichvisarea).testtrialinds(:,isplit);
ii=1;
itrial = testtrialinds(ii);
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
tempsvmcumpsth = squeeze( SVMcumpsthall.(whichvisarea).score(:,itrial,:,isplit) );
subplot(2,3,a)
plot(SVMcumpsthall.(whichvisarea).psthbinTends, tempsvmcumpsth)
title(whichvisarea)
end

% compare score between correct vs incorrect trials
% incorrect trials do have lower scores, but the distribution is not
% exactly separate
testtrialtypes = SVMcumpsthall.(whichvisarea).trialorder(SVMcumpsthall.(whichvisarea).testtrialinds);
Nsplits = size(testtrialtypes,2);
testlabel = NaN(size(testtrialtypes));
testscore = NaN(size(testtrialtypes));
for isplit = 1:Nsplits
testtrialinds = SVMcumpsthall.(whichvisarea).testtrialinds(:,isplit);
testlabel(:,isplit) = squeeze( SVMcumpsthall.(whichvisarea).label(end,testtrialinds,isplit) );
tempscore = squeeze( SVMcumpsthall.(whichvisarea).score(end,testtrialinds,:,isplit) );
[sv,si]=max(tempscore,[],2);
testscore(:,isplit) = squeeze(sv);
end
figure;
hold all
histogram(testscore(testlabel==testtrialtypes))
histogram(testscore(testlabel~=testtrialtypes))


cumpsthtl = SVMcumpsthall.(whichvisarea).psthbinTends;
Ntesttrials = size(SVMcumpsthall.(whichvisarea).testtrialinds,1);
Nsplits = size(SVMcumpsthall.(whichvisarea).testtrialinds,2);
testscorecumpsth = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
for isplit = 1:Nsplits
testtrialinds = SVMcumpsthall.(whichvisarea).testtrialinds(:,isplit);

testscorecumpsth(:,:,isplit)
end

% IC test trials
% first, focus on V1 and LM, and on trials where both areas' decoders had correct predictions


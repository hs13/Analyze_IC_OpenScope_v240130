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

cumpsthtl = SVMcumpsthall.VISp.psthbinTends;
[ttindsordered,~]=sort(SVMcumpsthall.VISp.testtrialinds(:));
testtrialstt = SVMcumpsthall.VISp.trialorder(ttindsordered);
testlabelfinal = struct();
testscorecumpsth = struct();
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    finallabel = squeeze(SVMcumpsthall.(whichvisarea).label(end,:,:));
    
    Ntesttrials = size(SVMcumpsthall.(whichvisarea).testtrialinds,1);
    Nsplits = size(SVMcumpsthall.(whichvisarea).testtrialinds,2);
    origtestlabel = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
    origtestscore = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
    for isplit = 1:Nsplits
        testtrialinds = SVMcumpsthall.(whichvisarea).testtrialinds(:,isplit);
        origtestlabel(:,:,isplit) = SVMcumpsthall.(whichvisarea).label(:,testtrialinds,isplit);
        for itt = 1:numel(traintrialtypes)
            temptrials = finallabel(testtrialinds,isplit)==traintrialtypes(itt);
            origtestscore(:,temptrials,isplit) = squeeze(SVMcumpsthall.(whichvisarea).score(:,testtrialinds(temptrials),itt,isplit));
        end
    end
    origtestlabel = reshape(origtestlabel, length(cumpsthtl), Ntesttrials*Nsplits);
    origtestscore = reshape(origtestscore, length(cumpsthtl), Ntesttrials*Nsplits);
    
    finaltestlabel = origtestlabel(end,:);
    
    [tti,testtrialord]=sort(SVMcumpsthall.(whichvisarea).testtrialinds(:));
    if ~isequal(ttindsordered, tti)
        error('test trial set different for every area?')
    end
    testlabelfinal.(whichvisarea) = finaltestlabel(testtrialord);
    testscorecumpsth.(whichvisarea) = origtestscore(:,testtrialord);    
end

normtestscorecumpsth = struct();
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempts = testscorecumpsth.(whichvisarea);
    normtestscorecumpsth.(whichvisarea) = (tempts-tempts(1,:))./(tempts(end,:)-tempts(1,:));
end

% IC test trials
% first, focus on V1 and LM, and on trials where both areas' decoders had correct predictions
itt =1;
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    trialsoi = testtrialstt==traintrialtypes(itt);
    subplot(1,numel(visareas),a)
    imagesc(cumpsthtl,1:nnz(trialsoi), testscorecumpsth.(whichvisarea)(:,trialsoi)')
    colorbar
end

whichvisareaA = 'VISp';
whichvisareaB = 'VISal';
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    testlabelfinal.(whichvisareaA)==traintrialtypes(itt) & testlabelfinal.(whichvisareaB)==traintrialtypes(itt) );
figure
for ii = 1:12
    subplot(3,4,ii)
    hold all
    plot(cumpsthtl, normtestscorecumpsth.(whichvisareaA)(:,trialsoind(ii)))
    plot(cumpsthtl, normtestscorecumpsth.(whichvisareaB)(:,trialsoind(ii)))
end

% find the first timepoint at which normtestscorecumpsth crosses 
% 0.25, 0.5, 0.75 (called T25, T50, T75 respsectively)
% divide into 4 trial types
% 1. areaA faster than areaB if all three timeponts (T25, T50, T75) are earlier for A than B
% 2. areaB faster than areaA if all three timeponts (T25, T50, T75) are earlier for A than B
% 3. areaA ramping starts earlier and finishes later than area B
% 4. areaB ramping starts earlier and finishes later than area A
% 5. crisscrossing if T25 and T75 go in one direction and T50 goes in the opposite direction

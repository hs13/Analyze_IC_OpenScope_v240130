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

SVMcumpsthall = struct();
for a = 1:numel(visareas)
    clearvars SVMcumpsth
    whichvisarea = visareas{a};
    load([pathsvm, 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
    SVMcumpsthall.(whichvisarea) = SVMcumpsth;
end

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

milestonesvec = [0.25; 0.50; 0.75];
Tmilestones = struct();
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempnorm = normtestscorecumpsth.(whichvisarea);
    Nalltest = size(tempnorm,2);
    Tmilestones.(whichvisarea) = NaN(length(milestonesvec), Nalltest);
    for t = 1:length(milestonesvec)
        % first timepoint to cross milestone
        [~,mi]=max(tempnorm>=milestonesvec(t),[],1);
        % v50 = tempnorm(sub2ind(size(tempnorm), mi, 1:size(tempnorm,2)));
        % v49 = tempnorm(sub2ind(size(tempnorm), mi-1, 1:size(tempnorm,2)));
        % if ~all(v50>=0.5 & v49<0.5)
        %     error('check mi')
        % end
        Tmilestones.(whichvisarea)(t,:) = cumpsthtl(mi);
    end
    if ~isequal(sort(Tmilestones.(whichvisarea),1), Tmilestones.(whichvisarea))
        error('Tmilestones does not follow order -- check algorithm')
    end
end


%% check correspondence between SVMall and final timepoint of SVMcumpsthall
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
whichvisarea = 'VISp';
figure; plot(SVMall.(whichvisarea).spkcnt.all.score(:), reshape( SVMcumpsthall.(whichvisarea).score(end,:,:,:),[],1), '.')

[sv,si]=max(SVMall.(whichvisarea).spkcnt.all.score,[],2);
temp = squeeze(traintrialtypes(si));
temp1 = SVMall.(whichvisarea).spkcnt.all.label;
isequal( squeeze(traintrialtypes(si)), SVMall.(whichvisarea).spkcnt.all.label)

isequal(SVMall.(whichvisarea).spkcnt.traintrialinds, SVMcumpsthall.(whichvisarea).traintrialinds)
isequal(SVMall.(whichvisarea).spkcnt.testtrialinds, SVMcumpsthall.(whichvisarea).testtrialinds)

%% compare score between correct vs incorrect trials
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

%% IC test trials
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
% example trials
figure
for ii = 1:12
    subplot(3,4,ii)
    hold all
    plot(cumpsthtl, normtestscorecumpsth.(whichvisareaA)(:,trialsoind(ii)))
    plot(cumpsthtl, normtestscorecumpsth.(whichvisareaB)(:,trialsoind(ii)))
end

% average across trials
figure
for itt = 1:numel(traintrialtypes)
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    testlabelfinal.(whichvisareaA)==traintrialtypes(itt) & testlabelfinal.(whichvisareaB)==traintrialtypes(itt) );
subplot(2,2,itt)
hold all
shadedErrorBar(cumpsthtl, mean(normtestscorecumpsth.(whichvisareaA)(:,trialsoind),2), ...
    std(normtestscorecumpsth.(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
shadedErrorBar(cumpsthtl, mean(normtestscorecumpsth.(whichvisareaB)(:,trialsoind),2), ...
    std(normtestscorecumpsth.(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
title(traintrialtypes(itt))
end

% find the first timepoint at which normtestscorecumpsth crosses 
% 0.25, 0.5, 0.75 (called T25, T50, T75 respsectively)
% divide into 4 trial types
% 1. areaA ramping starts earlier and finishes later than area B
% 2. areaB ramping starts earlier and finishes later than area A
% 3. areaA faster than areaB if all three timeponts (T25, T50, T75) are earlier for A than B
% 4. areaB faster than areaA if all three timeponts (T25, T50, T75) are earlier for A than B
% 5. simultaneous: all three timeponts are identical
% 6. crisscrossing if T25 and T75 go in one direction and T50 goes in the opposite direction
whichvisareaA = 'VISp';
whichvisareaB = 'VISal';
Ntestall = size(normtestscorecumpsth.(whichvisareaA),2);
comparedynamics = zeros(1,Ntestall);
for dyn = 1:6
switch dyn
    case 1
        trialsinclass = Tmilestones.(whichvisareaA)(1,:)<=Tmilestones.(whichvisareaB)(1,:) ...
            & Tmilestones.(whichvisareaA)(3,:)>=Tmilestones.(whichvisareaB)(3,:);
    case 2
        trialsinclass = Tmilestones.(whichvisareaA)(1,:)>=Tmilestones.(whichvisareaB)(1,:) ...
            & Tmilestones.(whichvisareaA)(3,:)<=Tmilestones.(whichvisareaB)(3,:);
    case 3
trialsinclass = all(Tmilestones.(whichvisareaA)<=Tmilestones.(whichvisareaB),1);
    case 4
trialsinclass = all(Tmilestones.(whichvisareaA)>=Tmilestones.(whichvisareaB),1);
    case 5
trialsinclass = all(Tmilestones.(whichvisareaA)==Tmilestones.(whichvisareaB),1);
    case 6
        temp = Tmilestones.(whichvisareaA)>Tmilestones.(whichvisareaB);
        trialsinclass = temp(1,:)==~temp(2,:) & temp(3,:)==~temp(2,:);
end
comparedynamics(trialsinclass) = dyn;
end

itt = 4;
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    testlabelfinal.(whichvisareaA)==traintrialtypes(itt) & testlabelfinal.(whichvisareaB)==traintrialtypes(itt) );
[v,c]=uniquecnt(comparedynamics(trialsoind));
disp([v',c'])

% compare ramp time: longer means slower
rampA = Tmilestones.(whichvisareaA)(3,:)-Tmilestones.(whichvisareaA)(1,:);
rampB = Tmilestones.(whichvisareaB)(3,:)-Tmilestones.(whichvisareaB)(1,:);
figure; hold all
plot(rampA(trialsoind), rampB(trialsoind), 'o')
xl = xlim;
plot(xl,xl, '-')
p = signrank(rampA(trialsoind), rampB(trialsoind));
fprintf('ramp delay (ms) %d trials %s vs %s p=%.4f\n', traintrialtypes(itt), ...
    whichvisareaA, whichvisareaB, p)
fprintf('median: %.2f vs %.2f, mean: %.2f vs %.2f\n', ...
    median(rampA(trialsoind)), median(rampB(trialsoind)), ...
    mean(rampA(trialsoind)), mean(rampB(trialsoind)))


T50A = Tmilestones.(whichvisareaA)(2,:);
T50B = Tmilestones.(whichvisareaB)(2,:);
figure; hold all
plot(T50A(trialsoind), T50B(trialsoind), 'o')
xl = xlim;
plot(xl,xl, '-')
p = signrank(T50A(trialsoind), T50B(trialsoind));
fprintf('T50 %d trials %s vs %s p=%.4f\n', traintrialtypes(itt), whichvisareaA, whichvisareaB, p)
fprintf('median: %.2f vs %.2f, mean: %.2f vs %.2f\n', ...
    median(T50A(trialsoind)), median(T50B(trialsoind)), ...
    mean(T50A(trialsoind)), mean(T50B(trialsoind)))


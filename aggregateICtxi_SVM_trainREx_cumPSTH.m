datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

Twin = 5;
neuopt = 'RS';
svmdesc = 'trainREx';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

pathsv = [datadir 'SVM_' svmdesc '_selectareas' filesep];
whichICblock = 'ICwcfg1';
whichblock = [whichICblock '_presentations'];
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};

switch svmdesc
    case 'trainICRC'
        traintrialtypes = [106, 107, 110, 111];
    case 'trainREx'
        traintrialtypes = [1201, 1299];
        ICtrialtypes = [106, 111];
    otherwise
        error([svmdesc ' not recognized'])
end

SVMtrainRExcumpsthagg = struct();
RExtestlabelfinalagg = struct();
RExtestscorecumpsthagg = struct();
RExtestnormscorecumpsthagg = struct();
RExtestTmilestonesagg = struct();

RExinfvalidmodeagg = struct(); % valid if there is consensus among cross-validations
% i.e., valid if mode prediction across Nsplits is not tied
RExinfmodefreqagg = struct();
RExinfvalidmaxscoreagg = struct();
RExinfmaxscoresplitagg = struct();
RExinflabelfinalagg = struct();
RExinfscorecumpsthagg = struct();
RExinfnormscorecumpsthagg = struct();
RExinfTmilestonesagg = struct();

for ises = 1:Nsessions
    tic
    pathsvm = [pathsv nwbsessions{ises} filesep];

    for a = 1:numel(visareas)
        clearvars SVMcumpsth
        whichvisarea = visareas{a};
        svmcumpsthfn = [pathsvm, 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'];
        if exist(svmcumpsthfn, 'file')
            load(svmcumpsthfn)
            SVMtrainRExcumpsthagg(ises).(whichvisarea) = SVMcumpsth;
        else
            SVMtrainRExcumpsthagg(ises).(whichvisarea) = [];
        end
    end

    cumpsthtl = SVMtrainRExcumpsthagg(ises).VISp.psthbinTends;
    [ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
    testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        finallabel = squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).label(end,:,:));

        Ntesttrials = size(SVMtrainRExcumpsthagg(ises).(whichvisarea).testtrialinds,1);
        Nsplits = size(SVMtrainRExcumpsthagg(ises).(whichvisarea).testtrialinds,2);
        origtestlabel = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
        origtestscore = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
        for isplit = 1:Nsplits
            testtrialinds = SVMtrainRExcumpsthagg(ises).(whichvisarea).testtrialinds(:,isplit);
            origtestlabel(:,:,isplit) = SVMtrainRExcumpsthagg(ises).(whichvisarea).label(:,testtrialinds,isplit);
            for itt = 1:numel(traintrialtypes)
                temptrials = finallabel(testtrialinds,isplit)==traintrialtypes(itt);
                origtestscore(:,temptrials,isplit) = squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).score(:,testtrialinds(temptrials),itt,isplit));
            end
        end
        origtestlabel = reshape(origtestlabel, length(cumpsthtl), Ntesttrials*Nsplits);
        origtestscore = reshape(origtestscore, length(cumpsthtl), Ntesttrials*Nsplits);

        finaltestlabel = origtestlabel(end,:);

        [tti,testtrialord]=sort(SVMtrainRExcumpsthagg(ises).(whichvisarea).testtrialinds(:));
        if ~isequal(ttindsordered, tti)
            error('test trial set different for every area?')
        end
        RExtestlabelfinalagg(ises).(whichvisarea) = finaltestlabel(testtrialord);
        RExtestscorecumpsthagg(ises).(whichvisarea) = origtestscore(:,testtrialord);
    end

    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempts = RExtestscorecumpsthagg(ises).(whichvisarea);
        trialsvalramp = tempts(end,:)-tempts(1,:) >0;
        RExtestnormscorecumpsth = (tempts-tempts(1,:))./(tempts(end,:)-tempts(1,:));
        RExtestnormscorecumpsth(:,~trialsvalramp)=NaN;

        RExtestnormscorecumpsthagg(ises).(whichvisarea) = RExtestnormscorecumpsth;
    end

    milestonesvec = [0.25; 0.50; 0.75];
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempnorm = RExtestnormscorecumpsthagg(ises).(whichvisarea);
        Nalltest = size(tempnorm,2);
        RExtestTmilestonesagg(ises).(whichvisarea) = NaN(length(milestonesvec), Nalltest);
        for t = 1:length(milestonesvec)
            % first timepoint to cross milestone
            [~,mi]=max(tempnorm>=milestonesvec(t),[],1);
            % v50 = tempnorm(sub2ind(size(tempnorm), mi, 1:size(tempnorm,2)));
            % v49 = tempnorm(sub2ind(size(tempnorm), mi-1, 1:size(tempnorm,2)));
            % if ~all(v50>=0.5 & v49<0.5)
            %     error('check mi')
            % end
            RExtestTmilestonesagg(ises).(whichvisarea)(t,:) = cumpsthtl(mi);
        end
        if ~isequal(sort(RExtestTmilestonesagg(ises).(whichvisarea),1), RExtestTmilestonesagg(ises).(whichvisarea))
            error('RExtestTmilestones does not follow order -- check algorithm')
        end
    end

    % for inference, chose mode prediction as the final label
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        finallabel = squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).label(end,:,:));

        ICtrialinds = find( ismember(SVMtrainRExcumpsthagg(ises).(whichvisarea).trialorder, ICtrialtypes) );
        ICtrialstt = SVMtrainRExcumpsthagg(ises).(whichvisarea).trialorder(ICtrialinds);
        NICtrials = numel(ICtrialinds);
        originflabel = SVMtrainRExcumpsthagg(ises).(whichvisarea).label(:,ICtrialinds,:);
        originfscore = SVMtrainRExcumpsthagg(ises).(whichvisarea).score(:,ICtrialinds,:,:);

        % isequal( squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).label(end,:,:)), SVMtrainRExagg.(whichICblock)(ises).(whichvisarea).spkcnt.all.label)
        [M,F,C]=mode(originflabel,3);
        Cnumel = cellfun(@numel, C);

        finalscore = squeeze(originfscore(end,:,:,:));
        finalmodescore = NaN(NICtrials,Nsplits);
        originfmodecore = NaN(length(cumpsthtl), NICtrials, Nsplits);
        for typi = 1:numel(traintrialtypes)
            temptrials = M(end,:)==traintrialtypes(typi);
            finalmodescore(temptrials,:) = squeeze(finalscore(temptrials,typi,:));
            originfmodecore(:,temptrials,:) = squeeze(originfscore(:,temptrials,typi,:));
        end
        [mv,maxscoresplit]=max(finalmodescore,[],2);
        mmv=max( reshape(finalscore, NICtrials,numel(traintrialtypes)*Nsplits),[],2);

        RExinfscorecumpsth = NaN(length(cumpsthtl), NICtrials);
        for isplit = 1:Nsplits
            temptrials = maxscoresplit==isplit;
            RExinfscorecumpsth(:,temptrials) = squeeze( originfmodecore(:,temptrials,isplit) );
        end

        trialsvalramp = RExinfscorecumpsth(end,:)-RExinfscorecumpsth(1,:)>0;
        RExinfnormscorecumpsth = (RExinfscorecumpsth-RExinfscorecumpsth(1,:))./(RExinfscorecumpsth(end,:)-RExinfscorecumpsth(1,:));
        RExinfnormscorecumpsth(:,~trialsvalramp)=NaN;

        RExinflabelfinalagg(ises).(whichvisarea) = M(end,:);
        % RExinflabelcumpsthagg(ises).(whichvisarea) = M;
        RExinfmodefreqagg(ises).(whichvisarea) = F;
        RExinfvalidmodeagg(ises).(whichvisarea) = Cnumel==1;
        % trial is valid if there is consensus among cross-validations
        % i.e., trial is valid if mode prediction across Nsplits is not tied
        RExinfvalidmaxscoreagg(ises).(whichvisarea) = mv==mmv;
        RExinfmaxscoresplitagg(ises).(whichvisarea) = maxscoresplit;
        RExinfscorecumpsthagg(ises).(whichvisarea) = RExinfscorecumpsth;
        RExinfnormscorecumpsthagg(ises).(whichvisarea) = RExinfnormscorecumpsth;
    end

    milestonesvec = [0.25; 0.50; 0.75];
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempnorm = RExinfnormscorecumpsthagg(ises).(whichvisarea);
        NICtrials = size(tempnorm,2);
        RExinfTmilestones = NaN(length(milestonesvec), NICtrials);
        for t = 1:length(milestonesvec)
            % first timepoint to cross milestone
            [~,mi]=max(tempnorm>=milestonesvec(t),[],1);
            % v50 = tempnorm(sub2ind(size(tempnorm), mi, 1:size(tempnorm,2)));
            % v49 = tempnorm(sub2ind(size(tempnorm), mi-1, 1:size(tempnorm,2)));
            % if ~all(v50>=0.5 & v49<0.5)
            %     error('check mi')
            % end
            RExinfTmilestones(t,:) = cumpsthtl(mi);
        end
        if ~isequal(sort(RExinfTmilestones,1), RExinfTmilestones)
            error('RExinfTmilestones does not follow order -- check algorithm')
        end

        RExinfTmilestonesagg(ises).(whichvisarea) = RExinfTmilestones;
    end

    toc
end

save([pathsv 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc '_' neuopt 'agg.mat'], 'SVMtrainRExcumpsthagg', ...
    'RExtestlabelfinalagg', 'RExtestscorecumpsthagg', 'RExtestnormscorecumpsthagg', 'RExtestTmilestonesagg', ...
    'RExinfvalidmodeagg', 'RExinfmodefreqagg', 'RExinfvalidmaxscoreagg', 'RExinfmaxscoresplitagg', ...
    'RExinflabelfinalagg', 'RExinfscorecumpsthagg', 'RExinfnormscorecumpsthagg', 'RExinfTmilestonesagg', '-v7.3')

%%
load('G:\My Drive\RESEARCH\ICexpts_revision23\openscope_HR_SVMtrainREx_zscore_agg.mat', 'Nneuronsperarea')
discardbelowNneurons = 50;

%% valid ramp trials: define based on spearman correlation between timepoints and scores
% set an arbitrary threshold at 0.5 to determine validity.
% alternatively, could use a more inclusive threshold value of 0
threshspear = 0.5;
RExtestrhorampxtimeagg = struct();
RExinfrhorampxtimeagg = struct();
for ises = 1:Nsessions
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempts = RExtestscorecumpsthagg(ises).(whichvisarea);
        RExtestrhorampxtimeagg(ises).(whichvisarea) = corr(tempts, (1:length(cumpsthtl))', 'type', 'Spearman');

        tempts = RExinfscorecumpsthagg(ises).(whichvisarea);
        RExinfrhorampxtimeagg(ises).(whichvisarea) = corr(tempts, (1:length(cumpsthtl))', 'type', 'Spearman');
    end
end

rhorampxtime = cat(2,RExtestrhorampxtimeagg.(whichvisarea));
figure; histogram(rhorampxtime(:), 'binwidth', 0.01)

tempts = RExtestscorecumpsthagg(ises).(whichvisarea);
rhoampxtimevec = RExtestrhorampxtimeagg(ises).(whichvisarea);
figure
hold all
plot(cumpsthtl, tempts(:,rhoampxtimevec>0.5), 'b-')
plot(cumpsthtl, tempts(:,rhoampxtimevec<=0.5), 'r-')

tempts = RExinfscorecumpsthagg(ises).(whichvisarea);
rhoampxtimevec = RExinfrhorampxtimeagg(ises).(whichvisarea);
figure
hold all
plot(cumpsthtl, tempts(:,rhoampxtimevec>0.5), 'b-')
plot(cumpsthtl, tempts(:,rhoampxtimevec<=0.5), 'r-')

%% check correspondence between SVMall and final timepoint of SVMcumpsthall
SVMcumpsthall = SVMtrainRExcumpsthagg(ises);
pathsvm = [pathsv nwbsessions{ises} filesep];

SVMall = struct();
whichvisarea = 'VISp';
% for a = 1:numel(visareas)
%     whichvisarea = visareas{a};
clearvars SVMtrainICRC SVMtrainREx
load([pathsvm, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
switch svmdesc
    case 'trainICRC'
        SVMall.(whichvisarea) = SVMtrainICRC;
    case 'trainREx'
        SVMall.(whichvisarea) = SVMtrainREx;
    otherwise
        error([svmdesc ' not recognized'])
end
% end

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
SVMcumpsthall = SVMtrainRExcumpsthagg(ises);
whichvisarea = 'VISp';

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

%% XRE test trials
% first, focus on V1 and LM, and on trials where both areas' decoders had correct predictions
ises = 11;
[ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);

figure
for itt = 1:numel(traintrialtypes)
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        trialsoi = testtrialstt==traintrialtypes(itt);
        subplot(numel(traintrialtypes),numel(visareas),(itt-1)*numel(visareas)+a)
        imagesc(cumpsthtl,1:nnz(trialsoi), RExtestscorecumpsthagg(ises).(whichvisarea)(:,trialsoi)')
        colorbar
        title(sprintf('X_R_E_%d test %s', itt, visareas{a}))
    end
end
ises = 11;
figure
for itt = 1:numel(traintrialtypes)
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        trialsoi = testtrialstt==traintrialtypes(itt);
        subplot(numel(traintrialtypes),numel(visareas),(itt-1)*numel(visareas)+a)
        imagesc(cumpsthtl,1:nnz(trialsoi), RExtestnormscorecumpsthagg(ises).(whichvisarea)(:,trialsoi)')
        colorbar
        caxis([0 1])
        title(sprintf('X_R_E_%d test %s', itt, visareas{a}))
    end
end

% check that reordering went well
Ntt2p = 4;
[ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);

tt2p = randperm(length(testtrialstt),Ntt2p);
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for t = 1:Ntt2p
        itrial = tt2p(t);
        subplot(Ntt2p, numel(visareas), (t-1)*numel(visareas)+a)
        hold all
        plot(cumpsthtl, RExtestscorecumpsthagg(ises).(whichvisarea)(:,itrial), 'k--', 'linewidth', 1)
        [r,c]=find(SVMtrainRExcumpsthagg(ises).(whichvisarea).testtrialinds==ttindsordered(itrial));
        tempcumpsth = squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).score(:,ttindsordered(itrial),:,c));
        plot(cumpsthtl, tempcumpsth)
        title(sprintf('%s test trial #%d', whichvisarea, itrial))
    end
end

whichvisareaA = 'VISp';
whichvisareaB = 'VISal';
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );

% example trials
figure
for ii = 1:12
    subplot(3,4,ii)
    hold all
    plot(cumpsthtl, RExtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind(ii)))
    plot(cumpsthtl, RExtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind(ii)))
end

% average across trials
whichvisareaA = 'VISp';
whichvisareaB = 'VISl';
figure
for itt = 1:numel(traintrialtypes)
    trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
        RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
    subplot(1,2,itt)
    hold all
    shadedErrorBar(cumpsthtl, nanmean(RExtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
        nanstd(RExtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
    shadedErrorBar(cumpsthtl, nanmean(RExtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
        nanstd(RExtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
    title(traintrialtypes(itt))
end

%% test mean ramppsth
ab = 2;
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

figure
for ises = 1:Nsessions
    a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
    if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
        continue
    end
    [ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
    testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);
    validramptrials = RExtestrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExtestrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

    for itt = 1:numel(traintrialtypes)+1
        if itt<=numel(traintrialtypes)
            trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
                RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
        else
            trialsoind = find( validramptrials' & ...
                RExtestlabelfinalagg(ises).(whichvisareaA)==testtrialstt & RExtestlabelfinalagg(ises).(whichvisareaB)==testtrialstt );
        end
        subplot( (numel(traintrialtypes)+1)*2,Nsessions, (itt-1)*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(RExtestscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(RExtestscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(RExtestscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(RExtestscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('test score')
        if itt<=numel(traintrialtypes)
            title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))
        else
            title(sprintf('Ses%d Pool XRE', ises ))
        end

        subplot( (numel(traintrialtypes)+1)*2,Nsessions, (itt+numel(traintrialtypes))*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(RExtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(RExtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(RExtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(RExtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('norm test score')
        if itt<=numel(traintrialtypes)
            title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))
        else
            title(sprintf('Ses%d Pool XRE', ises ))
        end
    end
end

% %% inf mean ramppsth
% ab = 2;
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

figure
for ises = 1:Nsessions
    a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
    if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
        continue
    end
    ICtrialinds = find( ismember(SVMtrainRExcumpsthagg(ises).VISp.trialorder, ICtrialtypes) );
    ICtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ICtrialinds);
    validramptrials = RExinfrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExinfrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

    ICtrialsttconvert = ICtrialstt;
    ICtrialsttconvert(ICtrialstt==106) = 1201;
    ICtrialsttconvert(ICtrialstt==111) = 1299;
    if ~isequal( unique(ICtrialsttconvert), traintrialtypes )
        error('check ICtrialsttconvert')
    end

    for itt = 1:numel(traintrialtypes)+1
        if itt<=numel(traintrialtypes)
            trialsoind = find( validramptrials' & ICtrialstt==ICtrialtypes(itt) & ...
                RExinflabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExinflabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
        else
            trialsoind = find( validramptrials' & ...
                RExinflabelfinalagg(ises).(whichvisareaA)==ICtrialsttconvert & RExinflabelfinalagg(ises).(whichvisareaB)==ICtrialsttconvert );
        end
        subplot( (numel(traintrialtypes)+1)*2,Nsessions, (itt-1)*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(RExinfscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(RExinfscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(RExinfscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(RExinfscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('inf score')
        if itt<=numel(traintrialtypes)
            title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))
        else
            title(sprintf('Ses%d Pool XRE', ises ))
        end

        subplot( (numel(traintrialtypes)+1)*2,Nsessions, (itt+numel(traintrialtypes))*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(RExinfnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(RExinfnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(RExinfnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(RExinfnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('norm inf score')
        if itt<=numel(traintrialtypes)
            title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))
        else
            title(sprintf('Ses%d Pool XRE', ises ))
        end
    end
end

%% COMPARE DYNAMICS OF TWO PAREAS TRIAL-BY-TRIAL: test trials, when both areas had the *correct* prediction
% find the first timepoint at which normtestscorecumpsth crosses
% 0.25, 0.5, 0.75 (called T25, T50, T75 respsectively)
% divide into 6 trial types
% 1. areaA ramping starts earlier and finishes later than area B
% 2. areaB ramping starts earlier and finishes later than area A
% 3. areaA faster than areaB if all three timeponts (T25, T50, T75) are earlier for A than B
% 4. areaB faster than areaA if all three timeponts (T25, T50, T75) are earlier for A than B
% 5. simultaneous: all three timeponts are identical
% 6. crisscrossing if T25 and T75 go in one direction and T50 goes in the opposite direction
dynamicslabels = 0:6;
RExtestcomparedynamicsprob = struct();
RExtestcomparedynamicscnt = struct();
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
    RExtestcomparedynamicsprob.(ABfield) = zeros(numel(traintrialtypes)+1, length(dynamicslabels), Nsessions);
    RExtestcomparedynamicscnt.(ABfield) = zeros(numel(traintrialtypes)+1, length(dynamicslabels), Nsessions);
    for ises = 1:Nsessions
        a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
        if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
            RExtestcomparedynamicsprob.(ABfield)(:,:,ises) = NaN;
            RExtestcomparedynamicscnt.(ABfield)(:,:,ises) = NaN;
            continue
        end
        [ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
        testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);
        validramptrials = RExtestrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExtestrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

        Ntestall = size(RExtestnormscorecumpsthagg(ises).(whichvisareaA),2);
        comparedynamics = NaN(1,Ntestall);
        for dyn = 1:6
            switch dyn
                case 1
                    trialsinclass = RExtestTmilestonesagg(ises).(whichvisareaA)(1,:)<=RExtestTmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & RExtestTmilestonesagg(ises).(whichvisareaA)(3,:)>=RExtestTmilestonesagg(ises).(whichvisareaB)(3,:);
                case 2
                    trialsinclass = RExtestTmilestonesagg(ises).(whichvisareaA)(1,:)>=RExtestTmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & RExtestTmilestonesagg(ises).(whichvisareaA)(3,:)<=RExtestTmilestonesagg(ises).(whichvisareaB)(3,:);
                case 3
                    trialsinclass = all(RExtestTmilestonesagg(ises).(whichvisareaA)<=RExtestTmilestonesagg(ises).(whichvisareaB),1);
                case 4
                    trialsinclass = all(RExtestTmilestonesagg(ises).(whichvisareaA)>=RExtestTmilestonesagg(ises).(whichvisareaB),1);
                case 5
                    trialsinclass = all(RExtestTmilestonesagg(ises).(whichvisareaA)==RExtestTmilestonesagg(ises).(whichvisareaB),1);
                case 6
                    temp = RExtestTmilestonesagg(ises).(whichvisareaA)>RExtestTmilestonesagg(ises).(whichvisareaB);
                    trialsinclass = temp(1,:)==~temp(2,:) & temp(3,:)==~temp(2,:);
            end
            comparedynamics(trialsinclass) = dyn;
        end

        for itt = 1:numel(traintrialtypes)+1
            if itt<=numel(traintrialtypes)
                trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
                    RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
            else
                trialsoind = find( validramptrials' & ...
                    RExtestlabelfinalagg(ises).(whichvisareaA)==testtrialstt & RExtestlabelfinalagg(ises).(whichvisareaB)==testtrialstt );
            end
            [v,c]=uniquecnt(comparedynamics(trialsoind));
            % disp([v',c'])
            RExtestcomparedynamicsprob.(ABfield)(itt, ismember(dynamicslabels,v), ises) = c/numel(trialsoind);
            RExtestcomparedynamicscnt.(ABfield)(itt, ismember(dynamicslabels,v), ises) = c;
        end
    end
end

figure
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
    for itt = 1:numel(traintrialtypes)+1
        subplot(2,numel(traintrialtypes)+1,(ab-1)*(numel(traintrialtypes)+1)+itt)
        imagesc(squeeze(RExtestcomparedynamicsprob.(ABfield)(itt,:,:))')
        set(gca,'Xtick',1:length(dynamicslabels), 'XtickLabel', dynamicslabels)
        if itt<=numel(traintrialtypes)
            title(sprintf('%s vs %s Trial %d', whichvisareaA, whichvisareaB, traintrialtypes(itt) ))
        else
            title(sprintf('%s vs %s Pool XRE Test Trials', whichvisareaA, whichvisareaB ))
        end
        colorbar
    end
end
colormap redblue

%{
threshspear == 0.5
rampAB Trial 1201
VISp vs VISl N=8 z-score 0.01 (-0.20~0.36) (p=0.8438, pright=0.4219, pleft=0.6289)
VISp vs VISal N=9 z-score -0.05 (-0.21~0.17) (p=1.0000, pright=0.5449, pleft=0.5000)
rampAB Trial 1299
VISp vs VISl N=7 z-score -0.24 (-0.50~-0.05) (p=0.2969, pright=0.8906, pleft=0.1484)
VISp vs VISal N=9 z-score -0.05 (-0.25~0.08) (p=0.4961, pright=0.7871, pleft=0.2480)
rampAB Pool XRE Test Trials
VISp vs VISl N=8 z-score -0.07 (-0.32~-0.00) (p=0.1484, pright=0.9453, pleft=0.0742)
VISp vs VISal N=10 z-score -0.12 (-0.19~-0.03) (p=0.0273, pright=0.9902, pleft=0.0137)
T50AB Trial 1201
VISp vs VISl N=8 z-score 0.15 (-0.29~0.23) (p=0.9453, pright=0.5781, pleft=0.4727)
VISp vs VISal N=9 z-score -0.14 (-0.24~0.12) (p=0.4961, pright=0.7871, pleft=0.2480)
T50AB Trial 1299
VISp vs VISl N=7 z-score -0.32 (-0.51~0.03) (p=0.3750, pright=0.8516, pleft=0.1875)
VISp vs VISal N=9 z-score 0.02 (-0.11~0.36) (p=0.8203, pright=0.4102, pleft=0.6328)
T50AB Pool XRE Test Trials
VISp vs VISl N=8 z-score -0.04 (-0.27~0.04) (p=0.5469, pright=0.7695, pleft=0.2734)
VISp vs VISal N=10 z-score -0.08 (-0.26~0.05) (p=0.3223, pright=0.8623, pleft=0.1611)
%}

%% compare ramp time: longer means slower
% calculate AUROC for rampA vs rampB for each session, then see if the
% distribution is significantly different from 0.5 across sessions
comparisonmetrics = {'Ntrials', 'zscore', 'Zwsr', 'TSwsr', 'CohenD', 'Zmww', 'TSmww', 'AUC'};
comparisonvectors = {'rampAB', 'T50AB'};
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};
    tempAB = struct();
    tempvecdiffAB = struct();
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
        tic
        for imet = 1:numel(comparisonmetrics)
            tempAB.(ABfield).(comparisonmetrics{imet}) = NaN(numel(traintrialtypes)+1, Nsessions);
        end
        tempvecdiffAB.(ABfield) = cell(numel(traintrialtypes)+1, Nsessions);
        for ises = 1:Nsessions
            a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
            if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
                continue
            end
            [ttindsordered,~]=sort(SVMtrainRExcumpsthagg(ises).VISp.testtrialinds(:));
            testtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ttindsordered);
            validramptrials = RExtestrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExtestrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

            switch vec2compare
                case 'rampAB'
                    vecA = RExtestTmilestonesagg(ises).(whichvisareaA)(3,:)-RExtestTmilestonesagg(ises).(whichvisareaA)(1,:);
                    vecB = RExtestTmilestonesagg(ises).(whichvisareaB)(3,:)-RExtestTmilestonesagg(ises).(whichvisareaB)(1,:);
                case 'T50AB'
                    vecA = RExtestTmilestonesagg(ises).(whichvisareaA)(2,:);
                    vecB = RExtestTmilestonesagg(ises).(whichvisareaB)(2,:);
                otherwise
                    error([vec2compare ' not recognized'])
            end

            for itt = 1:numel(traintrialtypes)+1
                if itt<=numel(traintrialtypes)
                    trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
                        RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
                else
                    trialsoind = find( validramptrials' & ...
                        RExtestlabelfinalagg(ises).(whichvisareaA)==testtrialstt & RExtestlabelfinalagg(ises).(whichvisareaB)==testtrialstt );
                end
                if numel(trialsoind)<=0
                    continue
                end
                [pwsr,tblwsr,statswsr] = signrank(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                zscore = mean(vecA(trialsoind)-vecB(trialsoind))/std(vecA(trialsoind)-vecB(trialsoind));

                [pmww,tblmww,statsmww] = ranksum(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                [~,~,~,AUC] = perfcurve([ones(size(trialsoind)) zeros(size(trialsoind))], [vecA(trialsoind), vecB(trialsoind)], '1');
                stdpooled = sqrt( ( std(vecA(trialsoind))^2+std(vecB(trialsoind))^2 )/2 );
                CohenD = mean(vecA(trialsoind)-vecB(trialsoind))/stdpooled;

                tempAB.(ABfield).Ntrials(itt,ises) = numel(trialsoind);
                tempAB.(ABfield).zscore(itt,ises) = zscore;
                tempAB.(ABfield).TSwsr(itt,ises) = statswsr.signedrank;
                if isfield(statswsr, 'zval')
                    tempAB.(ABfield).Zwsr(itt,ises) = statswsr.zval;
                end

                tempAB.(ABfield).CohenD(itt,ises) = CohenD;
                tempAB.(ABfield).TSmww(itt,ises) = statsmww.ranksum;
                if isfield(statsmww, 'zval')
                    tempAB.(ABfield).Zmww(itt,ises) = statsmww.zval;
                end
                tempAB.(ABfield).AUC(itt,ises) = AUC;

                tempvecdiffAB.(ABfield){itt,ises} = vecA(trialsoind)-vecB(trialsoind);
            end
        end
        toc
    end
    switch vec2compare
        case 'rampAB'
            RExtestrampAB = tempAB;
            RExtestalldifframpAB = tempvecdiffAB;
        case 'T50AB'
            RExtestalldiffT50AB = tempvecdiffAB;
        otherwise
            error([vec2compare ' not recognized'])
    end
end

figure
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};
    switch vec2compare
        case 'rampAB'
            tempvecdiffAB = RExtestalldifframpAB;
        case 'T50AB'
            tempvecdiffAB = RExtestalldiffT50AB;
        otherwise
            error([vec2compare ' not recognized'])
    end
    for itt = 1:numel(traintrialtypes)+1
        subplot(numel(comparisonvectors), numel(traintrialtypes)+1, (v-1)*(numel(traintrialtypes)+1)+itt)

        tempvec = cat(2,tempvecdiffAB.VISp_VISl{itt,:});
        histogram(tempvec)
        p = signrank(tempvec);
        if itt<=numel(traintrialtypes)
            title(sprintf('%s V1 vs LM Trial %d p=%.4f\n', vec2compare,traintrialtypes(itt),p) )
        else
            title(sprintf('%s V1 vs LM  Pool XRE Test Trials p=%.4f\n', vec2compare,p) )
        end
    end
end

compmetric = 'z-score';
% vec2compare = 'T50AB';
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};

    switch vec2compare
        case 'rampAB'
            tempAB = RExtestrampAB;
        case 'T50AB'
            tempAB = RExtestT50AB;
        otherwise
            error([vec2compare ' not recognized'])
    end
    for itt = 1:numel(traintrialtypes)+1
        if itt<=numel(traintrialtypes)
            fprintf('%s Trial %d\n', vec2compare, traintrialtypes(itt))
        else
            fprintf('%s Pool XRE Test Trials\n', vec2compare)
        end
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
            switch compmetric
                case 'z-score'
                    tempvec = tempAB.(ABfield).zscore(itt,:);
                case 'Znorm'
                    tempvec = tempAB.(ABfield).Zwsr(itt,:)./sqrt( tempAB.(ABfield).Ntrials(itt,:) );
                case 'Zwsr'
                    tempvec = tempAB.(ABfield).Zwsr(itt,:);
                otherwise
                    error([compmetric ' not recognized'])
            end
            p = signrank(tempvec,0);
            pright = signrank(tempvec,0, 'tail', 'right');
            pleft = signrank(tempvec,0, 'tail', 'left');
            fprintf('%s vs %s N=%d %s %.2f (%.2f~%.2f) (p=%.4f, pright=%.4f, pleft=%.4f)\n', whichvisareaA,whichvisareaB, ...
                nnz(~isnan(tempvec)), compmetric, prctile(tempvec,50), prctile(tempvec,25), prctile(tempvec,75), p, pright, pleft)
        end
    end
end

%%
figure
histogram(RExtestrampAB.(ABfield).zscore(itt,:))

n = RExtestrampAB.(ABfield).Ntrials;
figure; plot(RExtestrampAB.(ABfield).AUC, (RExtestrampAB.(ABfield).TSmww -n.*(n+1)/2)./(n.^2), 'o') % exact match
figure; plot(RExtestrampAB.(ABfield).AUC, RExtestrampAB.(ABfield).Zmww, 'o')
corr(RExtestrampAB.(ABfield).AUC(:), RExtestrampAB.(ABfield).Zmww(:), 'rows','complete') % 0.9786
figure; plot(RExtestrampAB.(ABfield).AUC, RExtestrampAB.(ABfield).CohenD, 'o')
corr(RExtestrampAB.(ABfield).AUC(:), RExtestrampAB.(ABfield).CohenD(:), 'rows','complete') % 0.9351

figure; plot(RExtestrampAB.(ABfield).TSmww, RExtestrampAB.(ABfield).TSwsr, 'o')
corr( reshape(RExtestrampAB.(ABfield).TSmww,[],1), reshape(RExtestrampAB.(ABfield).TSwsr,[],1), 'rows','complete') % 0.9351
figure; plot(RExtestrampAB.(ABfield).Zmww, RExtestrampAB.(ABfield).Zwsr, 'o')
corr( reshape(RExtestrampAB.(ABfield).Zmww,[],1), reshape(RExtestrampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9885

figure; plot(RExtestrampAB.(ABfield).zscore, RExtestrampAB.(ABfield).Zwsr, 'o')
corr( reshape(RExtestrampAB.(ABfield).zscore,[],1), reshape(RExtestrampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9706
figure; plot(RExtestrampAB.(ABfield).zscore, RExtestrampAB.(ABfield).Zwsr./sqrt(n), 'o')
corr( reshape(RExtestrampAB.(ABfield).zscore,[],1), reshape(RExtestrampAB.(ABfield).Zwsr./sqrt(n),[],1), 'rows','complete') % 0.9892

figure; plot(RExtestrampAB.(ABfield).zscore, RExtestrampAB.(ABfield).TSwsr, 'o') % not correlated
corr( reshape(RExtestrampAB.(ABfield).zscore,[],1), reshape( RExtestrampAB.(ABfield).TSwsr,[],1), 'rows','complete') % 0.5393
figure; plot(RExtestrampAB.(ABfield).zscore, (RExtestrampAB.(ABfield).TSwsr -n.*(n+1)/2)./(n.^2), 'o')
corr( reshape(RExtestrampAB.(ABfield).zscore,[],1), reshape( (RExtestrampAB.(ABfield).TSwsr -n.*(n+1)/2)./(n.^2),[],1), 'rows','complete') % 0.9869

corr(reshape(RExtestrampAB.(ABfield).AUC,[],1), reshape(RExtestrampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9785

figure; plot(RExtestrampAB.(ABfield).AUC, RExtestrampAB.(ABfield).Zwsr./sqrt(n), 'o')
corr(reshape(RExtestrampAB.(ABfield).AUC,[],1), reshape(RExtestrampAB.(ABfield).Zwsr./sqrt(n),[],1), 'rows','complete') % 0.9886
figure; plot(RExtestrampAB.(ABfield).zscore, RExtestrampAB.(ABfield).CohenD, 'o')
corr(RExtestrampAB.(ABfield).zscore(:), RExtestrampAB.(ABfield).CohenD(:), 'rows','complete') % 0.9969

figure; plot(RExtestrampAB.(ABfield).Zwsr, RExtestrampAB.(ABfield).TSwsr, 'o') % correlated


% % relationship between AUROC (AUC) and ranksum test statistic (TSmww)
% tempauc = SP_gratings(d).(whichR).(preds{ii}).AUC{1};
% n0 = SP_gratings(d).(whichR).(preds{ii}).Ntrials0{1};
% n1 = SP_gratings(d).(whichR).(preds{ii}).Ntrials1{1};
% tempu = SP_gratings(d).(whichR).(preds{ii}).TSmww{1};
% if max(abs((tempu-n1.*(n1+1)/2)-(tempauc.*n0.*n1)))>2^-32
%     error('unexpected mismatch between perfcurve and ranksum results')
% end

%% example session ramp and T50 comparison
ises = 11;
rampA = RExtestTmilestonesagg(ises).(whichvisareaA)(3,:)-RExtestTmilestonesagg(ises).(whichvisareaA)(1,:);
rampB = RExtestTmilestonesagg(ises).(whichvisareaB)(3,:)-RExtestTmilestonesagg(ises).(whichvisareaB)(1,:);
itt = 1;
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    RExtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
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


T50A = RExtestTmilestonesagg(ises).(whichvisareaA)(2,:);
T50B = RExtestTmilestonesagg(ises).(whichvisareaB)(2,:);
figure; hold all
plot(T50A(trialsoind), T50B(trialsoind), 'o')
xl = xlim;
plot(xl,xl, '-')
p = signrank(T50A(trialsoind), T50B(trialsoind));
fprintf('T50 %d trials %s vs %s p=%.4f\n', traintrialtypes(itt), whichvisareaA, whichvisareaB, p)
fprintf('median: %.2f vs %.2f, mean: %.2f vs %.2f\n', ...
    median(T50A(trialsoind)), median(T50B(trialsoind)), ...
    mean(T50A(trialsoind)), mean(T50B(trialsoind)))

%% XRE inference trials
% first, focus on V1 and LM, and on trials where both areas' decoders had correct predictions
ises = 11;
Ninftrials = size(RExinfscorecumpsthagg(ises).(whichvisarea),2);
ICtrialinds = find( ismember(SVMtrainRExcumpsthagg(ises).(whichvisarea).trialorder, ICtrialtypes) );
ICtrialstt = SVMtrainRExcumpsthagg(ises).(whichvisarea).trialorder(ICtrialinds);

figure
% for itt = 1:numel(traintrialtypes)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    subplot(1,numel(visareas),a)
    imagesc(cumpsthtl,1:Ninftrials, RExinfscorecumpsthagg(ises).(whichvisarea)')
    colorbar
    title(sprintf('I_C to X_R_E inference %s', visareas{a}))
end
% end
ises = 11;
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    subplot(1,numel(visareas),a)
    imagesc(cumpsthtl,1:Ninftrials, RExinfnormscorecumpsthagg(ises).(whichvisarea)')
    colorbar
    caxis([0 1])
    title(sprintf('I_C to X_R_E inference %s', visareas{a}))
end

% average across trials
whichvisareaA = 'VISp';
whichvisareaB = 'VISl';
figure
for itt = 1:numel(ICtrialtypes)
    trialsoind = find( ICtrialstt==ICtrialtypes(itt) & ...
        RExinflabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExinflabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
    subplot(1,2,itt)
    hold all
    shadedErrorBar(cumpsthtl, nanmean(RExinfnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
        nanstd(RExinfnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
    shadedErrorBar(cumpsthtl, nanmean(RExinfnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
        nanstd(RExinfnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
    title(sprintf('%d as %d', ICtrialtypes(itt), traintrialtypes(itt) ))
end

%% COMPARE DYNAMICS OF TWO PAREAS TRIAL-BY-TRIAL: inf trials, when both areas had the *correct* prediction
% find the first timepoint at which norminfscorecumpsth crosses
% 0.25, 0.5, 0.75 (called T25, T50, T75 respsectively)
% divide into 6 trial types
% 1. areaA ramping starts earlier and finishes later than area B
% 2. areaB ramping starts earlier and finishes later than area A
% 3. areaA faster than areaB if all three timeponts (T25, T50, T75) are earlier for A than B
% 4. areaB faster than areaA if all three timeponts (T25, T50, T75) are earlier for A than B
% 5. simultaneous: all three timeponts are identical
% 6. crisscrossing if T25 and T75 go in one direction and T50 goes in the opposite direction
dynamicslabels = 0:6;
RExinfcomparedynamicsprob = struct();
RExinfcomparedynamicscnt = struct();
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
    RExinfcomparedynamicsprob.(ABfield) = zeros(numel(traintrialtypes)+1, length(dynamicslabels), Nsessions);
    RExinfcomparedynamicscnt.(ABfield) = zeros(numel(traintrialtypes)+1, length(dynamicslabels), Nsessions);
    for ises = 1:Nsessions
        a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
        if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
            RExinfcomparedynamicsprob.(ABfield)(:,:,ises) = NaN;
            RExinfcomparedynamicscnt.(ABfield)(:,:,ises) = NaN;
            continue
        end
        ICtrialinds = find( ismember(SVMtrainRExcumpsthagg(ises).VISp.trialorder, ICtrialtypes) );
        ICtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ICtrialinds);
        validramptrials = RExinfrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExinfrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

        ICtrialsttconvert = ICtrialstt;
        ICtrialsttconvert(ICtrialstt==106) = 1201;
        ICtrialsttconvert(ICtrialstt==111) = 1299;
        if ~isequal( unique(ICtrialsttconvert), traintrialtypes )
            error('check ICtrialsttconvert')
        end

        Ninfall = size(RExinfnormscorecumpsthagg(ises).(whichvisareaA),2);
        comparedynamics = NaN(1,Ninfall);
        for dyn = 1:6
            switch dyn
                case 1
                    trialsinclass = RExinfTmilestonesagg(ises).(whichvisareaA)(1,:)<=RExinfTmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & RExinfTmilestonesagg(ises).(whichvisareaA)(3,:)>=RExinfTmilestonesagg(ises).(whichvisareaB)(3,:);
                case 2
                    trialsinclass = RExinfTmilestonesagg(ises).(whichvisareaA)(1,:)>=RExinfTmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & RExinfTmilestonesagg(ises).(whichvisareaA)(3,:)<=RExinfTmilestonesagg(ises).(whichvisareaB)(3,:);
                case 3
                    trialsinclass = all(RExinfTmilestonesagg(ises).(whichvisareaA)<=RExinfTmilestonesagg(ises).(whichvisareaB),1);
                case 4
                    trialsinclass = all(RExinfTmilestonesagg(ises).(whichvisareaA)>=RExinfTmilestonesagg(ises).(whichvisareaB),1);
                case 5
                    trialsinclass = all(RExinfTmilestonesagg(ises).(whichvisareaA)==RExinfTmilestonesagg(ises).(whichvisareaB),1);
                case 6
                    temp = RExinfTmilestonesagg(ises).(whichvisareaA)>RExinfTmilestonesagg(ises).(whichvisareaB);
                    trialsinclass = temp(1,:)==~temp(2,:) & temp(3,:)==~temp(2,:);
            end
            comparedynamics(trialsinclass) = dyn;
        end

        for itt = 1:numel(traintrialtypes)+1
            if itt<=numel(traintrialtypes)
                trialsoind = find( validramptrials' & ICtrialstt==ICtrialtypes(itt) & ...
                    RExinflabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExinflabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
            else
                trialsoind = find( validramptrials' & ...
                    RExinflabelfinalagg(ises).(whichvisareaA)==ICtrialsttconvert& RExinflabelfinalagg(ises).(whichvisareaB)==ICtrialsttconvert );
            end
            if numel(trialsoind)<1
                continue
            end
            [v,c]=uniquecnt(comparedynamics(trialsoind));
            % disp([v',c'])
            RExinfcomparedynamicsprob.(ABfield)(itt, ismember(dynamicslabels,v), ises) = c/numel(trialsoind);
            RExinfcomparedynamicscnt.(ABfield)(itt, ismember(dynamicslabels,v), ises) = c;
        end

    end
end

if ~isequaln(sum(RExinfcomparedynamicscnt.(ABfield)(1:numel(traintrialtypes), :,:),1), ...
        RExinfcomparedynamicscnt.(ABfield)(numel(traintrialtypes)+1, :,:) )
    error('check pooling')
end
[0:6; nansum(RExinfcomparedynamicscnt.(ABfield),3)]

IC1dyn4vec = squeeze(RExinfcomparedynamicsprob.VISp_VISl(1,4,:));
IC1dyn3vec = squeeze(RExinfcomparedynamicsprob.VISp_VISl(1,3,:));
fprintf('IC1->XRE1 V1 vs LM dyn4(LM faster) vs dyn3(V1 faster) p=%.4f\n', signrank(IC1dyn4vec, IC1dyn3vec))

ICpooldyn4vec = squeeze(RExinfcomparedynamicsprob.VISp_VISl(3,4,:));
ICpooldyn3vec = squeeze(RExinfcomparedynamicsprob.VISp_VISl(3,3,:));
fprintf('IC->XRE pool V1 vs LM dyn4(LM faster) vs dyn3(V1 faster) p=%.4f\n', signrank(ICpooldyn4vec, ICpooldyn3vec))

% more 4 than 3
% IC1->XRE1 V1 vs LM dyn4(LM faster) vs dyn3(V1 faster) p=0.0156
% IC->XRE pool V1 vs LM dyn4(LM faster) vs dyn3(V1 faster) p=0.0156

ICpooldynmat = squeeze(RExinfcomparedynamicsprob.VISp_VISl(3,2:end,:));
ICpooldynmat(:, any(isnan(ICpooldynmat),1) ) = [];
[p,tbl,stats]=friedman(ICpooldynmat');
figure; multcompare(stats)

figure
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
    for itt = 1:numel(traintrialtypes)+1
        subplot(2,numel(traintrialtypes)+1,(ab-1)*(numel(traintrialtypes)+1)+itt)
        imagesc(squeeze(RExinfcomparedynamicsprob.(ABfield)(itt,:,:))')
        set(gca,'Xtick',1:length(dynamicslabels), 'XtickLabel', dynamicslabels)
        if itt<=numel(traintrialtypes)
            title(sprintf('%s vs %s Trial %d', whichvisareaA, whichvisareaB, ICtrialtypes(itt) ))
        else
            title(sprintf('%s vs %s Pool IC Trials', whichvisareaA, whichvisareaB ))
        end
        colorbar
    end

end
colormap redblue

%% compare ramp time: longer means slower
% calculate AUROC for rampA vs rampB for each session, then see if the
% distribution is significantly different from 0.5 across sessions
comparisonmetrics = {'Ntrials', 'zscore', 'Zwsr', 'TSwsr', 'CohenD', 'Zmww', 'TSmww', 'AUC'};
comparisonvectors = {'rampAB', 'T50AB'};
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};
    tempAB = struct();
    tempvecdiffAB = struct();
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
        tic
        for imet = 1:numel(comparisonmetrics)
            tempAB.(ABfield).(comparisonmetrics{imet}) = NaN(numel(traintrialtypes)+1, Nsessions);
        end
        tempvecdiffAB.(ABfield) = cell(numel(traintrialtypes)+1, Nsessions);
        for ises = 1:Nsessions
            a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
            if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
                continue
            end
            ICtrialinds = find( ismember(SVMtrainRExcumpsthagg(ises).VISp.trialorder, ICtrialtypes) );
            ICtrialstt = SVMtrainRExcumpsthagg(ises).VISp.trialorder(ICtrialinds);
            validramptrials = RExinfrhorampxtimeagg(ises).(whichvisareaA)>threshspear & RExinfrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

            ICtrialsttconvert = ICtrialstt;
            ICtrialsttconvert(ICtrialstt==106) = 1201;
            ICtrialsttconvert(ICtrialstt==111) = 1299;
            if ~isequal( unique(ICtrialsttconvert), traintrialtypes )
                error('check ICtrialsttconvert')
            end

            switch vec2compare
                case 'rampAB'
                    vecA = RExinfTmilestonesagg(ises).(whichvisareaA)(3,:)-RExinfTmilestonesagg(ises).(whichvisareaA)(1,:);
                    vecB = RExinfTmilestonesagg(ises).(whichvisareaB)(3,:)-RExinfTmilestonesagg(ises).(whichvisareaB)(1,:);
                case 'T50AB'
                    vecA = RExinfTmilestonesagg(ises).(whichvisareaA)(2,:);
                    vecB = RExinfTmilestonesagg(ises).(whichvisareaB)(2,:);
                otherwise
                    error([vec2compare ' not recognized'])
            end

            for itt = 1:numel(traintrialtypes)+1
                if itt<=numel(traintrialtypes)
                    trialsoind = find( validramptrials' & ICtrialstt==ICtrialtypes(itt) & ...
                        RExinflabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & RExinflabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
                else
                    trialsoind = find( validramptrials' & ...
                        RExinflabelfinalagg(ises).(whichvisareaA)==ICtrialsttconvert & RExinflabelfinalagg(ises).(whichvisareaB)==ICtrialsttconvert );
                end
                if numel(trialsoind)<=3
                    continue
                end
                [pwsr,tblwsr,statswsr] = signrank(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                zscore = mean(vecA(trialsoind)-vecB(trialsoind))/std(vecA(trialsoind)-vecB(trialsoind));

                [pmww,tblmww,statsmww] = ranksum(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                [~,~,~,AUC] = perfcurve([ones(size(trialsoind)) zeros(size(trialsoind))], [vecA(trialsoind), vecB(trialsoind)], '1');
                stdpooled = sqrt( ( std(vecA(trialsoind))^2+std(vecB(trialsoind))^2 )/2 );
                CohenD = mean(vecA(trialsoind)-vecB(trialsoind))/stdpooled;

                tempAB.(ABfield).Ntrials(itt,ises) = numel(trialsoind);
                tempAB.(ABfield).zscore(itt,ises) = zscore;
                tempAB.(ABfield).TSwsr(itt,ises) = statswsr.signedrank;
                if isfield(statswsr, 'zval')
                    tempAB.(ABfield).Zwsr(itt,ises) = statswsr.zval;
                end

                tempAB.(ABfield).CohenD(itt,ises) = CohenD;
                tempAB.(ABfield).TSmww(itt,ises) = statsmww.ranksum;
                if isfield(statsmww, 'zval')
                    tempAB.(ABfield).Zmww(itt,ises) = statsmww.zval;
                end
                tempAB.(ABfield).AUC(itt,ises) = AUC;

                tempvecdiffAB.(ABfield){itt,ises} = vecA(trialsoind)-vecB(trialsoind);
            end
        end
        toc
    end
    switch vec2compare
        case 'rampAB'
            RExinframpAB = tempAB;
            RExinfalldifframpAB = tempvecdiffAB;
        case 'T50AB'
            RExinfT50AB = tempAB;
            RExinfalldiffT50AB = tempvecdiffAB;
        otherwise
            error([vec2compare ' not recognized'])
    end
end

figure
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};
    switch vec2compare
        case 'rampAB'
            tempvecdiffAB = RExinfalldifframpAB;
        case 'T50AB'
            tempvecdiffAB = RExinfalldiffT50AB;
        otherwise
            error([vec2compare ' not recognized'])
    end
    for itt = 1:numel(traintrialtypes)+1
        subplot(numel(comparisonvectors), numel(traintrialtypes)+1, (v-1)*(numel(traintrialtypes)+1)+itt)

        tempvec = cat(2,tempvecdiffAB.VISp_VISl{itt,:});
        histogram(tempvec)
        p = signrank(tempvec);
        if itt<=numel(traintrialtypes)
            title(sprintf('%s V1 vs LM Trial %d p=%.4f\n', vec2compare,ICtrialtypes(itt),p) )
        else
            title(sprintf('%s V1 vs LM  Pool IC Trials p=%.4f\n', vec2compare,p) )
        end
    end
end

compmetric = 'z-score';
% vec2compare = 'T50AB';
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};

    switch vec2compare
        case 'rampAB'
            tempAB = RExinframpAB;
        case 'T50AB'
            tempAB = RExinfT50AB;
        otherwise
            error([vec2compare ' not recognized'])
    end
    for itt = 1:numel(traintrialtypes)+1
        if itt<=numel(traintrialtypes)
            fprintf('%s Trial %d\n', vec2compare,ICtrialtypes(itt))
        else
            fprintf('%s Pool IC Trials\n', vec2compare)
        end
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
            switch compmetric
                case 'z-score'
                    tempvec = tempAB.(ABfield).zscore(itt,:);
                case 'Znorm'
                    tempvec = tempAB.(ABfield).Zwsr(itt,:)./sqrt( tempAB.(ABfield).Ntrials(itt,:) );
                case 'Zwsr'
                    tempvec = tempAB.(ABfield).Zwsr(itt,:);
                otherwise
                    error([compmetric ' not recognized'])
            end
            p = signrank(tempvec,0);
            pright = signrank(tempvec,0, 'tail', 'right');
            pleft = signrank(tempvec,0, 'tail', 'left');
            fprintf('%s vs %s N=%d %s %.2f (%.2f~%.2f) (p=%.4f, pright=%.4f, pleft=%.4f)\n', whichvisareaA,whichvisareaB, ...
                nnz(~isnan(tempvec)), compmetric, prctile(tempvec,50), prctile(tempvec,25), prctile(tempvec,75), p, pright, pleft)
        end
    end
end

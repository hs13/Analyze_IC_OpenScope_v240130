datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

Twin = 5;
neuopt = 'RS';
svmdesc = 'trainICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

pathsv = [datadir 'SVM_' svmdesc '_selectareas' filesep];
whichICblock = 'ICwcfg1';
whichblock = [whichICblock '_presentations'];
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};

switch svmdesc
    case 'trainICRC'
        traintrialtypes = [106, 107, 110, 111];
        TREtrialtypes =[1105, 1109];
    case 'trainREx'
        traintrialtypes = [1201, 1299];
    otherwise
        error([svmdesc ' not recognized'])
end

SVMtrainICRCcumpsthagg = struct();
ICRCtestlabelfinalagg = struct();
ICRCtestscorecumpsthagg = struct();
ICRCtestnormscorecumpsthagg = struct();
Tmilestonesagg = struct();

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
            SVMtrainICRCcumpsthagg(ises).(whichvisarea) = SVMcumpsth;
        else
            SVMtrainICRCcumpsthagg(ises).(whichvisarea) = [];
        end
    end

    cumpsthtl = SVMtrainICRCcumpsthagg(ises).VISp.psthbinTends;
    [ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
    testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainICRCcumpsthagg(ises).(whichvisarea))
            continue
        end
        finallabel = squeeze(SVMtrainICRCcumpsthagg(ises).(whichvisarea).label(end,:,:));

        Ntesttrials = size(SVMtrainICRCcumpsthagg(ises).(whichvisarea).testtrialinds,1);
        Nsplits = size(SVMtrainICRCcumpsthagg(ises).(whichvisarea).testtrialinds,2);
        origtestlabel = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
        origtestscore = NaN(length(cumpsthtl), Ntesttrials, Nsplits);
        for isplit = 1:Nsplits
            testtrialinds = SVMtrainICRCcumpsthagg(ises).(whichvisarea).testtrialinds(:,isplit);
            origtestlabel(:,:,isplit) = SVMtrainICRCcumpsthagg(ises).(whichvisarea).label(:,testtrialinds,isplit);
            for itt = 1:numel(traintrialtypes)
                temptrials = finallabel(testtrialinds,isplit)==traintrialtypes(itt);
                origtestscore(:,temptrials,isplit) = squeeze(SVMtrainICRCcumpsthagg(ises).(whichvisarea).score(:,testtrialinds(temptrials),itt,isplit));
            end
        end
        origtestlabel = reshape(origtestlabel, length(cumpsthtl), Ntesttrials*Nsplits);
        origtestscore = reshape(origtestscore, length(cumpsthtl), Ntesttrials*Nsplits);

        finaltestlabel = origtestlabel(end,:);

        [tti,testtrialord]=sort(SVMtrainICRCcumpsthagg(ises).(whichvisarea).testtrialinds(:));
        if ~isequal(ttindsordered, tti)
            error('test trial set different for every area?')
        end
        ICRCtestlabelfinalagg(ises).(whichvisarea) = finaltestlabel(testtrialord);
        ICRCtestscorecumpsthagg(ises).(whichvisarea) = origtestscore(:,testtrialord);
    end

    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainICRCcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempts = ICRCtestscorecumpsthagg(ises).(whichvisarea);
        trialsvalramp = tempts(end,:)-tempts(1,:) >0;
        testnormscorecumpsth = (tempts-tempts(1,:))./(tempts(end,:)-tempts(1,:));
        testnormscorecumpsth(:,~trialsvalramp)=NaN;

        ICRCtestnormscorecumpsthagg(ises).(whichvisarea) = testnormscorecumpsth;
    end

    milestonesvec = [0.25; 0.50; 0.75];
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainICRCcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempnorm = ICRCtestnormscorecumpsthagg(ises).(whichvisarea);
        Nalltest = size(tempnorm,2);
        Tmilestonesagg(ises).(whichvisarea) = NaN(length(milestonesvec), Nalltest);
        for t = 1:length(milestonesvec)
            % first timepoint to cross milestone
            [~,mi]=max(tempnorm>=milestonesvec(t),[],1);
            % v50 = tempnorm(sub2ind(size(tempnorm), mi, 1:size(tempnorm,2)));
            % v49 = tempnorm(sub2ind(size(tempnorm), mi-1, 1:size(tempnorm,2)));
            % if ~all(v50>=0.5 & v49<0.5)
            %     error('check mi')
            % end
            Tmilestonesagg(ises).(whichvisarea)(t,:) = cumpsthtl(mi);
        end
        if ~isequal(sort(Tmilestonesagg(ises).(whichvisarea),1), Tmilestonesagg(ises).(whichvisarea))
            error('Tmilestones does not follow order -- check algorithm')
        end
    end

    %%%%%%%%%%%%%%%%%%%% RESUME HERE
    % for inference, chose mode prediction as the final label
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainRExcumpsthagg(ises).(whichvisarea))
            continue
        end
        finallabel = squeeze(SVMtrainICRCcumpsthagg(ises).(whichvisarea).label(end,:,:));

        TREtrialinds = find( ismember(SVMtrainICRCcumpsthagg(ises).(whichvisarea).trialorder, TREtrialtypes) );
        TREtrialstt = SVMtrainICRCcumpsthagg(ises).(whichvisarea).trialorder(TREtrialinds);
        NTREtrials = numel(TREtrialinds);
        originflabel = SVMtrainICRCcumpsthagg(ises).(whichvisarea).label(:,TREtrialinds,:);
        originfscore = SVMtrainICRCcumpsthagg(ises).(whichvisarea).score(:,TREtrialinds,:,:);

        % isequal( squeeze(SVMtrainRExcumpsthagg(ises).(whichvisarea).label(end,:,:)), SVMtrainRExagg.(whichICblock)(ises).(whichvisarea).spkcnt.all.label)
        [M,F,C]=mode(originflabel,3);
        Cnumel = cellfun(@numel, C);

        finalscore = squeeze(originfscore(end,:,:,:));
        finalmodescore = NaN(NTREtrials,Nsplits);
        originfmodecore = NaN(length(cumpsthtl), NTREtrials, Nsplits);
        for typi = 1:numel(traintrialtypes)
            temptrials = M(end,:)==traintrialtypes(typi);
            finalmodescore(temptrials,:) = squeeze(finalscore(temptrials,typi,:));
            originfmodecore(:,temptrials,:) = squeeze(originfscore(:,temptrials,typi,:));
        end
        [mv,maxscoresplit]=max(finalmodescore,[],2);
        mmv=max( reshape(finalscore, NTREtrials,numel(traintrialtypes)*Nsplits),[],2);

        RExinfscorecumpsth = NaN(length(cumpsthtl), NTREtrials);
        for isplit = 1:Nsplits
            temptrials = maxscoresplit==isplit;
            RExinfscorecumpsth(:,temptrials) = squeeze( originfmodecore(:,temptrials,isplit) );
        end

        trialsvalramp = RExinfscorecumpsth(end,:)-RExinfscorecumpsth(1,:)>0;
        RExinfnormscorecumpsth = (RExinfscorecumpsth-RExinfscorecumpsth(1,:))./(RExinfscorecumpsth(end,:)-RExinfscorecumpsth(1,:));
        RExinfnormscorecumpsth(:,~trialsvalramp)=NaN;

        RExinflabelfinalagg(ises).(whichvisarea) = M(end,:);
        % ICRCinflabelcumpsthagg(ises).(whichvisarea) = M;
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
        NTREtrials = size(tempnorm,2);
        RExinfTmilestones = NaN(length(milestonesvec), NTREtrials);
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

save([pathsv 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc '_' neuopt 'agg.mat'], 'SVMtrainICRCcumpsthagg', ...
    'ICRCtestlabelfinalagg', 'ICRCtestscorecumpsthagg', 'ICRCtestnormscorecumpsthagg', 'Tmilestonesagg', '-v7.3')

%%
load('G:\My Drive\RESEARCH\ICexpts_revision23\openscope_HR_SVMtrainICRC_zscore_agg.mat', 'Nneuronsperarea')
discardbelowNneurons = 50;

%% valid ramp trials: define based on spearman correlation between timepoints and scores
% set an arbitrary threshold at 0.5 to determine validity.
% alternatively, could use a more inclusive threshold value of 0
threshspear = 0.5;
testrhorampxtimeagg = struct();
for ises = 1:Nsessions
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        if isempty(SVMtrainICRCcumpsthagg(ises).(whichvisarea))
            continue
        end
        tempts = ICRCtestscorecumpsthagg(ises).(whichvisarea);
        testrhorampxtimeagg(ises).(whichvisarea) = corr(tempts, (1:length(cumpsthtl))', 'type', 'Spearman');
    end
end

rhorampxtime = cat(2,testrhorampxtimeagg.(whichvisarea));
figure; histogram(rhorampxtime(:), 'binwidth', 0.01)

tempts = ICRCtestscorecumpsthagg(ises).(whichvisarea);
rhoampxtimevec = testrhorampxtimeagg(ises).(whichvisarea);
figure
hold all
plot(cumpsthtl, tempts(:,rhoampxtimevec>0.5), 'b-')
plot(cumpsthtl, tempts(:,rhoampxtimevec<=0.5), 'r-')

%% check correspondence between SVMall and final timepoint of SVMcumpsthall
SVMcumpsthall = SVMtrainICRCcumpsthagg(ises);
pathsvm = [pathsv nwbsessions{ises} filesep];

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
SVMcumpsthall = SVMtrainICRCcumpsthagg(ises);
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

%% IC test trials
% first, focus on V1 and LM, and on trials where both areas' decoders had correct predictions
ises = 11;
[ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);

itt =1;
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    trialsoi = testtrialstt==traintrialtypes(itt);
    subplot(1,numel(visareas),a)
    imagesc(cumpsthtl,1:nnz(trialsoi), ICRCtestscorecumpsthagg(ises).(whichvisarea)(:,trialsoi)')
    colorbar
end

% check that reordering went well
Ntt2p = 4;
tt2p = randperm(length(testtrialstt),Ntt2p);
figure
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for t = 1:Ntt2p
        itrial = tt2p(t);
        subplot(Ntt2p, numel(visareas), (t-1)*numel(visareas)+a)
        hold all
        plot(cumpsthtl, ICRCtestscorecumpsthagg(ises).(whichvisarea)(:,itrial), 'k--', 'linewidth', 1)
        [r,c]=find(SVMtrainICRCcumpsthagg(ises).(whichvisarea).testtrialinds==ttindsordered(itrial));
        tempcumpsth = squeeze(SVMtrainICRCcumpsthagg(ises).(whichvisarea).score(:,ttindsordered(itrial),:,c));
        plot(cumpsthtl, tempcumpsth)
        title(sprintf('%s test trial #%d', whichvisarea, itrial))
    end
end

whichvisareaA = 'VISp';
whichvisareaB = 'VISal';
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );

% example trials
figure
for ii = 1:12
    subplot(3,4,ii)
    hold all
    plot(cumpsthtl, ICRCtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind(ii)))
    plot(cumpsthtl, ICRCtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind(ii)))
end

% average across trials
figure
for itt = 1:numel(traintrialtypes)
    trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
        ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
    subplot(2,2,itt)
    hold all
    shadedErrorBar(cumpsthtl, mean(ICRCtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
        std(ICRCtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
    shadedErrorBar(cumpsthtl, mean(ICRCtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
        std(ICRCtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
    title(traintrialtypes(itt))
end

%% mean ramppsth
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
    [ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
    testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);
    validramptrials = testrhorampxtimeagg(ises).(whichvisareaA)>threshspear & testrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

    for itt = 1:numel(traintrialtypes)
        trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
            ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
        subplot(numel(traintrialtypes)*2,Nsessions, (itt-1)*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(ICRCtestscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(ICRCtestscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(ICRCtestscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(ICRCtestscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('test score')
        title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))

        subplot(numel(traintrialtypes)*2,Nsessions, (itt-1+numel(traintrialtypes))*Nsessions+ises )
        hold all
        shadedErrorBar(cumpsthtl, nanmean(ICRCtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),2), ...
            nanstd(ICRCtestnormscorecumpsthagg(ises).(whichvisareaA)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'b', 'LineWidth', 2},1)
        shadedErrorBar(cumpsthtl, nanmean(ICRCtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),2), ...
            nanstd(ICRCtestnormscorecumpsthagg(ises).(whichvisareaB)(:,trialsoind),0,2)/sqrt(numel(trialsoind)), {'Color', 'r', 'LineWidth', 2},1)
        ylabel('norm test score')
        title(sprintf('Ses%d Trial%d', ises, traintrialtypes(itt) ))
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
comparedynamicsprob = struct();
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
    comparedynamicsprob.(ABfield) = zeros(numel(traintrialtypes), length(dynamicslabels), Nsessions);
    for ises = 1:Nsessions
        a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
        if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
            comparedynamicsprob.(ABfield)(:,:,ises) = NaN;
            continue
        end
        Ntestall = size(ICRCtestnormscorecumpsthagg(ises).(whichvisareaA),2);
        comparedynamics = NaN(1,Ntestall);
        for dyn = 1:6
            switch dyn
                case 1
                    trialsinclass = Tmilestonesagg(ises).(whichvisareaA)(1,:)<=Tmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & Tmilestonesagg(ises).(whichvisareaA)(3,:)>=Tmilestonesagg(ises).(whichvisareaB)(3,:);
                case 2
                    trialsinclass = Tmilestonesagg(ises).(whichvisareaA)(1,:)>=Tmilestonesagg(ises).(whichvisareaB)(1,:) ...
                        & Tmilestonesagg(ises).(whichvisareaA)(3,:)<=Tmilestonesagg(ises).(whichvisareaB)(3,:);
                case 3
                    trialsinclass = all(Tmilestonesagg(ises).(whichvisareaA)<=Tmilestonesagg(ises).(whichvisareaB),1);
                case 4
                    trialsinclass = all(Tmilestonesagg(ises).(whichvisareaA)>=Tmilestonesagg(ises).(whichvisareaB),1);
                case 5
                    trialsinclass = all(Tmilestonesagg(ises).(whichvisareaA)==Tmilestonesagg(ises).(whichvisareaB),1);
                case 6
                    temp = Tmilestonesagg(ises).(whichvisareaA)>Tmilestonesagg(ises).(whichvisareaB);
                    trialsinclass = temp(1,:)==~temp(2,:) & temp(3,:)==~temp(2,:);
            end
            comparedynamics(trialsinclass) = dyn;
        end

        [ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
        testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);
        validramptrials = testrhorampxtimeagg(ises).(whichvisareaA)>threshspear & testrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

        for itt = 1:numel(traintrialtypes)
            trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
                ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
            [v,c]=uniquecnt(comparedynamics(trialsoind));
            % disp([v',c'])
            comparedynamicsprob.(ABfield)(itt, ismember(dynamicslabels,v), ises) = c/numel(trialsoind);
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
    for itt = 1:numel(traintrialtypes)
        subplot(2,numel(traintrialtypes),(ab-1)*numel(traintrialtypes)+itt)
        imagesc(squeeze(comparedynamicsprob.(ABfield)(itt,:,:))')
        set(gca,'Xtick',1:length(dynamicslabels), 'XtickLabel', dynamicslabels)
        title(sprintf('%s vs %s Trial %d', whichvisareaA, whichvisareaB, traintrialtypes(itt) ))
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
            tempAB.(ABfield).(comparisonmetrics{imet}) = NaN(numel(traintrialtypes), Nsessions);
        end
        for ises = 1:Nsessions
            a = find(strcmp(visareas, whichvisareaA)); b = find(strcmp(visareas, whichvisareaB));
            if Nneuronsperarea(ises,a)<discardbelowNneurons || Nneuronsperarea(ises,b)<discardbelowNneurons
                continue
            end
            switch vec2compare
                case 'rampAB'
                    vecA = Tmilestonesagg(ises).(whichvisareaA)(3,:)-Tmilestonesagg(ises).(whichvisareaA)(1,:);
                    vecB = Tmilestonesagg(ises).(whichvisareaB)(3,:)-Tmilestonesagg(ises).(whichvisareaB)(1,:);
                case 'T50AB'
                    vecA = Tmilestonesagg(ises).(whichvisareaA)(2,:);
                    vecB = Tmilestonesagg(ises).(whichvisareaB)(2,:);
                otherwise
                    error([vec2compare ' not recognized'])
            end

            [ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
            testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);
            validramptrials = testrhorampxtimeagg(ises).(whichvisareaA)>threshspear & testrhorampxtimeagg(ises).(whichvisareaB)>threshspear;

            for itt = 1:numel(traintrialtypes)
                trialsoind = find( validramptrials' & testtrialstt==traintrialtypes(itt) & ...
                    ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
                [pwsr,tblwsr,statswsr] = signrank(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                zscore = mean(vecA(trialsoind)-vecB(trialsoind))/std(vecA(trialsoind)-vecB(trialsoind));

                [pmww,tblmww,statsmww] = ranksum(vecA(trialsoind), vecB(trialsoind), 'method', 'approximate');
                [~,~,~,AUC] = perfcurve([ones(size(trialsoind)) zeros(size(trialsoind))], [vecA(trialsoind), vecB(trialsoind)], '1');
                stdpooled = sqrt( ( std(vecA(trialsoind))^2+std(vecB(trialsoind))^2 )/2 );
                CohenD = mean(vecA(trialsoind)-vecB(trialsoind))/stdpooled;

                tempAB.(ABfield).Ntrials(itt,ises) = numel(trialsoind);
                tempAB.(ABfield).zscore(itt,ises) = zscore;
                tempAB.(ABfield).TSwsr(itt,ises) = statswsr.signedrank;
                tempAB.(ABfield).Zwsr(itt,ises) = statswsr.zval;

                tempAB.(ABfield).CohenD(itt,ises) = CohenD;
                tempAB.(ABfield).TSmww(itt,ises) = statsmww.ranksum;
                tempAB.(ABfield).Zmww(itt,ises) = statsmww.zval;
                tempAB.(ABfield).AUC(itt,ises) = AUC;
            end
        end
        toc
    end
    switch vec2compare
        case 'rampAB'
            rampAB = tempAB;
        case 'T50AB'
            T50AB = tempAB;
        otherwise
            error([vec2compare ' not recognized'])
    end
end

compmetric = 'z-score';
% vec2compare = 'T50AB';
for v = 1:numel(comparisonvectors)
    vec2compare = comparisonvectors{v};
    switch vec2compare
        case 'rampAB'
            tempAB = rampAB;
        case 'T50AB'
            tempAB = T50AB;
        otherwise
            error([vec2compare ' not recognized'])
    end
    for itt = 1:4
        fprintf('%s Trial %d\n', vec2compare, traintrialtypes(itt))
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
histogram(rampAB.(ABfield).zscore(itt,:))

n = rampAB.(ABfield).Ntrials;
figure; plot(rampAB.(ABfield).AUC, (rampAB.(ABfield).TSmww -n.*(n+1)/2)./(n.^2), 'o') % exact match
figure; plot(rampAB.(ABfield).AUC, rampAB.(ABfield).Zmww, 'o')
corr(rampAB.(ABfield).AUC(:), rampAB.(ABfield).Zmww(:), 'rows','complete') % 0.9786
figure; plot(rampAB.(ABfield).AUC, rampAB.(ABfield).CohenD, 'o')
corr(rampAB.(ABfield).AUC(:), rampAB.(ABfield).CohenD(:), 'rows','complete') % 0.9351

figure; plot(rampAB.(ABfield).TSmww, rampAB.(ABfield).TSwsr, 'o')
corr( reshape(rampAB.(ABfield).TSmww,[],1), reshape(rampAB.(ABfield).TSwsr,[],1), 'rows','complete') % 0.9351
figure; plot(rampAB.(ABfield).Zmww, rampAB.(ABfield).Zwsr, 'o')
corr( reshape(rampAB.(ABfield).Zmww,[],1), reshape(rampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9885

figure; plot(rampAB.(ABfield).zscore, rampAB.(ABfield).Zwsr, 'o')
corr( reshape(rampAB.(ABfield).zscore,[],1), reshape(rampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9706
figure; plot(rampAB.(ABfield).zscore, rampAB.(ABfield).Zwsr./sqrt(n), 'o')
corr( reshape(rampAB.(ABfield).zscore,[],1), reshape(rampAB.(ABfield).Zwsr./sqrt(n),[],1), 'rows','complete') % 0.9892

figure; plot(rampAB.(ABfield).zscore, rampAB.(ABfield).TSwsr, 'o') % not correlated
corr( reshape(rampAB.(ABfield).zscore,[],1), reshape( rampAB.(ABfield).TSwsr,[],1), 'rows','complete') % 0.5393
figure; plot(rampAB.(ABfield).zscore, (rampAB.(ABfield).TSwsr -n.*(n+1)/2)./(n.^2), 'o')
corr( reshape(rampAB.(ABfield).zscore,[],1), reshape( (rampAB.(ABfield).TSwsr -n.*(n+1)/2)./(n.^2),[],1), 'rows','complete') % 0.9869

corr(reshape(rampAB.(ABfield).AUC,[],1), reshape(rampAB.(ABfield).Zwsr,[],1), 'rows','complete') % 0.9785

figure; plot(rampAB.(ABfield).AUC, rampAB.(ABfield).Zwsr./sqrt(n), 'o')
corr(reshape(rampAB.(ABfield).AUC,[],1), reshape(rampAB.(ABfield).Zwsr./sqrt(n),[],1), 'rows','complete') % 0.9886
figure; plot(rampAB.(ABfield).zscore, rampAB.(ABfield).CohenD, 'o')
corr(rampAB.(ABfield).zscore(:), rampAB.(ABfield).CohenD(:), 'rows','complete') % 0.9969

figure; plot(rampAB.(ABfield).Zwsr, rampAB.(ABfield).TSwsr, 'o') % correlated

% % relationship between AUROC (AUC) and ranksum test statistic (TSmww)
% tempauc = SP_gratings(d).(whichR).(preds{ii}).AUC{1};
% n0 = SP_gratings(d).(whichR).(preds{ii}).Ntrials0{1};
% n1 = SP_gratings(d).(whichR).(preds{ii}).Ntrials1{1};
% tempu = SP_gratings(d).(whichR).(preds{ii}).TSmww{1};
% if max(abs((tempu-n1.*(n1+1)/2)-(tempauc.*n0.*n1)))>2^-32
%     error('unexpected mismatch between perfcurve and ranksum results')
% end

%%
ises = 11;
[ttindsordered,~]=sort(SVMtrainICRCcumpsthagg(ises).VISp.testtrialinds(:));
testtrialstt = SVMtrainICRCcumpsthagg(ises).VISp.trialorder(ttindsordered);

rampA = Tmilestonesagg(ises).(whichvisareaA)(3,:)-Tmilestonesagg(ises).(whichvisareaA)(1,:);
rampB = Tmilestonesagg(ises).(whichvisareaB)(3,:)-Tmilestonesagg(ises).(whichvisareaB)(1,:);
itt = 4;
trialsoind = find( testtrialstt==traintrialtypes(itt) & ...
    ICRCtestlabelfinalagg(ises).(whichvisareaA)==traintrialtypes(itt) & ICRCtestlabelfinalagg(ises).(whichvisareaB)==traintrialtypes(itt) );
figure; hold all
plot(rampA(trialsoind), rampB(trialsoind), 'o')
xlabel(['ramp ' whichvisareaA])
ylabel(['ramp ' whichvisareaB])
xl = xlim;
plot(xl,xl, '-')
p = signrank(rampA(trialsoind), rampB(trialsoind));
fprintf('ramp delay (ms) %d trials %s vs %s p=%.4f\n', traintrialtypes(itt), ...
    whichvisareaA, whichvisareaB, p)
fprintf('median: %.2f vs %.2f, mean: %.2f vs %.2f\n', ...
    median(rampA(trialsoind)), median(rampB(trialsoind)), ...
    mean(rampA(trialsoind)), mean(rampB(trialsoind)))


T50A = Tmilestonesagg(ises).(whichvisareaA)(2,:);
T50B = Tmilestonesagg(ises).(whichvisareaB)(2,:);
figure; hold all
plot(T50A(trialsoind), T50B(trialsoind), 'o')
xlabel(['T50 ' whichvisareaA])
ylabel(['T50 ' whichvisareaB])
xl = xlim;
plot(xl,xl, '-')
p = signrank(T50A(trialsoind), T50B(trialsoind));
fprintf('T50 %d trials %s vs %s p=%.4f\n', traintrialtypes(itt), whichvisareaA, whichvisareaB, p)
fprintf('median: %.2f vs %.2f, mean: %.2f vs %.2f\n', ...
    median(T50A(trialsoind)), median(T50B(trialsoind)), ...
    mean(T50A(trialsoind)), mean(T50B(trialsoind)))


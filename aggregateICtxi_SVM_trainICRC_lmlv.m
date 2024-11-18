if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};
% nwbsessions = nwbsessions(7:end);
Nsessions = numel(nwbsessions);

% if 0, keep all neurons; if 1, exclude zero variance neurons in train trial
% types; if 2 exclude zero variance neurons in all trial types
excludeneuvar0 = 0;
fprintf('neuron exclusion criterion %d\n', excludeneuvar0)

silrandagg = struct();
similagg = struct();
for ises = 1:numel(nwbsessions)
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    pltses = false;
    switch excludeneuvar0
        case 0
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes_incl.mat');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes_incl.mat');
        case 1
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes.mat_excltt');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes.mat_excltt');
        case 2
            silrandfn = strcat(pathpp, 'randomsilencing_SVM_trainICRC_lmlvslopes_excl.mat');
            similfn = strcat(pathpp, 'scoresimilarity_SVM_trainICRC_lmlvslopes_excl.mat');
        otherwise
            error('excludeneuvar0 option not recognized')
    end
    silrand = load(silrandfn);
    simil = load(similfn);
    if isempty(fieldnames(silrandagg))
        silrandagg = silrand;
        similagg = simil;
    else
        silrandagg(ises) = silrand;
        similagg(ises) = simil;
    end
end

%% SVM test and inference performance across LMLV slopes
disperses = silrand.disperses;
testperflmlvsagg = cat(1,silrandagg.testperflmlvs);
infperflmlvsagg = cat(1,silrandagg.infperflmlvs);
figure('Position', [100 100 1200 600])
subplot(2,3,1)
hold all
plot(disperses, testperflmlvsagg)
errorbar(disperses, mean(testperflmlvsagg,1), std(testperflmlvsagg,0,1)/sqrt(Nsessions), 'ko-', 'linewidth',2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(cat(1,silrandagg.testperfasis))*[1 1], 'c-')
xlabel('log(mean) vs log(var) slopes')
ylabel('test accuracy')
subplot(2,3,4)
hold all
plot(disperses, infperflmlvsagg)
errorbar(disperses, mean(infperflmlvsagg,1), std(infperflmlvsagg,0,1)/sqrt(Nsessions), 'ko-', 'linewidth',2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(cat(1,silrandagg.infperfasis))*[1 1], 'c-')
xlabel('log(mean) vs log(var) slopes')
ylabel('inference performance')

for ii = 1:4
    switch ii
        case 1
            ylab = 'P(IC1|TRE1)';
            r=1; c=1;
            isp = 5;
        case 2
            ylab = 'P(LC1|TRE1)';
            r=1; c=2;
            isp = 2;
        case 3
            ylab = 'P(LC2|TRE2)';
            r=2; c=3;
            isp = 3;
        case 4
            ylab = 'P(IC2|TRE2)';
            r=2; c=4;
            isp = 6;
    end
    subplot(2,3,isp)
    hold all
    temp = cat(4,silrandagg.infcmlmlvs);
    tempagg = squeeze(temp(r,c,:,:))';
    plot(disperses, tempagg)
    errorbar( disperses, mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), 'ko-', 'linewidth',2)
    xl = [disperses(1) disperses(end)];
    tempasis = cat(3,silrandagg.infcmasis);
    plot(xl, squeeze(mean(tempasis(r,c,:),3))*[1 1], 'c-')
    xlabel('log(mean) vs log(var) slopes')
    ylabel(ylab)
end

%% performance change when silencing random portion of neurons
propneusilvec = silrand.propneusilvec;
lmlvcols = cool(numel(disperses));
lmlvlegs = cell(1,numel(disperses)+1);

figure
hold all
for islope = 1:numel(disperses)
    if disperses(islope)==1
        lw = 2;
    else
        lw = 1;
    end
    tempagg = NaN(Nsessions, numel(propneusilvec));
    for ises = 1:Nsessions
        tempagg(ises,:) = nanmean(silrandagg(ises).silrandtestperflmlvs{islope},1);
    end
    errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
        'Color', lmlvcols(islope,:), 'LineWidth', lw)
    lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
end
tempagg = NaN(Nsessions, numel(propneusilvec));
for ises = 1:Nsessions
    tempagg(ises,:) = nanmean(silrandagg(ises).silrandtestperfasis,1);
end
errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
    'Color', 'k', 'LineWidth', 2)
lmlvlegs{end} = 'as-is';
legend(lmlvlegs)
xlabel('%Neurons')
ylabel('test accuracy')

figure
hold all
for islope = 1:numel(disperses)
    if disperses(islope)==1
        lw = 2;
    else
        lw = 1;
    end
    tempagg = NaN(Nsessions, numel(propneusilvec));
    for ises = 1:Nsessions
        tempagg(ises,:) = nanmean(silrandagg(ises).silrandinfperflmlvs{islope},1);
    end
    errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
        'Color', lmlvcols(islope,:), 'LineWidth', lw)
    lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
end
tempagg = NaN(Nsessions, numel(propneusilvec));
for ises = 1:Nsessions
    tempagg(ises,:) = nanmean(silrandagg(ises).silrandinfperfasis,1);
end
errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
    'Color', 'k', 'LineWidth', 2)
lmlvlegs{end} = 'as-is';
legend(lmlvlegs)
xlabel('%Neurons')
ylabel('inference performance')

figure
for isp = 1:4
    switch isp
        case 1
            ylab = 'P(IC1|TRE1)';
            r=1; c=1;
        case 2
            ylab = 'P(LC1|TRE1)';
            r=1; c=2;
        case 3
            ylab = 'P(LC2|TRE2)';
            r=2; c=3;
        case 4
            ylab = 'P(IC2|TRE2)';
            r=2; c=4;
    end
    subplot(2,2,isp)
    hold all
    for islope = 1:numel(disperses)
        if disperses(islope)==1
            lw = 2;
        else
            lw = 1;
        end
        tempagg = NaN(Nsessions, numel(propneusilvec));
        for ises = 1:Nsessions
            tempagg(ises,:) = nanmean(silrandagg(ises).silrandinfcmlmlvs{islope}(:,:,r,c),1);
        end
        errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
            'Color', lmlvcols(islope,:), 'LineWidth', lw)
        lmlvlegs{islope} = sprintf('Slope%.1f',disperses(islope));
    end
    tempagg = NaN(Nsessions, numel(propneusilvec));
    for ises = 1:Nsessions
        tempagg(ises,:) = nanmean(silrandagg(ises).silrandinfcmasis(:,:,r,c),1);
    end
    errorbar( 100*(1-propneusilvec), mean(tempagg,1), std(tempagg,0,1)/sqrt(Nsessions), ...
        'Color', 'k', 'LineWidth', 2)
    lmlvlegs{end} = 'as-is';
    legend(lmlvlegs)
    xlabel('%Neurons')
    ylabel(ylab)
end

%% SVM score Spearman correlation (rank similarity) each cross-validation fold with trial mean vec
hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];
testt = [106,107,110,111];
inferencett = [1105 1109];
Nhireptt = numel(hireptt);
Ntt = numel(testt);
disperses = simil.disperses;
propneusilvec = simil.propneusilvec;
ttoind = find(~ismember(hireptt, testt));

whichstat = 'avg';

meanvecscorerholmlvsagg = cat(1, similagg.meanvecscorerholmlvs);
meanvecscorerholmlvsaggses = squeeze( mean(cat(3,meanvecscorerholmlvsagg.(whichstat)),1) )';
meanvecscorerhoasisagg = cat(1, similagg.meanvecscorerhoasis);
meanvecscorerhoasisaggses = squeeze( mean(cat(2,meanvecscorerhoasisagg.(whichstat)),1) )';

rhoscorelmlvsagg = cat(1, similagg.rhoscorelmlvs);
rhoscorelmlvstestagg = cat(1, rhoscorelmlvsagg.test);
rhoscorelmlvstestaggstat = cat(5, rhoscorelmlvstestagg.(whichstat));
rhoscorelmlvstestaggses = squeeze(mean(rhoscorelmlvstestaggstat(1,propneusilvec==0,:,:,:),3))';
rhoscoreasisagg = cat(1, similagg.rhoscoreasis);
rhoscoreasistestagg = cat(1, rhoscoreasisagg.test);
rhoscoreasistestaggstat = cat(4, rhoscoreasistestagg.(whichstat));
rhoscoreasistestaggses = squeeze(mean(rhoscoreasistestaggstat(1,propneusilvec==0,:,:),3))';
rhoscorelmlvssimilagg = cat(1, rhoscorelmlvsagg.simil);
rhoscorelmlvssimilaggstat = cat(5, rhoscorelmlvssimilagg.(whichstat));
rhoscorelmlvssimilaggses = squeeze(rhoscorelmlvssimilaggstat(1,propneusilvec==0,:,:,:));
rhoscoreasissimilagg = cat(1, rhoscoreasisagg.simil);
rhoscoreasissimilaggstat = cat(4, rhoscoreasissimilagg.(whichstat));
rhoscoreasissimilaggses = squeeze(rhoscoreasissimilaggstat(1,propneusilvec==0,:,:,:));

figure
subplot(3,4,1)
hold all
pl = plot(disperses,  meanvecscorerholmlvsaggses);
errorbar(disperses, mean(meanvecscorerholmlvsaggses,1), std(meanvecscorerholmlvsaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(meanvecscorerhoasisaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('trial mean vector score consistency')

subplot(3,4,2)
hold all
pl = plot(disperses,  rhoscorelmlvstestaggses);
errorbar(disperses, mean(rhoscorelmlvstestaggses,1), std(rhoscorelmlvstestaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoscoreasistestaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('test trials vs trial mean')

for ii = 1:numel(ttoind)
    rhoscorelmlvsttaggses = squeeze(rhoscorelmlvssimilaggses(ttoind(ii),:,:))';
    rhoscoreasisttaggses = squeeze(rhoscoreasissimilaggses(ttoind(ii),:))';
subplot(3,4,2+ii)
hold all
pl = plot(disperses,  rhoscorelmlvsttaggses);
errorbar(disperses, mean(rhoscorelmlvsttaggses,1), std(rhoscorelmlvsttaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoscoreasisttaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title(sprintf('Trial%d vs trial mean', hireptt(ttoind(ii))))
end

%% SVM score Spearman correlation (rank similarity) across cross-validation folds
whichstat = 'avg';

meanvecscorerholmlvsagg = cat(1, similagg.meanvecscorerholmlvs);
meanvecscorerholmlvsaggses = squeeze( mean(cat(3,meanvecscorerholmlvsagg.(whichstat)),1) )';
meanvecscorerhoasisagg = cat(1, similagg.meanvecscorerhoasis);
meanvecscorerhoasisaggses = squeeze( mean(cat(2,meanvecscorerhoasisagg.(whichstat)),1) )';

rhoscorelmlvsagg = cat(1, similagg.rhoscorelmlvs);
rhoscorelmlvstestagg = cat(1, rhoscorelmlvsagg.testpair);
rhoscorelmlvstestaggstat = cat(5, rhoscorelmlvstestagg.(whichstat));
rhoscorelmlvstestaggses = squeeze(mean(rhoscorelmlvstestaggstat(1,propneusilvec==0,:,:,:),3))';
rhoscoreasisagg = cat(1, similagg.rhoscoreasis);
rhoscoreasistestagg = cat(1, rhoscoreasisagg.testpair);
rhoscoreasistestaggstat = cat(4, rhoscoreasistestagg.(whichstat));
rhoscoreasistestaggses = squeeze(mean(rhoscoreasistestaggstat(1,propneusilvec==0,:,:),3))';
rhoscorelmlvssimilagg = cat(1, rhoscorelmlvsagg.similpair);
rhoscorelmlvssimilaggstat = cat(5, rhoscorelmlvssimilagg.(whichstat));
rhoscorelmlvssimilaggses = squeeze(rhoscorelmlvssimilaggstat(1,propneusilvec==0,:,:,:));
rhoscoreasissimilagg = cat(1, rhoscoreasisagg.similpair);
rhoscoreasissimilaggstat = cat(4, rhoscoreasissimilagg.(whichstat));
rhoscoreasissimilaggses = squeeze(rhoscoreasissimilaggstat(1,propneusilvec==0,:,:));

figure
subplot(3,4,1)
hold all
pl = plot(disperses,  meanvecscorerholmlvsaggses);
errorbar(disperses, mean(meanvecscorerholmlvsaggses,1), std(meanvecscorerholmlvsaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(meanvecscorerhoasisaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('trial mean vector score consistency')

subplot(3,4,2)
hold all
pl = plot(disperses,  rhoscorelmlvstestaggses);
errorbar(disperses, mean(rhoscorelmlvstestaggses,1), std(rhoscorelmlvstestaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoscoreasistestaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('test trial pairs')

for ii = 1:numel(ttoind)
    rhoscorelmlvsttaggses = squeeze(rhoscorelmlvssimilaggses(ttoind(ii),:,:))';
    rhoscoreasisttaggses = squeeze(rhoscoreasissimilaggses(ttoind(ii),:))';
subplot(3,4,2+ii)
hold all
pl = plot(disperses,  rhoscorelmlvsttaggses);
errorbar(disperses, mean(rhoscorelmlvsttaggses,1), std(rhoscorelmlvsttaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoscoreasisttaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title(sprintf('Trial%d pairs', hireptt(ttoind(ii))))
end

%% SVM score Spearman correlation (rank similarity) each random subset (per fixed proportion) with trial mean vec
hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];
testt = [106,107,110,111];
inferencett = [1105 1109];
Nhireptt = numel(hireptt);
Ntt = numel(testt);
ttoind = find(~ismember(hireptt, testt));
disperses = simil.disperses;
propneusilvec = simil.propneusilvec;
iprop = propneusilvec==0.5;

whichstat = 'prct';

rhoxneusublmlvsagg = cat(1, similagg.rhoxneusublmlvs);
rhoxneusublmlvstestagg = cat(1, rhoxneusublmlvsagg.test);
rhoxneusublmlvstestaggstat = cat(5, rhoxneusublmlvstestagg.(whichstat));
rhoxneusublmlvstestaggses = squeeze(mean(rhoxneusublmlvstestaggstat(:,iprop,:,:,:),[1,3]))';
rhoxneusubasisagg = cat(1, similagg.rhoxneusubasis);
rhoxneusubasistestagg = cat(1, rhoxneusubasisagg.test);
rhoxneusubasistestaggstat = cat(4, rhoxneusubasistestagg.(whichstat));
rhoxneusubasistestaggses = squeeze(mean(rhoxneusubasistestaggstat(:,iprop,:,:),[1,3]))';

rhoxneusublmlvssimilagg = cat(1, rhoxneusublmlvsagg.simil);
rhoxneusublmlvssimilaggstat = cat(5, rhoxneusublmlvssimilagg.(whichstat));
rhoxneusublmlvssimilaggses = squeeze(mean(rhoxneusublmlvssimilaggstat(:,iprop,:,:,:),1));
rhoxneusubasissimilagg = cat(1, rhoxneusubasisagg.simil);
rhoxneusubasissimilaggstat = cat(4, rhoxneusubasissimilagg.(whichstat));
rhoxneusubasissimilaggses = squeeze(mean(rhoxneusubasissimilaggstat(:,iprop,:,:),1));

figure
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', sprintf('randomly silence %.0f%% of neurons: SVM score consistency between trial mean vector and neuronal subsets', 100*propneusilvec(iprop)), 'edgecolor', 'none')
subplot(3,3,1)
hold all
pl = plot(disperses,  rhoxneusublmlvstestaggses);
errorbar(disperses, mean(rhoxneusublmlvstestaggses,1), std(rhoxneusublmlvstestaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoxneusubasistestaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('test trials vs trial mean')

for ii = 1:numel(ttoind)
    rhoxneusublmlvsttaggses = (reshape(rhoxneusublmlvssimilaggses(ttoind(ii),:,:),numel(disperses),Nsessions))';
    rhoxneusubasisttaggses = (reshape(rhoxneusubasissimilaggses(ttoind(ii),:),1,Nsessions))';
subplot(3,3,1+ii)
hold all
pl = plot(disperses,  rhoxneusublmlvsttaggses);
errorbar(disperses, mean(rhoxneusublmlvsttaggses,1), std(rhoxneusublmlvsttaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoxneusubasisttaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title(sprintf('Trial%d vs trial mean', hireptt(ttoind(ii))))
end

%% SVM score Spearman correlation (rank similarity) across random subsets (per fixed proportion)
iprop = 3;
whichstat = 'medianpool';

rhoxneusublmlvsagg = cat(1, similagg.rhoxneusublmlvs);
rhoxneusublmlvstestagg = cat(1, rhoxneusublmlvsagg.testpair);
rhoxneusublmlvstestaggstat = cat(5, rhoxneusublmlvstestagg.(whichstat));
rhoxneusublmlvstestaggses = squeeze(mean(rhoxneusublmlvstestaggstat(1,iprop,:,:,:),3))';
rhoxneusubasisagg = cat(1, similagg.rhoxneusubasis);
rhoxneusubasistestagg = cat(1, rhoxneusubasisagg.testpair);
rhoxneusubasistestaggstat = cat(4, rhoxneusubasistestagg.(whichstat));
rhoxneusubasistestaggses = squeeze(mean(rhoxneusubasistestaggstat(1,iprop,:,:),3))';
rhoxneusublmlvssimilagg = cat(1, rhoxneusublmlvsagg.similpair);
rhoxneusublmlvssimilaggstat = cat(5, rhoxneusublmlvssimilagg.(whichstat));
rhoxneusublmlvssimilaggses = squeeze(rhoxneusublmlvssimilaggstat(1,iprop,:,:,:));
rhoxneusubasissimilagg = cat(1, rhoxneusubasisagg.similpair);
rhoxneusubasissimilaggstat = cat(5, rhoxneusubasissimilagg.(whichstat));
rhoxneusubasissimilaggses = squeeze(rhoxneusubasissimilaggstat(1,iprop,:,:,:));

figure
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', sprintf('randomly silence %.0f%% of neurons: SVM score consistency across neuronal subsets', 100*propneusilvec(iprop)), 'edgecolor', 'none')
subplot(3,3,1)
hold all
pl = plot(disperses,  rhoxneusublmlvstestaggses);
errorbar(disperses, mean(rhoxneusublmlvstestaggses,1), std(rhoxneusublmlvstestaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoxneusubasistestaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title('test trial pairs')

for ii = 1:numel(ttoind)
    rhoxneusublmlvsttaggses = (reshape(rhoxneusublmlvssimilaggses(ttoind(ii),:,:),numel(disperses),Nsessions))';
    rhoxneusubasisttaggses = (reshape(rhoxneusubasissimilaggses(ttoind(ii),:),1,Nsessions))';
subplot(3,3,1+ii)
hold all
pl = plot(disperses,  rhoxneusublmlvsttaggses);
errorbar(disperses, mean(rhoxneusublmlvsttaggses,1), std(rhoxneusublmlvsttaggses,0,1)/sqrt(Nsessions), 'ko-', 'LineWidth', 2)
xl = [disperses(1) disperses(end)];
plot(xl, mean(rhoxneusubasisttaggses)*[1 1], 'c-', 'LineWidth', 1)
xlabel('LMLV slopes')
ylabel('% rho(SVM score)==1')
title(sprintf('Trial%d pairs', hireptt(ttoind(ii))))
end


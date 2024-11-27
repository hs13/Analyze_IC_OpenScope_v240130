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

preproc = 'meancenter';
whichSVMkernel = 'Linear';
svmdesc = 'trainICRC';

excludeneuvar0 = 0;
optimizeSVM = 2;
fprintf('neuron exclusion criterion %d optimize option %d\n', excludeneuvar0, optimizeSVM)

% if 0, keep all neurons; if 1, exclude zero variance neurons in train trial
% types; if 2 exclude zero variance neurons in all trial types
switch excludeneuvar0
    case 0
        optexclneu = 'incl';
    case 1
        optexclneu = 'excltt';
    case 2
        optexclneu = 'excl';
    otherwise
        error('excludeneuvar0 option not recognized')
end

% optimizeSVM: 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
switch optimizeSVM
    case 0
        optoptim = '_nooptim';
    case 1
        optoptim = '_alloptim';
    case 2
        optoptim = '';
    case 3
        optoptim = '_onevsall';
    otherwise
        error('optimizeSVM option %d not recognized', SVMout.optimizeSVM)
end

consvmagg = struct();
for ises = 1:numel(nwbsessions)
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

    pltses = false;
    consvmfn = strcat(pathpp, 'consistency_SVM12_', svmdesc, optoptim, '_', preproc, '_lmlvslopes_', optexclneu, '.mat');

    consvm = load(consvmfn);
    if isempty(fieldnames(consvmagg))
        consvmagg = consvm;
    else
        consvmagg(ises) = consvm;
    end
end

disperses = consvm.disperses;
rhoXdivasisagg = cat(1, consvmagg.rhoXdivasis);
rhoXdivlmlvsagg = cat(1, consvmagg.rhoXdivlmlvs);

fs = 12;
%%
whichstat = 'prct';
switch whichstat
    case 'avg'
        ylab = 'mean rho(SVM score)';
    case 'prct'
        ylab = '% rho(SVM score)==1';
    case 'medianpool'
        ylab = 'median rho(SVM score)';
    otherwise
        error(['specify ' whichstat])
end

testt = [106,107,110,111];
hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];

figure('Position', [1000 100 1000 1000])
for isp = 1:4
    switch isp
        case 1
            whichrhoscorefield = 'test';
            ttoind = find(true(size(testt)));
        case 2
            whichrhoscorefield = 'simil';
            ttoind = find(~ismember(hireptt,testt));
        case 3
            whichrhoscorefield = 'testpair';
            ttoind = find(true(size(testt)));
        case 4
            whichrhoscorefield = 'similpair';
            ttoind = find(~ismember(hireptt,testt));
    end
    rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
    rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
    rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);

    rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
    rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
    rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);

    subplot(2,2,isp)
    hold all
    pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind,:),1)), squeeze(std(rhoXdivlmlvsaggstat(:,ttoind,:),0,1))/sqrt(Nsessions) );
    yl = ylim;
    for ii = 1:numel(pl)
        if ismember(hireptt(ttoind(ii)), [1105 1109])
            lw = 2;
            %disp(ii)
        else
            lw = 0.2;
        end
        errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind(ii),:),1)), squeeze(std(rhoXdivlmlvsaggstat(:,ttoind(ii),:),0,1))/sqrt(Nsessions), 'Color', pl(ii).Color, 'LineWidth', lw)
        h = squeeze(mean(rhoXdivasisaggstat(:,ttoind(ii) ),1));
        plot([disperses(1) disperses(end)], h*[1 1], 'Color', pl(ii).Color, 'LineWidth', lw)
        text(disperses(1), yl(1)+(numel(pl)-ii-1)*0.08*range(yl), num2str(hireptt(ttoind(ii))), 'Color', pl(ii).Color, 'FontSize', 10, 'VerticalAlignment', 'bottom')
    end
    % errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind,:),[1,2])), squeeze(std(mean(rhoXdivlmlvsaggstat(:,ttoind,:),2), 0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2 );
    ylim(yl)
    set(gca, 'FontSize', fs)
    xlabel('LMLV slopes')
    ylabel(ylab)
    title(whichrhoscorefield)
end


%%
ttoind = find(~ismember(hireptt, testt));

whichstat = 'prct';
switch whichstat
    case 'avg'
        ylab = 'mean rho(SVM score)';
    case 'prct'
        ylab = '% rho(SVM score)==1';
    case 'medianpool'
        ylab = 'median rho(SVM score)';
    otherwise
        error(['specify ' whichstat])
end

testt = [106,107,110,111];
hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];

figure('Position', [1000 100 1000 1000])
whichrhoscorefield = 'test';
rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);
rhoXdivlmlvsaggstatavg = squeeze(mean(rhoXdivlmlvsaggstat,2));

rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);
rhoXdivasisaggstatavg = squeeze(mean(rhoXdivasisaggstat,2));

subplot(3,3,1)
hold all
plot(disperses, rhoXdivlmlvsaggstatavg)
pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstatavg,1)), squeeze(std(rhoXdivlmlvsaggstatavg,0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2  );
h = mean(rhoXdivasisaggstatavg);
plot([disperses(1) disperses(end)], h*[1 1], 'c-', 'linewidth', 1)
sigslope = false(size(disperses));
for islope = 1:numel(disperses)
    p = signrank(squeeze(rhoXdivlmlvsaggstatavg(:,islope)), rhoXdivasisaggstatavg);
    sigslope(islope) = p<0.05;
end
scatter(disperses(sigslope), squeeze(mean(rhoXdivlmlvsaggstatavg(:,sigslope),1)), 'ro', 'filled' )
    set(gca, 'FontSize', fs)
xlabel('LMLV slopes')
ylabel(ylab)
title(whichrhoscorefield)

for ii = 1:numel(ttoind)
    whichrhoscorefield = 'simil';
    rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
    rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
    rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);

    rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
    rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
    rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);

    subplot(3,3,ii+1)
    hold all
    plot(disperses, squeeze(rhoXdivlmlvsaggstat(:,ttoind(ii),:)))
    pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind(ii),:),1)), squeeze(std(rhoXdivlmlvsaggstat(:,ttoind(ii),:),0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2  );
    h = mean(rhoXdivasisaggstat(:,ttoind(ii)) );
    plot([disperses(1) disperses(end)], h*[1 1], 'c-', 'linewidth', 1)
    % ylim(yl)
    sigslope = false(size(disperses));
    for islope = 1:numel(disperses)
        p = signrank(squeeze(rhoXdivlmlvsaggstat(:,ttoind(ii),islope)), rhoXdivasisaggstat(:,ttoind(ii)));
        sigslope(islope) = p<0.05;
    end
    scatter(disperses(sigslope), squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind(ii),sigslope),1)), 'ro', 'filled' )

    set(gca, 'FontSize', fs)
    xlabel('LMLV slopes')
    ylabel(ylab)
    title([whichrhoscorefield num2str(hireptt(ttoind(ii)))])
end

%%
ttoind = find(~ismember(hireptt, testt));

whichstat = 'prct';
switch whichstat
    case 'avg'
        ylab = 'mean rho(SVM score)';
    case 'prct'
        ylab = '% rho(SVM score)==1';
    case 'medianpool'
        ylab = 'median rho(SVM score)';
    otherwise
        error(['specify ' whichstat])
end

testt = [106,107,110,111];
hireptt = [0, 101, 105, 106, 107, 109, 110, 111, 1105, 1109, 1201, 1299];

figure('Position', [1000 100 1000 1000])
whichrhoscorefield = 'testpair';
rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);
rhoXdivlmlvsaggstatavg = squeeze(mean(rhoXdivlmlvsaggstat,2));

rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);
rhoXdivasisaggstatavg = squeeze(mean(rhoXdivasisaggstat,2));

subplot(3,3,1)
hold all
plot(disperses, rhoXdivlmlvsaggstatavg)
pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstatavg,1)), squeeze(std(rhoXdivlmlvsaggstatavg,0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2  );
h = mean(rhoXdivasisaggstatavg);
plot([disperses(1) disperses(end)], h*[1 1], 'c-', 'linewidth', 1)
sigslope = false(size(disperses));
for islope = 1:numel(disperses)
    p = signrank(squeeze(rhoXdivlmlvsaggstatavg(:,islope)), rhoXdivasisaggstatavg);
    sigslope(islope) = p<0.05;
end
scatter(disperses(sigslope), squeeze(mean(rhoXdivlmlvsaggstatavg(:,sigslope),1)), 'ro', 'filled' )
    set(gca, 'FontSize', fs)
xlabel('LMLV slopes')
ylabel(ylab)
title(whichrhoscorefield)

for ii = 1:numel(ttoind)
    whichrhoscorefield = 'similpair';
    rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
    rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
    rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);

    rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
    rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
    rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);

    subplot(3,3,ii+1)
    hold all
    plot(disperses, squeeze(rhoXdivlmlvsaggstat(:,ttoind(ii),:)))
    pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind(ii),:),1)), squeeze(std(rhoXdivlmlvsaggstat(:,ttoind(ii),:),0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2  );
    h = mean(rhoXdivasisaggstat(:,ttoind(ii)) );
    plot([disperses(1) disperses(end)], h*[1 1], 'c-', 'linewidth', 1)
    % ylim(yl)
    sigslope = false(size(disperses));
    for islope = 1:numel(disperses)
        p = signrank(squeeze(rhoXdivlmlvsaggstat(:,ttoind(ii),islope)), rhoXdivasisaggstat(:,ttoind(ii)));
        sigslope(islope) = p<0.05;
    end
    scatter(disperses(sigslope), squeeze(mean(rhoXdivlmlvsaggstat(:,ttoind(ii),sigslope),1)), 'ro', 'filled' )

    set(gca, 'FontSize', fs)
    xlabel('LMLV slopes')
    ylabel(ylab)
    title([whichrhoscorefield num2str(hireptt(ttoind(ii)))])
end

%%
whichrhoscorefield = 'simil'; % simil/similpair
whichstat = 'avg';
switch whichstat
    case 'avg'
        ylab = 'mean rho(SVM score)';
    case 'prct'
        ylab = '% rho(SVM score)==1';
    case 'medianpool'
        ylab = 'median rho(SVM score)';
    otherwise
        error(['specify ' whichstat])
end
infttoind = find(ismember(hireptt, [1105 1109]));

rhoXdivlmlvsaggperf = cat(1, rhoXdivlmlvsagg.(whichrhoscorefield));
rhoXdivlmlvsaggses = squeeze(mean( cat(4,rhoXdivlmlvsaggperf.(whichstat)) ,1));
rhoXdivlmlvsaggstat = permute( rhoXdivlmlvsaggses, [3,1,2]);
rhoXdivlmlvsaggstatavg = squeeze(mean(rhoXdivlmlvsaggstat(:,infttoind,:),2));

rhoXdivasisaggperf = cat(1, rhoXdivasisagg.(whichrhoscorefield));
rhoXdivasisaggses = squeeze(mean( cat(3,rhoXdivasisaggperf.(whichstat)) ,1));
rhoXdivasisaggstat = permute( rhoXdivasisaggses, [2,1]);
rhoXdivasisaggstatavg = squeeze(mean(rhoXdivasisaggstat(:,infttoind),2));

figure
hold all
plot(disperses, rhoXdivlmlvsaggstatavg)
pl = errorbar(disperses, squeeze(mean(rhoXdivlmlvsaggstatavg,1)), squeeze(std(rhoXdivlmlvsaggstatavg,0,1))/sqrt(Nsessions), 'k-', 'LineWidth',2  );
h = mean(rhoXdivasisaggstatavg );
plot([disperses(1) disperses(end)], h*[1 1], 'c-', 'linewidth', 1)
% ylim(yl)
sigslope = false(size(disperses));
for islope = 1:numel(disperses)
    p = signrank( rhoXdivlmlvsaggstatavg(:,islope), rhoXdivasisaggstatavg );
    sigslope(islope) = p<0.05;
end
scatter(disperses(sigslope), squeeze(mean(rhoXdivlmlvsaggstatavg(:,sigslope),1)), 'ro', 'filled' )

    set(gca, 'FontSize', fs)
xlabel('LMLV slopes')
ylabel(ylab)
title([whichrhoscorefield  ' TRE trials'])

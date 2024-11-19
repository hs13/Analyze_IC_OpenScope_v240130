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
Nsessions = numel(nwbsessions);

load([drivepath 'RESEARCH/ICexpts_revision23/openscope_psthavgall.mat'])
load([drivepath 'RESEARCH/ICexpts_revision23/openscope_popavg_all.mat'])

%% reactivcoeffagg
preprocs = {'meancenter', 'zscore'};
whichSVMkernel = 'Linear';
svmdesc = 'trainBK';
% optimizeSVM: 0 no optimization, 1 optimize hyperparameters, 2 onevsone, 3 onevsall
optimizeSVM = 2;
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
        error('optimizeSVM option %d not recognized', optimizeSVM)
end

neuV1RSagg = [];
reactivcoeffagg = struct();
for pp = 1:numel(preprocs)
    preproc = preprocs{pp};
reactivcoeffagg.(preproc) = [];
for ises = 1:numel(nwbsessions)
    mousedate = nwbsessions{ises};
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    disp(mousedate)
    whichsponind=2;
    svmsponfn = strcat(pathpp, 'SVMspon', num2str(whichsponind), '_invcv_', svmdesc, optoptim, '_V1_', whichSVMkernel, '_', preproc, '_incl.mat');
    load(svmsponfn, 'neuV1RS', 'reactivcoeff') % 'neuV1RS', 'Twin', 'Tslide', 'Tstartind', 'Tctr', 'SVMtrainBK', 'reactivcoeff'
    if pp==1
        neuV1RSagg = cat(1, neuV1RSagg, neuV1RS);
    end
    reactivcoeffagg.(preproc) = cat(1, reactivcoeffagg.(preproc), reactivcoeff);
end
end
neuV1RSagg = logical(neuV1RSagg);

load(svmsponfn, 'SVMtrainBK')
testt = SVMtrainBK.trialtypes;
Ntt = numel(testt);

%%
figure;
for itt = 1:Ntt
    subplot(2,3,itt)
    scatter(reactivcoeffagg.meancenter(:,itt), reactivcoeffagg.zscore(:,itt))
    title(testt(itt))
end

preproc = 'zscore';
ij=nchoosek(1:Ntt,2);
figure;
for isp = 1:size(ij,1)
    subplot(2,5,isp)
    scatter(reactivcoeffagg.(preproc)(:,ij(isp,1)), reactivcoeffagg.(preproc)(:,ij(isp,2)))
    xlabel(testt(ij(isp,1)))
    ylabel(testt(ij(isp,2)))
end

%%
preproc = 'zscore';
figure;
for itt = 1:Ntt
    switch testt(itt)
        case 0
            neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
            neuoj = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)<0.05;
            neuok = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)>=0.05;
        case 106
            neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc13(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RSagg);
        case 107
            neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc14(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RSagg);
        case 110
            neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc23(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RSagg);
        case 111
            neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc24(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RSagg);
        otherwise
            error('check itt')
    end
    reactivcoeffcell = cell(4,1);
    reactivcoeffcell{1} = reactivcoeffagg.(preproc)(:,itt);
    reactivcoeffcell{2} = reactivcoeffagg.(preproc)(neuoi==1,itt);
    reactivcoeffcell{3} = reactivcoeffagg.(preproc)(neuoj==1,itt);
    reactivcoeffcell{4} = reactivcoeffagg.(preproc)(neuok==1,itt);
    rccellind = cell(4,1);
    rccellind{1} = ones(size(reactivcoeffagg.(preproc),1),1);
    rccellind{2} = 2*ones(nnz(neuoi),1);
    rccellind{3} = 3*ones(nnz(neuoj),1);
    rccellind{4} = 4*ones(nnz(neuok),1);

    subplot(2,3,itt)
    hold all
    nhist(reactivcoeffcell)
    % h=histogram(reactivcoeffagg(:,itt));
    % histogram(reactivcoeffagg(neuoi==1,itt), 'BinEdges', h.BinEdges)
    % histogram(reactivcoeffagg(neuoj==1,itt), 'BinEdges', h.BinEdges)
    title(testt(itt))
end

[P,ANOVATAB,STATS] = kruskalwallis(cat(1,reactivcoeffcell{:}), cat(1,rccellind{:}));
figure; multcompare(STATS)

for itt = 1:Ntt
    switch testt(itt)
        case 0
            neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
            neuoj = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)<0.05;
            neuok = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)>=0.05;
        case 106
            neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RSagg);
            indenc1n3 = ICsigall.ICwcfg1_presentations.indenc1 | ICsigall.ICwcfg1_presentations.indenc3;
            neuoj = indenc1n3(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RSagg);
        case 107
            neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RSagg);
            indenc1n4 = ICsigall.ICwcfg1_presentations.indenc1 | ICsigall.ICwcfg1_presentations.indenc4;
            neuoj = indenc1n4(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RSagg);
        case 110
            neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RSagg);
            indenc2n3 = ICsigall.ICwcfg1_presentations.indenc2 | ICsigall.ICwcfg1_presentations.indenc3;
            neuoj = indenc2n3(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RSagg);
        case 111
            neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RSagg);
            indenc2n4 = ICsigall.ICwcfg1_presentations.indenc2 | ICsigall.ICwcfg1_presentations.indenc4;
            neuoj = indenc2n4(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RSagg);
        otherwise
            error('check itt')
    end
p = ranksum(reactivcoeffagg.(preproc)(neuoj==1,itt), reactivcoeffagg.(preproc)(neuok==1,itt));
pleft = ranksum(reactivcoeffagg.(preproc)(neuoj==1,itt), reactivcoeffagg.(preproc)(neuok==1,itt), 'tail', 'left');
fprintf('Trial%d p=%.4f pleft=%.4f\n', testt(itt), p, pleft)
end



%%
tempSP = ICsigall.ICwcfg1_presentations.SP_ICvsRC(neuV1RSagg);
figure;
for itt = 1:Ntt
    switch testt(itt)
        case 0
            neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
            neuoj = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)<0.05;
            neuok = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)>=0.05;
        case 106
            neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc13(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RSagg);
        case 107
            neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc14(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RSagg);
        case 110
            neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc23(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RSagg);
        case 111
            neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc24(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RSagg);
        otherwise
            error('check itt')
    end
    subplot(2,3,itt)
    hold all 
plot(tempSP, reactivcoeffagg.(preproc)(:,itt), 'k.')
plot(tempSP(neuoi), reactivcoeffagg.(preproc)(neuoi,itt), 'b.')
plot(tempSP(neuok), reactivcoeffagg.(preproc)(neuok,itt), 'go')
yl = ylim;
plot([0 1], [0 0], 'Color', 0.5*[1 1 1])
plot([0.5 0.5], yl, 'Color', 0.5*[1 1 1])
ylim(yl)
set(gca, 'XGRid', 'on', 'YGRid', 'on')
xlabel('IC vs RC AUROC')
ylabel('react. coeff.')
    title(testt(itt))
end

%%
tempSP = ICsigall.ICwcfg1_presentations.SP_ICvsRC(neuV1RSagg);

figure;
for itt = 1:Ntt
    switch testt(itt)
        case 0
            neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
            neuoj = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)<0.05;
            neuok = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)>=0.05;
        case 106
            neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc13(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RSagg);
        case 107
            neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc14(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RSagg);
        case 110
            neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc23(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RSagg);
        case 111
            neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc24(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RSagg);
        otherwise
            error('check itt')
    end
    neuhireact = reactivcoeffagg.(preproc)(:,itt)>=0.2;
    subplot(2,3,itt)
    hold all
    h=histogram(tempSP(neuoi==1));
    % h=histogram(tempSP);
    % histogram(tempSP(neuoi==1), 'BinEdges', h.BinEdges)
    histogram(tempSP(neuoi==1 & neuhireact), 'BinEdges', h.BinEdges)
    xlabel('IC vs RC AUROC')

p = ranksum(tempSP(neuoi==1), tempSP(neuoi==1 & neuhireact));
pleft = ranksum(tempSP(neuoi==1), tempSP(neuoi==1 & neuhireact), 'tail', 'left');
    title(sprintf('Trial%d pleft=%.4f', testt(itt), pleft))
end

%% scatterplot (mean(RIC)-mean(RLC))/std(Rblank) vs reactivation coefficient
Rgestalt = ( mean(Ronavgall.ICwcfg1_presentations(:,ismember(ICtrialtypes,[106 111])),2) ...
    - mean(Ronavgall.ICwcfg1_presentations(:,ismember(ICtrialtypes,[107 110])),2) ) ...
    ./ Ronstdall.ICwcfg1_presentations(:,ICtrialtypes==0);
Rgestalt = Rgestalt(neuV1RSagg);

figure;
for itt = 1:Ntt
    switch testt(itt)
        case 0
            neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
            neuoj = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)<0.05;
            neuok = neuoi & ICsigall.ICwcfg1_presentations.PkwBI(neuV1RSagg)>=0.05;
        case 106
            neuoi = ICsigall.ICwcfg1_presentations.ICresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc13(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder1(neuV1RSagg);
        case 107
            neuoi = ICsigall.ICwcfg1_presentations.RCresp1(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc14(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder1(neuV1RSagg);
        case 110
            neuoi = ICsigall.ICwcfg1_presentations.RCresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc23(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.RCencoder2(neuV1RSagg);
        case 111
            neuoi = ICsigall.ICwcfg1_presentations.ICresp2(neuV1RSagg);
            neuoj = ICsigall.ICwcfg1_presentations.indenc24(neuV1RSagg);
            neuok = ICsigall.ICwcfg1_presentations.ICencoder2(neuV1RSagg);
        otherwise
            error('check itt')
    end
    subplot(2,3,itt)
    hold all 
plot(Rgestalt, reactivcoeffagg.(preproc)(:,itt), 'k.')
plot(Rgestalt(neuoi), reactivcoeffagg.(preproc)(neuoi,itt), 'b.')
% plot(Rgestalt(neuok), reactivcoeffagg.(preproc)(neuok,itt), 'go')
yl = ylim;
plot([0 1], [0 0], 'Color', 0.5*[1 1 1])
plot([0.5 0.5], yl, 'Color', 0.5*[1 1 1])
ylim(yl)
set(gca, 'XGRid', 'on', 'YGRid', 'on')
xlabel('(mean(RIC)-mean(RLC))/std(Rblank)')
ylabel('react. coeff.')
    title(testt(itt))
end

disp( [testt' prctile(reactivcoeffagg.(preproc),[25 50 75],1)'] )
disp( [testt' prctile(reactivcoeffagg.(preproc)(Rgestalt>0.5,:),[25 50 75],1)'] )

itt = testt==106;
neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
%neuoi = true(nnz(neuV1RSagg),1 );
ranksum(reactivcoeffagg.(preproc)(neuoi,itt), reactivcoeffagg.(preproc)(neuoi & Rgestalt>0.5,itt))


itt = testt==106;
neuoi = ICsigall.ICwcfg1_presentations.PkwBK(neuV1RSagg)<0.05;
ranksum(Rgestalt(neuoi), Rgestalt(neuoi & reactivcoeffagg.(preproc)(:,itt)>0.2))
figure; hold all
h = histogram(Rgestalt(neuoi), 'normalization', 'probability');
histogram(Rgestalt(neuoi & reactivcoeffagg.(preproc)(:,itt)>0.2), h.BinEdges, 'normalization', 'probability')

%%
kerwinhalf = 12; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

temppsth = psthavgall.ICwcfg1_presentations(:,:,neuV1RSagg);
temppsthavg = convn(temppsth, kergauss, 'same');
neu2plt = find(neuV1RSagg);
sesneu2plt = sesneuall(neuV1RSagg);

tt2p = [106   107   110   111];
ttcol = [0 .4 0; .5 0.25 0; 1 0.5 0; 0 1 0];
Ntop = 10;
figure('Position', [0 0 2000 1500])
for itt = 1:Ntt
    [sv,si]=sort(reactivcoeffagg.(preproc)(:,itt), 'descend', 'MissingPlacement', 'last');
    for ii = 1:Ntop
        ci = si(ii);
        subplot(Ntt,Ntop,Ntop*(itt-1)+ii)
        hold all
        for typi = 1:numel(tt2p)
            plot(psthtli, temppsthavg(:,ICtrialtypes==tt2p(typi),ci), 'Color', ttcol(typi,:), 'LineWidth', 1)
        end
        yl =ylim;
        plot([0 0], yl, '-', 'Color', 0.5*[1 1 1])
        plot(400+[0 0], yl, '-', 'Color', 0.5*[1 1 1])
        ylim(yl)
        xlim([-200 600])
        title(sprintf('Trial%d Ses%d Cell%d\nreact. coeff. %.4f', testt(itt), sesneu2plt(ci), neu2plt(ci), reactivcoeffagg.(preproc)(ci,itt) ))
    end
end

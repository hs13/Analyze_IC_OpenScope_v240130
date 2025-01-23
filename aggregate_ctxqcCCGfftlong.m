addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
% nwbdir = dir(datadir);
% nwbsessions = {nwbdir.name};
% nwbsessions = nwbsessions(( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', 'sub-631510', ...
    'sub-631570', 'sub-633229', 'sub-637484'};
Nsessions = numel(nwbsessions);
ses2anal = 1:Nsessions;

%{
ccgpath = 'S:\OpenScopeData\00248_v240130\CCG\';
for ises = ses2anal
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    clearvars -except ises ses2anal Nsessions nwbsessions ccgpath
    % pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathccg = [ccgpath nwbsessions{ises} filesep];
    load([pathccg 'ctxCCGfftnew.mat'])

    tli25 = CCGtli_fft>=-12 & CCGtli_fft<=12;
    CCGtli = CCGtli_fft(tli25);
    ctxCCG = ctxCCGfftnew(:,:,tli25);

    ctxCCGm0 = mean(ctxCCG,3);
    ctxCCGmc = ctxCCG - ctxCCGm0;

    tic
    [excpeak, inhtrough] = computeCCGxconn(CCGtli, ctxCCGmc, 10);
    toc
    save([pathccg 'ctxCCGjcxconn.mat'], 'excpeak', 'inhtrough')
end
%}
%%
% A sharp peak was deemed significant if the maximum of jitter-corrected CCG 
% amplitude within a ±10 ms window had a magnitude larger than sevenfold of the 
% standard deviation of the CCG flanks (between ±50–100 ms from zero). 
% All subsequent analysis was based on significant CCG sharp peaks.

ccgpath = 'S:\OpenScopeData\00248_v240130\CCG\';

% 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli', 'ctrctxCCG', 'ctrctxCCGweight'
neuctxagg = [];
neulocctxagg = [];
Nneurons_visarea = NaN(Nsessions,6);
visarealabelagg = [];
sesctxagg = [];
ctxCCGjcpostpkagg = cell(Nsessions,1);
ctxCCGjcpostTpkagg = cell(Nsessions,1);
ctxCCGjcprepkagg = cell(Nsessions,1);
ctxCCGjcpreTpkagg = cell(Nsessions,1);
ctxCCGflankstdagg = cell(Nsessions,1);

tic
for ises = ses2anal
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    clearvars ctxCCGfft ctxCCGweight visarealabels neuctx neulocctx

    pathccg = [ccgpath nwbsessions{ises} filesep];
    load([pathccg 'ctxCCGjcxconn.mat'])
    ctxCCGjcpostpkagg{ises} = squeeze(excpeak(:,:,1));
    ctxCCGjcpostTpkagg{ises} = squeeze(excpeak(:,:,2));

    ctxCCGjcprepkagg{ises} = squeeze(excpeak(:,:,1))';
    ctxCCGjcpreTpkagg{ises} = squeeze(excpeak(:,:,2))';

    load([pathccg 'ctxCCGfftnew.mat'])
    [v,c]=uniquecnt(visarealabels);
    Nneurons_visarea(ises,ismember(1:6, v)) = c(ismember(v,1:6));
    neuctxagg = cat(1, neuctxagg, neuctx);
    neulocctxagg = cat(1, neulocctxagg, neulocctx);
    visarealabelagg = cat(1, visarealabelagg, visarealabels);
    sesctxagg = cat(1, sesctxagg, ises*ones(size(visarealabels)) );

    flanktli = abs(CCGtli_fft)>=50;
    ctxCCGflankstdagg{ises} = std(ctxCCGfftnew(:,:,flanktli),0,3);

    toc
end

save('G:\My Drive\RESEARCH\ICexpts_revision23\openscope_ctxCCGfftagg.mat', ...
    'Nneurons_visarea', 'neuctxagg', 'neulocctxagg', 'visarealabelagg', 'sesctxagg', ...
    'ctxCCGjcpostpkagg', 'ctxCCGjcpostTpkagg', 'ctxCCGjcprepkagg', 'ctxCCGjcpreTpkagg', 'ctxCCGflankstdagg', '-v7.3')

%%
ccgpath = 'S:\OpenScopeData\00248_v240130\CCG\';

% 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli', 'ctrctxCCG', 'ctrctxCCGweight'
neuctxagg = [];
neulocctxagg = [];
Nneurons_visarea = NaN(Nsessions,6);
visarealabelagg = [];
sesctxagg = [];
ctxCCGmean4msagg = cell(Nsessions,1);

tic
for ises = ses2anal
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    clearvars ctxCCGfft ctxCCGweight visarealabels neuctx neulocctx

    pathccg = [ccgpath nwbsessions{ises} filesep];
    load([pathccg 'ctxCCGfftnew.mat'])
    [v,c]=uniquecnt(visarealabels);
    Nneurons_visarea(ises,ismember(1:6, v)) = c(ismember(v,1:6));
    neuctxagg = cat(1, neuctxagg, neuctx);
    neulocctxagg = cat(1, neulocctxagg, neulocctx);
    visarealabelagg = cat(1, visarealabelagg, visarealabels);
    sesctxagg = cat(1, sesctxagg, ises*ones(size(visarealabels)) );

    tli4ms = CCGtli_fft>0 & CCGtli_fft<=4;
    ctxCCGmean4msagg{ises} = mean(ctxCCGfftnew(:,:,tli4ms),3);

    toc
end


ctxCCGflankstd = std(ctxCCGfftnew(:,:,abs(CCGtli_fft)>=50),0,3);
ctxCCGflankstdpost = std(ctxCCGfftnew(:,:,CCGtli_fft>=50),0,3);
ctxCCGflankstdpre = std(ctxCCGfftnew(:,:,CCGtli_fft<=-50),0,3);
figure; 
for isp = 1:4
    switch isp
        case 1
            tempx = repmat(spkcntvec,1,length(spkcntvec));
            tempy = ctxCCGflankstdpost;
        case 2
            tempx = repmat(spkcntvec,1,length(spkcntvec))';
            tempy = ctxCCGflankstdpost;
        case 3
            tempx = repmat(spkcntvec,1,length(spkcntvec));
            tempy = ctxCCGflankstd;
        case 4
            tempx = repmat(spkcntvec,1,length(spkcntvec))';
            tempy = ctxCCGflankstd;
    end
subplot(2,2,isp)
plot(tempx, tempy, 'o')
tempR = corr(tempx(:),tempy(:), 'rows', 'complete');
title( sprintf('R=%.4f', tempR))
end


save('G:\My Drive\RESEARCH\ICexpts_revision23\openscope_ctxCCGmean4ms_geonormagg.mat', ...
    'Nneurons_visarea', 'neuctxagg', 'neulocctxagg', 'visarealabelagg', 'sesctxagg', 'ctxCCGmean4msagg', '-v7.3')

%% report number of neurons
fprintf('V1 neurons mean %.4f sem %.4f\n', nanmean(Nneurons_visarea(:,1)), ...
    nanstd(Nneurons_visarea(:,1))/sqrt(size(Nneurons_visarea,1)) )

fprintf('HVA neurons mean %.4f sem %.4f\n', nanmean(sum(Nneurons_visarea(:,2:6),2)), ...
nanstd(sum(Nneurons_visarea(:,2:6),2))/sqrt(size(Nneurons_visarea,1)) )

%% sanity check
tli25 = CCGtli_fft>=-12 & CCGtli_fft<=12;
CCGtli = CCGtli_fft(tli25);
ctxCCG = ctxCCGfftnew(:,:,tli25);

ctxCCGm0 = mean(ctxCCG,3);
ctxCCGmc = ctxCCG - ctxCCGm0;

posttli = CCGtli(CCGtli>0 & CCGtli<=10);
[postmv,postmi]= max(ctxCCGmc(:,:,CCGtli>0 & CCGtli<=10), [],3);

mean(postmv==squeeze(excpeak(:,:,1)), 'all')
%nnz(isnan(postmv) & ~isnan(excpk))

pretli = CCGtli(CCGtli<0 & CCGtli>=-10);
[premv,premi]= max(ctxCCGmc(:,:,CCGtli<0 & CCGtli>=-10), [],3);

temppost = postmv;
temppost(eye(size(temppost))==1)= NaN;
temppre = premv;
temppre(eye(size(temppre))==1)= NaN;

isequaln(temppost, temppre')

mv = max(abs((temppost-temppre')./temppost),[],'all');
[r,c]=find(abs((temppost-temppre')./temppost)==mv, 1,'first');

figure; hold all
plot(CCGtli, squeeze(ctxCCGmc(c,r,:)), 'k-')
plot(CCGtli, flip(squeeze(ctxCCGmc(r,c,:))), 'r--')
plot(posttli(postmi(c,r)), postmv(c,r), 'ko')
plot(pretli(premi(c,r)), premv(c,r), 'ro')

c=100; r= 161;
figure; hold all
plot(CCGtli_fft, squeeze(ctxCCGfftnew(c,r,:)), 'k-')
plot(CCGtli_fft, flip(squeeze(ctxCCGfftnew(r,c,:))), 'r--')

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
load([datadir 'postprocessed\openscope_popavg_all.mat'])
% 'nwbsessions', 'neuoindall', 'probeneuall', 'neulocall', 'neupeakchall', 'sesneuall', 'neuctxall', ...
% 'meanFRvecall', 'sponFRvecall', 'vis', ...
% 'unit_amplitude_all', 'unit_isi_violations_all', 'unit_wfdur_all', 'unit_amplitude_cutoff_all', 'unit_presence_ratio_all', ...
% 'ICsigfields', 'ICsigall', 'RFCIfields', 'RFCIall', ...
% 'RFCIspinfields', 'RFCIspinall', 'sizeCIfields', 'sizeCIall', ...
% 'oriparamsfields', 'oriparamsall', 'ori4paramsfields', 'ori4paramsall', '-v7.3')

%% Siegle et al's single unit filter criteria: isi_violations < 0.5 & amplitude_cutoff < 0.1 & presence_ratio > 0.9
neuV1 = contains(neulocall, 'VISp') & ~contains(neulocall, 'VISpm');
neuRS = unit_wfdur_all>0.4;
neufilt = (unit_isi_violations_all<0.5 & unit_amplitude_cutoff_all<0.5 & unit_presence_ratio_all>0.9);

neuinarea = neuV1 & neuRS;% & neufilt;

neuqcrs = neuRS;% & neufilt;
neuctxqcrs = neuqcrs(neuctxall);

if isequal(sesctxagg, sesneuall(neuctxall)) && isequal(neulocctxagg, neulocall(neuctxall))
else
    error('check correspondence between openscope_ctxCCGfftagg and openscope_popavg_all')
end

%% USE FLANK STANDARD DEVIATION AS THRESHOLD
thrflankstd = 7;
neuctx2incl = neuctxqcrs & ICsigall.ICwcfg1_presentations.PkwBK(neuctxall)<0.05;

Nctxqcrs = nnz(neuctxqcrs);
areas = {'all', 'V1', 'HVA', 'LM', 'sigkwBK', 'V1sigkwBK', 'HVAsigkwBK', 'LMsigkwBK'};
xflankstdthr = 0:0.5:20;
Pout_flank = struct();
Pin_flank = struct();
for a = 1:numel(areas)
    Pout_flank.(areas{a}) = NaN(Nctxqcrs,length(xflankstdthr));
    Pin_flank.(areas{a}) = NaN(Nctxqcrs,length(xflankstdthr));
end
for ises = ses2anal
    neuses = sesctxagg(neuctxqcrs)==ises;
    noteye = ~eye(size(ctxCCGjcpostpkagg{ises}));
    for a = 1:numel(areas)
        switch areas{a}
            case 'all'
                sesneuinarea = neuctxqcrs(sesctxagg==ises);
            case 'V1'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==1 & neuctxqcrs(sesctxagg==ises);
            case 'HVA'
                sesneuinarea = visarealabelagg(sesctxagg==ises)>1 & neuctxqcrs(sesctxagg==ises);
            case 'LM'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==2 & neuctxqcrs(sesctxagg==ises);
            case 'sigkwBK'
                sesneuinarea = neuctx2incl(sesctxagg==ises);
            case 'V1sigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==1 & neuctx2incl(sesctxagg==ises);
            case 'HVAsigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)>1 & neuctx2incl(sesctxagg==ises);
            case 'LMsigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==2 & neuctx2incl(sesctxagg==ises);
        end
        for ix = 1:length(xflankstdthr)
            tempCCG = ctxCCGjcpostpkagg{ises} > xflankstdthr(ix)*ctxCCGflankstdagg{ises};
            Pout_flank.(areas{a})(neuses,ix) = sum(tempCCG(neuctxqcrs(sesctxagg==ises), sesneuinarea),2)./sum(noteye(neuctxqcrs(sesctxagg==ises), sesneuinarea),2) ;
            
            tempCCG = ctxCCGjcprepkagg{ises} > xflankstdthr(ix)*ctxCCGflankstdagg{ises};
            Pin_flank.(areas{a})(neuses,ix) = sum(tempCCG(neuctxqcrs(sesctxagg==ises), sesneuinarea),2)./sum(noteye(neuctxqcrs(sesctxagg==ises), sesneuinarea),2) ;
        end
    end
end

%% firing rate confound...
figure
subplot(2,2,1)
plot(meanFRvecall(neuqcrs & neuctxagg==1), Pin_flank.HVA(:,xflankstdthr==thrflankstd), '.')
xlabel('Firing Rate')
ylabel('Input from HVA')

subplot(2,2,2)
plot(meanFRvecall(neuqcrs & neuctxagg==1), Pin_flank.V1(:,xflankstdthr==thrflankstd), '.')
xlabel('Firing Rate')
ylabel('Input from V1')

subplot(2,2,3)
histogram2(meanFRvecall(neuqcrs & neuctxagg==1), Pin_flank.HVA(:,xflankstdthr==thrflankstd), 'DisplayStyle', 'tile')
xlabel('Firing Rate')
ylabel('Input from HVA')

subplot(2,2,4)
histogram2(meanFRvecall(neuqcrs & neuctxagg==1), Pin_flank.V1(:,xflankstdthr==thrflankstd), 'DisplayStyle', 'tile')
xlabel('Firing Rate')
ylabel('Input from V1')

%% CCG input/output # of connections relationship

ICenc = ICsigall.ICwcfg1_presentations.ICencoder==1;
indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin2==1 ...
    | ICsigall.ICwcfg1_presentations.indin3==1 | ICsigall.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigall.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
neuctrl = indin(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
% neuctrl = neuctx2incl & visarealabelagg==1;

% positive relationship
figure
hold all
plot(Pin_flank.V1(:,xflankstdthr==thrflankstd), Pout_flank.V1(:,xflankstdthr==thrflankstd), '.', 'Color', [0.5 0.5 0.5])
plot(Pin_flank.V1(neuctrl,xflankstdthr==thrflankstd), Pout_flank.V1(neuctrl,xflankstdthr==thrflankstd), 'x', 'Color', [0.5 0 1])
plot(Pin_flank.V1(neuoi,xflankstdthr==thrflankstd), Pout_flank.V1(neuoi,xflankstdthr==thrflankstd), 'o', 'Color', [0 0.7 0], 'linewidth', 1)
xlabel('V1 inputs')
ylabel('V1 outputs')

figure
hold all
plot(Pin_flank.V1(:,xflankstdthr==thrflankstd), Pin_flank.HVA(:,xflankstdthr==thrflankstd), '.', 'Color', [0.5 0.5 0.5])
plot(Pin_flank.V1(neuctrl,xflankstdthr==thrflankstd), Pin_flank.HVA(neuctrl,xflankstdthr==thrflankstd), 'x', 'Color', [0.5 0 1])
plot(Pin_flank.V1(neuoi,xflankstdthr==thrflankstd), Pin_flank.HVA(neuoi,xflankstdthr==thrflankstd), 'o', 'Color', [0 0.7 0], 'linewidth', 1)
xlabel('V1 inputs')
ylabel('HVA inputs')

%%
fs = 18;

ICenc = ICsigall.ICwcfg1_presentations.ICencoder==1;
indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin2==1 ...
    | ICsigall.ICwcfg1_presentations.indin3==1 | ICsigall.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigall.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
neuctrl = indin(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
% neuctrl = neuctx2incl & visarealabelagg==1;

V1neugroups = {'ICenc', 'indin'};
Pout_flank_med = struct();
Pin_flank_med = struct();
Pout_flank_uq = struct();
Pin_flank_uq = struct();
Pout_flank_lq = struct();
Pin_flank_lq = struct();
Pout_flank_avg = struct();
Pin_flank_avg = struct();
Pout_flank_sem = struct();
Pin_flank_sem = struct();
for n = 1:2
    switch V1neugroups{n}
        case 'ICenc'
            tempneu =  ICenc(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
        case 'indin'
            tempneu = indin(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
    end
for a = 1:numel(areas)
    Pout_flank_med.(areas{a}).(V1neugroups{n}) = median(Pout_flank.(areas{a})(tempneu,:),1);
    Pin_flank_med.(areas{a}).(V1neugroups{n}) = median(Pin_flank.(areas{a})(tempneu,:),1);
    Pout_flank_uq.(areas{a}).(V1neugroups{n}) = prctile(Pout_flank.(areas{a})(tempneu,:),75,1);
    Pin_flank_uq.(areas{a}).(V1neugroups{n}) = prctile(Pin_flank.(areas{a})(tempneu,:),75,1);
    Pout_flank_lq.(areas{a}).(V1neugroups{n}) = prctile(Pout_flank.(areas{a})(tempneu,:),25,1);
    Pin_flank_lq.(areas{a}).(V1neugroups{n}) = prctile(Pin_flank.(areas{a})(tempneu,:),25,1);
    Pout_flank_avg.(areas{a}).(V1neugroups{n}) = mean(Pout_flank.(areas{a})(tempneu,:),1);
    Pin_flank_avg.(areas{a}).(V1neugroups{n}) = mean(Pin_flank.(areas{a})(tempneu,:),1);
    Pout_flank_sem.(areas{a}).(V1neugroups{n}) = std(Pout_flank.(areas{a})(tempneu,:),0,1)/sqrt(nnz(tempneu));
    Pin_flank_sem.(areas{a}).(V1neugroups{n}) = std(Pin_flank.(areas{a})(tempneu,:),0,1)/sqrt(nnz(tempneu));
end
end

%%
ICenc = ICsigall.ICwcfg1_presentations.ICencoder==1;
indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin2==1 ...
    | ICsigall.ICwcfg1_presentations.indin3==1 | ICsigall.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigall.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
neuctrl = indin(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
% neuctrl = neuctx2incl & visarealabelagg==1;
V1neugroups = {'ICenc', 'indin'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 300 240]);
hold all
for n = 1:2
errorbar(xflankstdthr, Pin_flank_avg.HVA.(V1neugroups{n}), Pin_flank_sem.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
end
xlim([0 8])
legend({'IC-enc.', 'Seg.'})
legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Factor of Flank Std', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)

tempP = Pin_flank.HVA;
tempPcat = cat(1, tempP(neuoi,:), tempP(neuctrl,:));
tempgroup = [ones(nnz(neuoi),1); 0.5+ones(nnz(neuctrl),1)];

figure; 
a = boxchart(tempP(neuoi,:));
hold on
b = boxchart(tempP(neuctrl,:));

% figure; 
% boxchart( reshape(repmat(xsm0thr,length(tempgroup),1),[],1), reshape(tempPin,[],1), 'GroupByColor', reshape(repmat(tempgroup,1,length(xsm0thr)),[],1));

figure; 
b = boxchart( reshape(repmat(1:length(xflankstdthr),length(tempgroup),1),[],1), reshape(tempPcat,[],1), 'GroupByColor', reshape(repmat(tempgroup,1,length(xflankstdthr)),[],1));
b(1).BoxFaceColor = [0 0.7 0];
b(1).MarkerColor = [0 0.7 0];
b(2).BoxFaceColor = [0.42 0 0.42];
b(2).MarkerColor = [0.42 0 0.42];
set(gca, 'XTick', 1:length(xflankstdthr), 'XTickLabel', xflankstdthr)

%%
yl = [0 1];
xtl = {'IC-enc.', 'Seg. Resp.'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 1:2
errorbar(xflankstdthr, Pin_flank_avg.HVA.(V1neugroups{n}), Pin_flank_sem.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
text(4, yl(2)-(n-1)*0.1*range(yl), xtl{n}, 'Color', tempCData(n,:), 'VerticalAlignment', 'top', 'FontSize', fs)
end
plot(thrflankstd*[1 1], yl, 'k--')
ylim(yl)
xlim([0 8])
% legend({'IC-enc.', 'Seg.'})
% legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 1:2
errorbar(xflankstdthr, Pin_flank_avg.V1.(V1neugroups{n}), Pin_flank_sem.V1.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
text(4, yl(2)-(n-1)*0.1*range(yl), xtl{n}, 'Color', tempCData(n,:), 'VerticalAlignment', 'top', 'FontSize', fs)
end
plot(thrflankstd*[1 1], yl, 'k--')
ylim(yl)
xlim([0 8])
% legend({'IC-enc.', 'Seg.'})
% legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Input from V1', 'FontSize', fs)

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 1:2
errorbar(xflankstdthr, Pout_flank_avg.V1.(V1neugroups{n}), Pout_flank_sem.V1.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
text(4, yl(2)-(n-1)*0.1*range(yl), xtl{n}, 'Color', tempCData(n,:), 'VerticalAlignment', 'top', 'FontSize', fs)
end
plot(thrflankstd*[1 1], yl, 'k--')
ylim(yl)
xlim([0 8])
% legend({'IC-enc.', 'Seg.'})
% legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Output in V1', 'FontSize', fs)

%{
tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 2:-1:1
errorbar(xsm0thr, Pin_med.HVA.(V1neugroups{n}), Pin_med.HVA.(V1neugroups{n})-Pin_lq.HVA.(V1neugroups{n}), Pin_uq.HVA.(V1neugroups{n})-Pin_med.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
end
xlim([0 8])
legend({'IC-enc.', 'Seg.'})
legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)
%}

%% bar plot 
fs = 12;

ICenc = ICsigall.ICwcfg1_presentations.ICencoder==1;
indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin2==1 ...
    | ICsigall.ICwcfg1_presentations.indin3==1 | ICsigall.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigall.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;
neuctrl = indin(neuqcrs & neuctxagg==1) & visarealabelagg(neuctxqcrs)==1;

for ifg = 1:6%[1 3]
switch ifg
    case 1
tempP = Pin_flank.LM(:,xflankstdthr==thrflankstd);
ylab = 'Input from LM';
    case 2
tempP = Pout_flank.LM(:,xflankstdthr==thrflankstd);
ylab = 'Output to LM';
    case 3
tempP = Pin_flank.HVA(:,xflankstdthr==thrflankstd);
ylab = 'Input from HVA';
    case 4
tempP = Pout_flank.HVA(:,xflankstdthr==thrflankstd);
ylab = 'Output to HVA';
    case 5
tempP = Pin_flank.V1(:,xflankstdthr==thrflankstd);
ylab = 'Input in V1';
    case 6
tempP = Pout_flank.V1(:,xflankstdthr==thrflankstd);
ylab = 'Output in V1';
end
yl = [0 0.02];
tempPcat = cat(1, tempP(neuoi,:), tempP(neuctrl,:));
tempgroup = [ones(nnz(neuoi),1); 2*ones(nnz(neuctrl),1)];

[X,Y,T,AUC] =perfcurve(tempgroup, tempPcat, 1, 'NBoot', 1000);
disp(AUC)

xtl = {'IC-enc.', 'Seg.'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [200*ifg-100 100 180 240])
hold all
for ii = 1:2
b = boxchart(tempgroup(tempgroup==ii), tempPcat(tempgroup==ii), 'notch' , 'on', 'linewidth', 2);%, 'BoxFaceColor', tempCData(ii,:)); %, 'FontName', 'Arial')
b.BoxFaceColor = tempCData(ii,:);
b.MarkerColor = 0.15+tempCData(ii,:);
end
% yl = ylim;
p=ranksum(tempP(neuoi,:), tempP(neuctrl,:));
if p<0.05
    scatter(1.5, yl(2), 50, 'k*', 'LineWidth', 1)
end
set(gca, 'FontSize', fs, 'XTickLabelRotation', 0, 'XTick', 1:2, ...
    'XTickLabel', xtl);%, 'YTick', yl(1):ytd:yl(2))%, 'YTickLabelRotation', 45)
xlim([0.5 2.5])
ylim(yl)
xlabel('V1 Neurons', 'FontSize', fs)
ylabel(ylab, 'FontSize', fs)
title(sprintf('p=%.4f', p), 'FontSize', fs, 'Fontweight', 'normal')
% annotation('textbox', [0 0.92 1 0.1], 'string', figtit, 'fontsize', fs, 'edgecolor', 'none', 'horizontalalignment', 'center')
end

figure
hold all
h = histogram(tempP(neuctrl,:));
histogram(tempP(neuoi,:), 'BinEdges', h.BinEdges);

%%
figure; boxplot(degconverge, visarealabelagg)


tempP = Pout_flank.V1(:,xflankstdthr==thrflankstd);
ylab = 'Output in V1';
tempPcat = cat(1, tempP(neuoi,:), tempP(neuctrl,:));
tempgroup = [ones(nnz(neuoi),1); 2*ones(nnz(neuctrl),1)];
[X,Y,T,AUC] =perfcurve(tempgroup, tempPcat, 1, 'NBoot', 1000);


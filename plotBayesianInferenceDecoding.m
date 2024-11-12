if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])

Nsessions = numel(nwbsessions);
Ntt = 4;
Ninfertt=2;
%%
pathpp = 'S:\OpenScopeData\00248_v240130\postprocessed\sub-620333\';
load(strcat(pathpp, 'bayesinferencedecoding_V1.mat'))
load(strcat(pathpp, 'bayesinferencedecoding_V1_lmlvslopes.mat'))
infdecoders = who('-file', strcat(pathpp, 'bayesinferencedecoding_V1.mat'));
asiscmagg = struct();
lmlvscmagg = struct();
asisperfagg = struct();
lmlvsperfagg = struct();
for d = 1:numel(infdecoders)
    asiscmagg.(infdecoders{d}).train = NaN(Ntt, Ntt, Nsessions);
    asiscmagg.(infdecoders{d}).test = NaN(Ntt, Ntt, Nsessions);
    asiscmagg.(infdecoders{d}).inference = NaN(Ninfertt, Ntt, Nsessions);
    lmlvscmagg.(infdecoders{d}).train = NaN(Ntt, Ntt, Nsessions, length(disperses));
    lmlvscmagg.(infdecoders{d}).test = NaN(Ntt, Ntt, Nsessions, length(disperses));
    lmlvscmagg.(infdecoders{d}).inference = NaN(Ninfertt, Ntt, Nsessions, length(disperses));

    asisperfagg.(infdecoders{d}).train = NaN(Nsessions,1);
    asisperfagg.(infdecoders{d}).test = NaN(Nsessions,1);
    asisperfagg.(infdecoders{d}).inference = NaN(Nsessions,1);
    lmlvsperfagg.(infdecoders{d}).train = NaN(Nsessions, length(disperses));
    lmlvsperfagg.(infdecoders{d}).test = NaN(Nsessions, length(disperses));
    lmlvsperfagg.(infdecoders{d}).inference = NaN(Nsessions, length(disperses));
end

for ises = 1:numel(nwbsessions)
mousedate = nwbsessions{ises};
pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];

load(strcat(pathpp, 'bayesinferencedecoding_V1.mat'))
load(strcat(pathpp, 'bayesinferencedecoding_V1_lmlvslopes.mat'))

infdecoders = who('-file', strcat(pathpp, 'bayesinferencedecoding_V1.mat'));
asiscm = struct();
lmlvscm = struct();
asisperf = struct();
lmlvsperf = struct();
for d = 1:numel(infdecoders)
    tempdecoder = eval(infdecoders{d});
    tempdecoderlmlv = eval([infdecoders{d} '_lmlvs']);

    asiscm.(infdecoders{d}).train = squeeze(mean(tempdecoder.trainacc,3));
    asiscm.(infdecoders{d}).test = squeeze(mean(tempdecoder.testacc,3));
    asiscm.(infdecoders{d}).inference = squeeze(mean(tempdecoder.infperf,3));
    
    asisperf.(infdecoders{d}).train = mean(diag(mean(tempdecoder.trainacc,3)));
    asisperf.(infdecoders{d}).test = mean(diag(mean(tempdecoder.testacc,3)));
    tempinfperf = squeeze(mean(tempdecoder.infperf,3));
    asisperf.(infdecoders{d}).inference = ( tempinfperf(1,1)-tempinfperf(1,2)-tempinfperf(2,3)+tempinfperf(2,4) )/2;

    lmlvscm.(infdecoders{d}).train = NaN(Ntt, Ntt, length(disperses));
    lmlvscm.(infdecoders{d}).test = NaN(Ntt, Ntt, length(disperses));
    lmlvscm.(infdecoders{d}).inference = NaN(Ninfertt, Ntt, length(disperses));   
    lmlvsperf.(infdecoders{d}).train = NaN(size(disperses));
    lmlvsperf.(infdecoders{d}).test = NaN(size(disperses));
    lmlvsperf.(infdecoders{d}).inference = NaN(size(disperses));
    for s= 1:numel(disperses)
        lmlvscm.(infdecoders{d}).train(:,:,s) = mean(tempdecoderlmlv(s).trainacc,3);
        lmlvscm.(infdecoders{d}).test(:,:,s) = mean(tempdecoderlmlv(s).testacc,3);
        lmlvscm.(infdecoders{d}).inference(:,:,s) = mean(tempdecoderlmlv(s).infperf,3);
        
        lmlvsperf.(infdecoders{d}).train(s) = mean(diag(mean(tempdecoderlmlv(s).trainacc,3)));
        lmlvsperf.(infdecoders{d}).test(s) = mean(diag(mean(tempdecoderlmlv(s).testacc,3)));
        tempinfperf = squeeze(mean(tempdecoderlmlv(s).infperf,3));
        lmlvsperf.(infdecoders{d}).inference(s) = ( tempinfperf(1,1)-tempinfperf(1,2)-tempinfperf(2,3)+tempinfperf(2,4) )/2;
    end

    asiscmagg.(infdecoders{d}).train(:,:,ises) = asiscm.(infdecoders{d}).train;
    asiscmagg.(infdecoders{d}).test(:,:,ises) = asiscm.(infdecoders{d}).test;
    asiscmagg.(infdecoders{d}).inference(:,:,ises) = asiscm.(infdecoders{d}).inference;

    lmlvscmagg.(infdecoders{d}).train(:,:,ises,:) = lmlvscm.(infdecoders{d}).train;
    lmlvscmagg.(infdecoders{d}).test(:,:,ises,:) = lmlvscm.(infdecoders{d}).test;
    lmlvscmagg.(infdecoders{d}).inference(:,:,ises,:) = lmlvscm.(infdecoders{d}).inference;

    asisperfagg.(infdecoders{d}).train(ises) = asisperf.(infdecoders{d}).train;
    asisperfagg.(infdecoders{d}).test(ises) = asisperf.(infdecoders{d}).test;
    asisperfagg.(infdecoders{d}).inference(ises) = asisperf.(infdecoders{d}).inference;

    lmlvsperfagg.(infdecoders{d}).train(ises,:) = lmlvsperf.(infdecoders{d}).train;
    lmlvsperfagg.(infdecoders{d}).test(ises,:) = lmlvsperf.(infdecoders{d}).test;
    lmlvsperfagg.(infdecoders{d}).inference(ises,:) = lmlvsperf.(infdecoders{d}).inference;
end
end

save([drivepath 'RESEARCH/logmean_logvar/OpenScope_bayesdecodespkcnt_ICwcfg1.mat'], ...
    'disperses', 'infdecoders', 'asiscmagg', 'lmlvscmagg', 'asisperfagg', 'lmlvsperfagg')

 
%%
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_bayesdecodespkcnt_ICwcfg1.mat'])

figure
for d = 1:numel(infdecoders)
    subplot(2,4,d)
    hold all
    plot(disperses, lmlvsperfagg.(infdecoders{d}).test, '-')
    plot(disperses, mean(lmlvsperfagg.(infdecoders{d}).test,1), 'k-', 'LineWidth', 2)
    xl = [disperses(1) disperses(end)];
    plot(xl, mean(asisperfagg.(infdecoders{d}).test,1)*[1 1], 'c-', 'LineWidth', 2)
    yl = ylim;
    plot([1 1], yl, 'k--')
    ylim(yl)
    xlim([0 xl(2)])
    xlabel('log(mean) vs log(var) slopes')
    ylabel('test accuracy')
    title(infdecoders{d})

    subplot(2,4,4+d)
    hold all
    plot(disperses, lmlvsperfagg.(infdecoders{d}).inference, '-')
    plot(disperses, mean(lmlvsperfagg.(infdecoders{d}).inference,1), 'k-', 'LineWidth', 2)
    xl = [disperses(1) disperses(end)];
    plot(xl, mean(asisperfagg.(infdecoders{d}).inference,1)*[1 1], 'c-', 'LineWidth', 2)
    yl = ylim;
    plot([1 1], yl, 'k--')
    ylim(yl)
    xlim([0 xl(2)])
    xlabel('log(mean) vs log(var) slopes')
    ylabel('inference performance')
    title(infdecoders{d})
end


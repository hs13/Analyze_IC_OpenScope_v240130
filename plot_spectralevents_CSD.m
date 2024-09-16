load('/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/IllusionOpenScope/sub-620333/LFP_CSD_probeC.mat')
load('/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/IllusionOpenScope/sub-620333/LFP_psth_probeC.mat')
load('/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/IllusionOpenScope/sub-620333/LFP_TFR_L23_probeC.mat')




kerwinhalf = 5; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));
lfpconv = convn(lfpresamp, kergauss, 'same');

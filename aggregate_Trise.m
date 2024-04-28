% find the first spike timing that exceeds 97.5 percentile
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);

Twin = 5;
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
whichblock = 'ICwcfg1_presentations';

psthblankconvpeakagg = cell(numel(probes), Nsessions);
psthblankconvTpeakagg = cell(numel(probes), Nsessions);
blankmaxdrateagg = cell(numel(probes), Nsessions);
blankTmaxriseagg = cell(numel(probes), Nsessions);
% psthconvavgagg
% psthconvprctagg
psthconvglobprctagg = cell(numel(probes), Nsessions);
psthconvavgpeakagg = cell(numel(probes), Nsessions);
psthconvavgTpeakagg = cell(numel(probes), Nsessions);
psthconvavgpeakprctagg = cell(numel(probes), Nsessions);
maxdrateagg = cell(numel(probes), Nsessions);
Tmaxriseagg = cell(numel(probes), Nsessions);
maxdrateprctagg = cell(numel(probes), Nsessions);

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    sesclk = tic;
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];


    for iprobe = 1:numel(probes)
        Trisefn = [pathpp 'Trise' num2str(Twin) 'msbins_probe' probes{iprobe} '.mat'];
        load(Trisefn)

psthblankconvpeakagg{iprobe, ises} = psthblankconvpeak;
psthblankconvTpeakagg{iprobe, ises} = psthblankconvTpeak;
blankmaxdrateagg{iprobe, ises} = blankmaxdrate;
blankTmaxriseagg{iprobe, ises} = blankTmaxrise;
% psthconvavgagg
% psthconvprctagg
% psthconvglobprctagg
psthconvavgpeakagg{iprobe, ises} = psthconvavgpeak;
psthconvavgTpeakagg{iprobe, ises} = psthconvavgTpeak;
psthconvavgpeakprctagg{iprobe, ises} = psthconvavgpeakprct;
maxdrateagg{iprobe, ises} = maxdrate;
Tmaxriseagg{iprobe, ises} = Tmaxrise;
maxdrateprctagg{iprobe, ises} = maxdrateprct;

    end
    toc(sesclk)
end

psthblankconvpeakall = cat(2, psthblankconvpeakagg{:});
psthblankconvTpeakall = cat(2, psthblankconvTpeakagg{:});
blankmaxdrateall = cat(2, blankmaxdrateagg{:});
blankTmaxriseall = cat(2, blankTmaxriseagg{:});
% psthconvavgagg
% psthconvprctagg
% psthconvglobprctagg
psthconvavgpeakall = cat(2, psthconvavgpeakagg{:});
psthconvavgTpeakall = cat(2, psthconvavgTpeakagg{:});
psthconvavgpeakprctall = cat(2, psthconvavgpeakprctagg{:});
maxdrateall = cat(2, maxdrateagg{:});
Tmaxriseall = cat(2, Tmaxriseagg{:});
maxdrateprctall = cat(2, maxdrateprctagg{:});

save([datadir 'Trise' num2str(Twin) 'msbinsall.mat'], 'vistrialrep', 'vistrialtypes', ...
    'psthblankconvpeakall', 'psthblankconvTpeakall', 'blankmaxdrateall', 'blankTmaxriseall', ...
    'psthconvavgpeakall', 'psthconvavgTpeakall', 'psthconvavgpeakprctall', ...
    'maxdrateall', 'Tmaxriseall', 'maxdrateprctall')
copyfile([datadir 'Trise' num2str(Twin) 'msbinsall.mat'], 'G:\My Drive\RESEARCH\ICexpts_revision23\')

%%
drivepath = 'G:\My Drive\RESEARCH\ICexpts_revision23\';
load([drivepath 'openscope_popavg_all.mat'])
load([drivepath 'Trise' num2str(Twin) 'msbinsall.mat'])

%% V1 vs LM or V1 vs AL
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
ICtrialtypedescription = {'Blank', 'X', 'T_C_1', 'I_C_1', 'L_C_1', 'T_C_2', 'L_C_2', 'I_C_2', ...
    'I_R_E_1', 'I_R_E_2', 'T_R_E_1', 'T_R_E_2', 'X_R_E_1', 'X_R_E_2', ...
    'In_B_R', 'In_B_L', 'In_T_L', 'In_T_R', 'Out_B_R', 'Out_B_L', 'Out_T_L', 'Out_T_L'};

neuV1 = contains(neulocall, 'VISp') & ~contains(neulocall, 'VISpm');
neuLM = contains(neulocall, 'VISl');
neuAL = contains(neulocall, 'VISal');
neuRS = unit_wfdur_all>0.4;
neufilt = (unit_isi_violations_all<0.5 & unit_amplitude_cutoff_all<0.5 & unit_presence_ratio_all>0.9);

visctx = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
neuvis = zeros(size(neulocall));
for a = 1:numel(visareas)
    if strcmp(visareas{a}, 'VISp')
neuinarea = contains(neulocall, 'VISp') & ~contains(neulocall, 'VISpm');        
    else
neuinarea = contains(neulocall, visareas{a});        
    end
    neuvis(neuinarea) = a;
end

validTpeak = psthconvavgpeakprctall>0.975;

typi = ICtrialtypes==1105;
figure; hold all
histogram( psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS), 'binwidth', 5 )
histogram( psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuLM & neuRS), 'binwidth', 5 )

validTpeak = psthconvavgpeakprctall>0.975;
ttoi = [106 107 110 111 506 511 1105 1109 1201 1299];
for t = 1:numel(ttoi)
    typi = ICtrialtypes==ttoi(t);
TpeakV1RS = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS);
TpeakHVARS = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuAL & neuRS);
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RS)) ones(1,numel(TpeakHVARS))], [TpeakV1RS TpeakHVARS], '1');
p = ranksum(TpeakV1RS, TpeakHVARS);
fprintf('Tpeak trial type %d V1RS vs ALRS AUROC=%.4f, p=%.4f\n', ttoi(t), AUC, p)
end
% Tpeak trial type 106 V1RS vs LMRS AUROC=0.4428, p=0.0105
% Tpeak trial type 107 V1RS vs LMRS AUROC=0.4215, p=0.0001
% Tpeak trial type 110 V1RS vs LMRS AUROC=0.4621, p=0.0982
% Tpeak trial type 111 V1RS vs LMRS AUROC=0.4305, p=0.0006
% Tpeak trial type 506 V1RS vs LMRS AUROC=0.4925, p=0.5483
% Tpeak trial type 511 V1RS vs LMRS AUROC=0.4838, p=0.1911
% Tpeak trial type 1105 V1RS vs LMRS AUROC=0.5206, p=0.2831
% Tpeak trial type 1109 V1RS vs LMRS AUROC=0.5280, p=0.1493
% Tpeak trial type 1201 V1RS vs LMRS AUROC=0.5146, p=0.4381
% Tpeak trial type 1299 V1RS vs LMRS AUROC=0.5432, p=0.0224

% Tpeak trial type 106 V1RS vs ALRS AUROC=0.3624, p=0.0000
% Tpeak trial type 107 V1RS vs ALRS AUROC=0.4118, p=0.0000
% Tpeak trial type 110 V1RS vs ALRS AUROC=0.3312, p=0.0000
% Tpeak trial type 111 V1RS vs ALRS AUROC=0.4013, p=0.0000
% Tpeak trial type 506 V1RS vs ALRS AUROC=0.4835, p=0.1633
% Tpeak trial type 511 V1RS vs ALRS AUROC=0.4819, p=0.1276
% Tpeak trial type 1105 V1RS vs ALRS AUROC=0.4768, p=0.2159
% Tpeak trial type 1109 V1RS vs ALRS AUROC=0.5172, p=0.3686
% Tpeak trial type 1201 V1RS vs ALRS AUROC=0.4941, p=0.7544
% Tpeak trial type 1299 V1RS vs ALRS AUROC=0.5094, p=0.6177

validTmaxrise = maxdrateprctall>0.975;
ttoi = [106 107 110 111 506 511 1105 1109 1201 1299];
for t = 1:numel(ttoi)
    typi = ICtrialtypes==ttoi(t);
TmaxriseV1RS = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS);
TmaxriseHVARS = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuLM & neuRS);
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RS)) ones(1,numel(TmaxriseHVARS))], [TmaxriseV1RS TmaxriseHVARS], '1');
p = ranksum(TmaxriseV1RS, TmaxriseHVARS);
fprintf('Tmaxrise trial type %d V1RS vs LMRS AUROC=%.4f, p=%.4f\n', ttoi(t), AUC, p)
end
% Tmaxrise trial type 106 V1RS vs LMRS AUROC=0.4763, p=0.4594
% Tmaxrise trial type 107 V1RS vs LMRS AUROC=0.4588, p=0.1209
% Tmaxrise trial type 110 V1RS vs LMRS AUROC=0.4899, p=0.7642
% Tmaxrise trial type 111 V1RS vs LMRS AUROC=0.4189, p=0.0033
% Tmaxrise trial type 506 V1RS vs LMRS AUROC=0.4858, p=0.2455
% Tmaxrise trial type 511 V1RS vs LMRS AUROC=0.4850, p=0.2169
% Tmaxrise trial type 1105 V1RS vs LMRS AUROC=0.5426, p=0.0879
% Tmaxrise trial type 1109 V1RS vs LMRS AUROC=0.5629, p=0.0128
% Tmaxrise trial type 1201 V1RS vs LMRS AUROC=0.5495, p=0.0399
% Tmaxrise trial type 1299 V1RS vs LMRS AUROC=0.5645, p=0.0095

% Tmaxrise trial type 106 V1RS vs ALRS AUROC=0.4288, p=0.0161
% Tmaxrise trial type 107 V1RS vs ALRS AUROC=0.4657, p=0.2177
% Tmaxrise trial type 110 V1RS vs ALRS AUROC=0.3725, p=0.0000
% Tmaxrise trial type 111 V1RS vs ALRS AUROC=0.4267, p=0.0119
% Tmaxrise trial type 506 V1RS vs ALRS AUROC=0.4809, p=0.1025
% Tmaxrise trial type 511 V1RS vs ALRS AUROC=0.4856, p=0.2179
% Tmaxrise trial type 1105 V1RS vs ALRS AUROC=0.5305, p=0.2192
% Tmaxrise trial type 1109 V1RS vs ALRS AUROC=0.5555, p=0.0307
% Tmaxrise trial type 1201 V1RS vs ALRS AUROC=0.5470, p=0.0622
% Tmaxrise trial type 1299 V1RS vs ALRS AUROC=0.5550, p=0.0319

validTmaxrise = maxdrateprctall>0.975;
% ttoi = [106 107 110 111 506 511 1105 1109 1201 1299];
ttoi = [106 111 1105 1109 1201 1299];
for t = 1:numel(ttoi)
    typi = ICtrialtypes==ttoi(t);
    neuoi = validTmaxrise(typi,:)' & neuvis>0 & neuRS;
    [p,tbl,stats]=kruskalwallis(Tmaxriseall(typi,neuoi), neuvis(neuoi));
    figure; C=multcompare(stats);
end


%% V1 IC-encoders vs V1 segment responders
validTpeak = psthconvavgpeakprctall>0.975;
ctrloptions = {'ICresp', 'segresp', 'inducerencoder', 'ICsigBK'};
for o = 1:numel(ctrloptions)
    ctrlopt = ctrloptions{o};
typi = ICtrialtypes==106;
TpeakV1RSICenc1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICencoder1);
switch ctrlopt
    case 'ICresp'
TpeakV1RSICctrl1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICresp1);
    case 'segresp'
TpeakV1RSICctrl1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3));
    case 'inducerencoder'
TpeakV1RSICctrl1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indenc1 | ICsigall.(whichblock).indenc3));
    case 'ICsigBK'
TpeakV1RSICctrl1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).indenc13);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSICenc1)) ones(1,numel(TpeakV1RSICctrl1))], [TpeakV1RSICenc1 TpeakV1RSICctrl1], '1');
p = ranksum(TpeakV1RSICenc1, TpeakV1RSICctrl1);
fprintf('Tpeak IC1 trials V1RS IC1-encoder vs IC1-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)

typi = ICtrialtypes==111;
TpeakV1RSICenc2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICencoder2);
switch ctrlopt
    case 'ICresp'
TpeakV1RSICctrl2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICresp2);
    case 'segresp'
TpeakV1RSICctrl2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4));
    case 'inducerencoder'
TpeakV1RSICctrl2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indenc2 | ICsigall.(whichblock).indenc4));
    case 'ICsigBK'
TpeakV1RSICctrl2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).indenc24);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSICenc2)) ones(1,numel(TpeakV1RSICctrl2))], [TpeakV1RSICenc2 TpeakV1RSICctrl2], '1');
p = ranksum(TpeakV1RSICenc2, TpeakV1RSICctrl2);
fprintf('Tpeak IC2 trials V1RS IC2-encoder vs IC2-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)

TpeakV1RSICenc = [TpeakV1RSICenc1 TpeakV1RSICenc2];
TpeakV1RSICctrl = [TpeakV1RSICctrl1 TpeakV1RSICctrl2];
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSICenc)) ones(1,numel(TpeakV1RSICctrl))], [TpeakV1RSICenc TpeakV1RSICctrl], '1');
p = ranksum(TpeakV1RSICenc, TpeakV1RSICctrl);
fprintf('Tpeak IC trials V1RS IC-encoder vs IC-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)
end
% Tpeak IC1 trials V1RS IC1-encoder vs IC1-ICresp AUROC=0.4315, p=0.4612
% Tpeak IC2 trials V1RS IC2-encoder vs IC2-ICresp AUROC=0.5485, p=0.5116
% Tpeak IC trials V1RS IC-encoder vs IC-ICresp AUROC=0.5060, p=0.9175
% Tpeak IC1 trials V1RS IC1-encoder vs IC1-segresp AUROC=0.4051, p=0.3146
% Tpeak IC2 trials V1RS IC2-encoder vs IC2-segresp AUROC=0.5272, p=0.7170
% Tpeak IC trials V1RS IC-encoder vs IC-segresp AUROC=0.4806, p=0.7409
% Tpeak IC1 trials V1RS IC1-encoder vs IC1-inducerencoder AUROC=0.4542, p=0.6459
% Tpeak IC2 trials V1RS IC2-encoder vs IC2-inducerencoder AUROC=0.5472, p=0.5536
% Tpeak IC trials V1RS IC-encoder vs IC-inducerencoder AUROC=0.5160, p=0.7957
% Tpeak IC1 trials V1RS IC1-encoder vs IC1-ICsigBK AUROC=0.4398, p=0.5183
% Tpeak IC2 trials V1RS IC2-encoder vs IC2-ICsigBK AUROC=0.5528, p=0.4752
% Tpeak IC trials V1RS IC-encoder vs IC-ICsigBK AUROC=0.5124, p=0.8302

validTmaxrise = maxdrateprctall>0.975;
ctrloptions = {'ICresp', 'segresp', 'inducerencoder', 'ICsigBK'};
for o = 1:numel(ctrloptions)
    ctrlopt = ctrloptions{o};
typi = ICtrialtypes==106;
TmaxriseV1RSICenc1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICencoder1);
switch ctrlopt
    case 'ICresp'
TmaxriseV1RSICctrl1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICresp1);
    case 'segresp'
TmaxriseV1RSICctrl1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3));
    case 'inducerencoder'
TmaxriseV1RSICctrl1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indenc1 | ICsigall.(whichblock).indenc3));
    case 'ICsigBK'
TmaxriseV1RSICctrl1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).indenc13);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSICenc1)) ones(1,numel(TmaxriseV1RSICctrl1))], [TmaxriseV1RSICenc1 TmaxriseV1RSICctrl1], '1');
p = ranksum(TmaxriseV1RSICenc1, TmaxriseV1RSICctrl1);
fprintf('Tmaxrise IC1 trials V1RS IC1-encoder vs IC1-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)

typi = ICtrialtypes==111;
TmaxriseV1RSICenc2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICencoder2);
switch ctrlopt
    case 'ICresp'
TmaxriseV1RSICctrl2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).ICresp2);
    case 'segresp'
TmaxriseV1RSICctrl2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4));
    case 'inducerencoder'
TmaxriseV1RSICctrl2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & (ICsigall.(whichblock).indenc2 | ICsigall.(whichblock).indenc4));
    case 'ICsigBK'
TmaxriseV1RSICctrl2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & ICsigall.(whichblock).indenc24);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSICenc2)) ones(1,numel(TmaxriseV1RSICctrl2))], [TmaxriseV1RSICenc2 TmaxriseV1RSICctrl2], '1');
p = ranksum(TmaxriseV1RSICenc2, TmaxriseV1RSICctrl2);
fprintf('Tmaxrise IC2 trials V1RS IC2-encoder vs IC2-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)

TmaxriseV1RSICenc = [TmaxriseV1RSICenc1 TmaxriseV1RSICenc2];
TmaxriseV1RSICctrl = [TmaxriseV1RSICctrl1 TmaxriseV1RSICctrl2];
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSICenc)) ones(1,numel(TmaxriseV1RSICctrl))], [TmaxriseV1RSICenc TmaxriseV1RSICctrl], '1');
p = ranksum(TmaxriseV1RSICenc, TmaxriseV1RSICctrl);
fprintf('Tmaxrise IC trials V1RS IC-encoder vs IC-%s AUROC=%.4f, p=%.4f\n', ctrlopt, AUC, p)
end
% Tmaxrise IC1 trials V1RS IC1-encoder vs IC1-ICresp AUROC=0.2942, p=0.1625
% Tmaxrise IC2 trials V1RS IC2-encoder vs IC2-ICresp AUROC=0.5216, p=0.8724
% Tmaxrise IC trials V1RS IC-encoder vs IC-ICresp AUROC=0.4248, p=0.4423
% Tmaxrise IC1 trials V1RS IC1-encoder vs IC1-segresp AUROC=0.2772, p=0.1350
% Tmaxrise IC2 trials V1RS IC2-encoder vs IC2-segresp AUROC=0.5147, p=0.9158
% Tmaxrise IC trials V1RS IC-encoder vs IC-segresp AUROC=0.4105, p=0.3643
% Tmaxrise IC1 trials V1RS IC1-encoder vs IC1-inducerencoder AUROC=0.3375, p=0.3097
% Tmaxrise IC2 trials V1RS IC2-encoder vs IC2-inducerencoder AUROC=0.5959, p=0.4923
% Tmaxrise IC trials V1RS IC-encoder vs IC-inducerencoder AUROC=0.4838, p=0.8796
% Tmaxrise IC1 trials V1RS IC1-encoder vs IC1-ICsigBK AUROC=0.2932, p=0.1609
% Tmaxrise IC2 trials V1RS IC2-encoder vs IC2-ICsigBK AUROC=0.5204, p=0.8795
% Tmaxrise IC trials V1RS IC-encoder vs IC-ICsigBK AUROC=0.4240, p=0.4374

%% V1 center-RF vs V1 segment-RF vs LM/AL center-RF
validTmaxrise = maxdrateprctall>0.975;
rfsigoptions = {'Pkw_rfclassic', 'pRFclassic'};
for o = 1:numel(rfsigoptions)
    neuopt = rfsigoptions{o};
    disp(neuopt)
typi = ICtrialtypes==106;
switch neuopt
    case 'Pkw_rfclassic'
TmaxriseV1RSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseV1RSsegRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]) );
TmaxriseLMRSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuLM & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseALRSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuAL & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
    case 'pRFclassic'
TmaxriseV1RSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseV1RSsegRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]) );
TmaxriseLMRSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuLM & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseALRSctrRFIC1 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuAL & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC1)) ones(1,numel(TmaxriseV1RSsegRFIC1))], ...
    [TmaxriseV1RSctrRFIC1 TmaxriseV1RSsegRFIC1], '1');
p = ranksum(TmaxriseV1RSctrRFIC1, TmaxriseV1RSsegRFIC1);
fprintf('Tmaxrise IC1 trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC1)) ones(1,numel(TmaxriseLMRSctrRFIC1))], ...
    [TmaxriseV1RSctrRFIC1 TmaxriseLMRSctrRFIC1], '1');
p = ranksum(TmaxriseV1RSctrRFIC1, TmaxriseLMRSctrRFIC1);
fprintf('Tmaxrise IC1 trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC1)) ones(1,numel(TmaxriseALRSctrRFIC1))], ...
    [TmaxriseV1RSctrRFIC1 TmaxriseALRSctrRFIC1], '1');
p = ranksum(TmaxriseV1RSctrRFIC1, TmaxriseALRSctrRFIC1);
fprintf('Tmaxrise IC1 trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

typi = ICtrialtypes==111;
switch neuopt
    case 'Pkw_rfclassic'
TmaxriseV1RSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseV1RSsegRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]) );
TmaxriseLMRSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuLM & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseALRSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuAL & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
    case 'pRFclassic'
TmaxriseV1RSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseV1RSsegRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]) );
TmaxriseLMRSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuLM & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TmaxriseALRSctrRFIC2 = Tmaxriseall(typi,validTmaxrise(typi,:)' & neuAL & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
end

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC2)) ones(1,numel(TmaxriseV1RSsegRFIC2))], ...
    [TmaxriseV1RSctrRFIC2 TmaxriseV1RSsegRFIC2], '1');
p = ranksum(TmaxriseV1RSctrRFIC2, TmaxriseV1RSsegRFIC2);
fprintf('Tmaxrise IC2 trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC2)) ones(1,numel(TmaxriseLMRSctrRFIC2))], ...
    [TmaxriseV1RSctrRFIC2 TmaxriseLMRSctrRFIC2], '1');
p = ranksum(TmaxriseV1RSctrRFIC2, TmaxriseLMRSctrRFIC2);
fprintf('Tmaxrise IC2 trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRFIC2)) ones(1,numel(TmaxriseALRSctrRFIC2))], ...
    [TmaxriseV1RSctrRFIC2 TmaxriseALRSctrRFIC2], '1');
p = ranksum(TmaxriseV1RSctrRFIC2, TmaxriseALRSctrRFIC2);
fprintf('Tmaxrise IC2 trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

TmaxriseV1RSctrRF = [TmaxriseV1RSctrRFIC1 TmaxriseV1RSctrRFIC2];
TmaxriseV1RSsegRF = [TmaxriseV1RSsegRFIC1 TmaxriseV1RSsegRFIC2];
TmaxriseLMRSctrRF = [TmaxriseLMRSctrRFIC1 TmaxriseLMRSctrRFIC2];
TmaxriseALRSctrRF = [TmaxriseALRSctrRFIC1 TmaxriseALRSctrRFIC2];
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRF)) ones(1,numel(TmaxriseV1RSsegRF))], ...
    [TmaxriseV1RSctrRF TmaxriseV1RSsegRF], '1');
p = ranksum(TmaxriseV1RSctrRF, TmaxriseV1RSsegRF);
fprintf('Tmaxrise IC trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRF)) ones(1,numel(TmaxriseLMRSctrRF))], ...
    [TmaxriseV1RSctrRF TmaxriseLMRSctrRF], '1');
p = ranksum(TmaxriseV1RSctrRF, TmaxriseLMRSctrRF);
fprintf('Tmaxrise IC trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TmaxriseV1RSctrRF)) ones(1,numel(TmaxriseALRSctrRF))], ...
    [TmaxriseV1RSctrRF TmaxriseALRSctrRF], '1');
p = ranksum(TmaxriseV1RSctrRF, TmaxriseALRSctrRF);
fprintf('Tmaxrise IC trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)
end


%% V1 center-RF vs V1 segment-RF vs LM/AL center-RF
validTpeak = psthconvavgpeakprctall>0.975;
rfsigoptions = {'Pkw_rfclassic', 'pRFclassic'};
for o = 1:numel(rfsigoptions)
    neuopt = rfsigoptions{o};
    disp(neuopt)
typi = ICtrialtypes==106;
switch neuopt
    case 'Pkw_rfclassic'
TpeakV1RSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TpeakV1RSsegRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]) );
TpeakLMRSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuLM & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TpeakALRSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuAL & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
    case 'pRFclassic'
TpeakV1RSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TpeakV1RSsegRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]) );
TpeakLMRSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuLM & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TpeakALRSctrRFIC1 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuAL & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
end
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC1)) ones(1,numel(TpeakV1RSsegRFIC1))], ...
    [TpeakV1RSctrRFIC1 TpeakV1RSsegRFIC1], '1');
p = ranksum(TpeakV1RSctrRFIC1, TpeakV1RSsegRFIC1);
fprintf('Tpeak IC1 trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC1)) ones(1,numel(TpeakLMRSctrRFIC1))], ...
    [TpeakV1RSctrRFIC1 TpeakLMRSctrRFIC1], '1');
p = ranksum(TpeakV1RSctrRFIC1, TpeakLMRSctrRFIC1);
fprintf('Tpeak IC1 trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC1)) ones(1,numel(TpeakALRSctrRFIC1))], ...
    [TpeakV1RSctrRFIC1 TpeakALRSctrRFIC1], '1');
p = ranksum(TpeakV1RSctrRFIC1, TpeakALRSctrRFIC1);
fprintf('Tpeak IC1 trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

typi = ICtrialtypes==111;
switch neuopt
    case 'Pkw_rfclassic'
TpeakV1RSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TpeakV1RSsegRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.Pkw_rfclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]) );
TpeakLMRSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuLM & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
TpeakALRSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuAL & neuRS & RFCIall.Pkw_rfclassic<0.05 & RFCIall.RFindclassic==1);
    case 'pRFclassic'
TpeakV1RSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TpeakV1RSsegRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuV1 & neuRS & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]) );
TpeakLMRSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuLM & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
TpeakALRSctrRFIC2 = psthconvavgTpeakall(typi,validTpeak(typi,:)' & neuAL & neuRS & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1);
end

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC2)) ones(1,numel(TpeakV1RSsegRFIC2))], ...
    [TpeakV1RSctrRFIC2 TpeakV1RSsegRFIC2], '1');
p = ranksum(TpeakV1RSctrRFIC2, TpeakV1RSsegRFIC2);
fprintf('Tpeak IC2 trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC2)) ones(1,numel(TpeakLMRSctrRFIC2))], ...
    [TpeakV1RSctrRFIC2 TpeakLMRSctrRFIC2], '1');
p = ranksum(TpeakV1RSctrRFIC2, TpeakLMRSctrRFIC2);
fprintf('Tpeak IC2 trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRFIC2)) ones(1,numel(TpeakALRSctrRFIC2))], ...
    [TpeakV1RSctrRFIC2 TpeakALRSctrRFIC2], '1');
p = ranksum(TpeakV1RSctrRFIC2, TpeakALRSctrRFIC2);
fprintf('Tpeak IC2 trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

TpeakV1RSctrRF = [TpeakV1RSctrRFIC1 TpeakV1RSctrRFIC2];
TpeakV1RSsegRF = [TpeakV1RSsegRFIC1 TpeakV1RSsegRFIC2];
TpeakLMRSctrRF = [TpeakLMRSctrRFIC1 TpeakLMRSctrRFIC2];
TpeakALRSctrRF = [TpeakALRSctrRFIC1 TpeakALRSctrRFIC2];
[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRF)) ones(1,numel(TpeakV1RSsegRF))], ...
    [TpeakV1RSctrRF TpeakV1RSsegRF], '1');
p = ranksum(TpeakV1RSctrRF, TpeakV1RSsegRF);
fprintf('Tpeak IC trials V1RS center-RF vs RF-on-segments AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRF)) ones(1,numel(TpeakLMRSctrRF))], ...
    [TpeakV1RSctrRF TpeakLMRSctrRF], '1');
p = ranksum(TpeakV1RSctrRF, TpeakLMRSctrRF);
fprintf('Tpeak IC trials V1RS vs LMRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)

[X,Y,T,AUC] = perfcurve([zeros(1,numel(TpeakV1RSctrRF)) ones(1,numel(TpeakALRSctrRF))], ...
    [TpeakV1RSctrRF TpeakALRSctrRF], '1');
p = ranksum(TpeakV1RSctrRF, TpeakALRSctrRF);
fprintf('Tpeak IC trials V1RS vs ALRS center-RF AUROC=%.4f, p=%.4f\n', AUC, p)
end


%% for each trial type, compare neuron groups
% 'visareas', 'V1LMAL', 'visareas_ctrRF', 'V1LMAL_ctrRF', 'V1ctrRF_V1onsegRF_LMctrRF_ALctrRF', 'V1ICenc_V1segresp_LMICenc_ALICenc'
whichT = 'Tmaxrise';
whichneugroup = 'V1ICenc_V1segresp_LMICenc'; % V1LM, V1ctrRF_V1onsegRF_LMctrRF, 
% V1ICenc_V1segresp, V1ICenc_V1segresp_LMICenc, V1ICenc_V1segresp_LMICenc_ALICenc all not significant

switch whichT
    case 'Tpeak'
validT = psthconvavgpeakprctall>0.975;
Tmat = psthconvavgTpeakall;
ylab = 'Peak Time (ms)';
    case 'Tmaxrise'
validT = maxdrateprctall>0.975;
Tmat = Tmaxriseall;
ylab = 'Rise Time (ms)';
end

yl = [0 250];
% typi = ICtrialtypes==1109;
ttoi = [106 111 1105 1109 1201 1299];
for t = 1:numel(ttoi)
    typi = ICtrialtypes==ttoi(t);

switch whichneugroup
    case 'visareas'
        neurogrouplabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        neuingroup = neuRS & neuvis>0;
        neurongroups = neuvis;
    case 'V1LMAL'
        neurogrouplabels = {'V1', 'LM', 'AL'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl', 'VISal'})));
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'V1LM'
        neurogrouplabels = {'V1', 'LM'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl'})));
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'visareas_ctrRF'
        neurogrouplabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        neuingroup = neuRS & neuvis>0 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = neuvis;
    case 'V1LMAL_ctrRF'
        neurogrouplabels = {'V1', 'LM', 'AL'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl', 'VISal'}))) & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'V1LM_ctrRF'
        neurogrouplabels = {'V1', 'LM'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl'}))) & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'V1ctrRF_V1onsegRF_LMctrRF_ALctrRF'
        neurogrouplabels = {'V1 center-RF', 'V1 RF-on-seg.', 'LM center-RF', 'AL center-RF'};
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 center-RF'
                    tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'V1 RF-on-seg.'
                    if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]);
                    elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]);
                    else
                        error('designate appropraite RFindclassic for image in IC image set')
                    end
                case 'LM center-RF'
                    tempneu = neuRS & neuLM & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'AL center-RF'
                    tempneu = neuRS & neuAL & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ctrRF_V1onsegRF_LMctrRF'
        neurogrouplabels = {'V1 center-RF', 'V1 RF-on-seg.', 'LM center-RF'};
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 center-RF'
                    tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'V1 RF-on-seg.'
                    if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]);
                    elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]);
                    else
                        error('designate appropraite RFindclassic for image in IC image set')
                    end
                case 'LM center-RF'
                    tempneu = neuRS & neuLM & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'AL center-RF'
                    tempneu = neuRS & neuAL & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ctrRF_V1onsegRF'
        neurogrouplabels = {'V1 center-RF', 'V1 RF-on-seg.'};
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 center-RF'
                    tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'V1 RF-on-seg.'
                    if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]);
                    elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]);
                    else
                        error('designate appropraite RFindclassic for image in IC image set')
                    end
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ICenc_V1segresp_LMICenc_ALICenc'
        neurogrouplabels = {'V1 IC-enc.', 'V1 seg.-resp.', 'LM IC-enc.', 'AL IC-enc.'};

        if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
            neuICenc = ICsigall.(whichblock).ICencoder1;
            neusegresp = ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3;
        elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
            neuICenc = ICsigall.(whichblock).ICencoder2;
            neusegresp = ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4;
        else
            error('designate appropraite RFindclassic for image in IC image set')
        end
        
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 IC-enc.'
                    tempneu = neuRS & neuV1 & neuICenc;
                case 'V1 seg.-resp.'
                    tempneu = neuRS & neuV1 & neusegresp;
                case 'LM IC-enc.'
                    tempneu = neuRS & neuLM & neuICenc;
                case 'AL IC-enc.'
                    tempneu = neuRS & neuAL & neuICenc;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ICenc_V1segresp_LMICenc'
        neurogrouplabels = {'V1 IC-enc.', 'V1 seg.-resp.', 'LM IC-enc.'};

        if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
            neuICenc = ICsigall.(whichblock).ICencoder1;
            neusegresp = ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3;
        elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
            neuICenc = ICsigall.(whichblock).ICencoder2;
            neusegresp = ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4;
        else
            error('designate appropraite RFindclassic for image in IC image set')
        end
        
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 IC-enc.'
                    tempneu = neuRS & neuV1 & neuICenc;
                case 'V1 seg.-resp.'
                    tempneu = neuRS & neuV1 & neusegresp;
                case 'LM IC-enc.'
                    tempneu = neuRS & neuLM & neuICenc;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ICenc_V1segresp'
        neurogrouplabels = {'V1 IC-enc.', 'V1 seg.-resp.'};

        if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
            neuICenc = ICsigall.(whichblock).ICencoder1;
            neusegresp = ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3;
        elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
            neuICenc = ICsigall.(whichblock).ICencoder2;
            neusegresp = ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4;
        else
            error('designate appropraite RFindclassic for image in IC image set')
        end
        
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 IC-enc.'
                    tempneu = neuRS & neuV1 & neuICenc;
                case 'V1 seg.-resp.'
                    tempneu = neuRS & neuV1 & neusegresp;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;
end
neuoi = validT(typi,:)' & neuingroup;

[p,tbl,stats]=kruskalwallis(Tmat(typi,neuoi), neurongroups(neuoi), 'off');
C = multcompare(stats,'display','off');
indsigcomp = find(C(:,end)<0.05);

figure('Position', [100+240*(t-1) 100 240 320])
hold all
b=boxchart(neurongroups(neuoi), Tmat(typi,neuoi), 'notch' , 'on', 'linewidth', 2);
% b.JitterOutliers = 'on';
% b.MarkerStyle = '.';
b.MarkerStyle = 'none';
b.BoxFaceColor = 0.5*[1 1 1];
b.BoxEdgeColor = 0.5*[1 1 1];
for c = 1:numel(indsigcomp)
    h = yl(2)+0.02*range(yl)*(numel(indsigcomp)-c);
    plot([C(indsigcomp(c),1) C(indsigcomp(c),2)], [h h], 'k-')
    plot( mean([C(indsigcomp(c),1) C(indsigcomp(c),2)]), h+0.01*range(yl), 'k*')
end
set(gca, 'XTick', 1:numel(neurogrouplabels), 'XTickLabel', neurogrouplabels, 'FontSize', fs)
ylim([yl(1) yl(2)+0.02*range(yl)*numel(indsigcomp)])
xlim([0.5 numel(neurogrouplabels)+0.5])
ylabel(ylab, 'FontSize', fs)
% title([ICtrialtypedescription{typi} ' Trials'], 'FontSize', fs)
title(sprintf('%s Trials p=%.4f', ICtrialtypedescription{typi}, p), 'FontSize', fs)
end

%% pool IC/TRE/XRE trial types
% 'visareas', 'V1LMAL', 'visareas_ctrRF', 'V1LMAL_ctrRF', 'V1ctrRF_V1onsegRF_LMctrRF_ALctrRF', 'V1ICenc_V1segresp_LMICenc_ALICenc'
whichT = 'Tmaxrise';
whichneugroup = 'V1LM'; % V1LM, V1ctrRF_V1onsegRF_LMctrRF, 
% V1ICenc_V1segresp, V1ICenc_V1segresp_LMICenc, V1ICenc_V1segresp_LMICenc_ALICenc all not significant
poolopt = 'poolIC';
fw=240; fh=320; xtr = 30;

switch whichT
    case 'Tpeak'
validT = psthconvavgpeakprctall>0.975;
Tmat = psthconvavgTpeakall;
ylab = 'Peak Time (ms)';
    case 'Tmaxrise'
validT = maxdrateprctall>0.975;
Tmat = Tmaxriseall;
ylab = 'Rise Time (ms)';
end


yl = [0 200];
switch poolopt
    case 'poolIC'
ttoi = [106 111];
    case 'poolTRE'
ttoi = [1105 1109];
    case 'poolXRE'
ttoi = [1201 1299];
end
Tpool = [];
labelpool = [];
for t = 1:numel(ttoi)
    typi = ICtrialtypes==ttoi(t);

switch whichneugroup
    case 'visareas'
        neurogrouplabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        neuingroup = neuRS & neuvis>0;
        neurongroups = neuvis;
    case 'V1LMAL'
        neurogrouplabels = {'V1', 'LM', 'AL'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl', 'VISal'})));
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'V1LM'
        neurogrouplabels = {'V1', 'LM'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl'})));
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
        % fw=180; fh=240;
        xtr = 0;
    case 'visareas_ctrRF'
        neurogrouplabels = {'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        neuingroup = neuRS & neuvis>0 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = neuvis;
    case 'V1LMAL_ctrRF'
        neurogrouplabels = {'V1', 'LM', 'AL'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl', 'VISal'}))) & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end
    case 'V1LM_ctrRF'
        neurogrouplabels = {'V1', 'LM'};
        neuingroup = neuRS & ismember(neuvis, find(ismember(visareas, {'VISp', 'VISl'}))) & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            neurongroups(neuvis==find(strcmp(visctx,neurogrouplabels{g})))=g;
        end

    case 'V1ctrRF_V1onsegRF_LMctrRF'
        neurogrouplabels = {'V1 center-RF', 'V1 RF-on-seg.', 'LM center-RF'};
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 center-RF'
                    tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'V1 RF-on-seg.'
                    if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[3,7]);
                    elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
                        tempneu = neuRS & neuV1 & RFCIall.pRFclassic<0.05 & ismember(RFCIall.RFindclassic,[5,9]);
                    else
                        error('designate appropraite RFindclassic for image in IC image set')
                    end
                case 'LM center-RF'
                    tempneu = neuRS & neuLM & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
                case 'AL center-RF'
                    tempneu = neuRS & neuAL & RFCIall.pRFclassic<0.05 & RFCIall.RFindclassic==1;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;

    case 'V1ICenc_V1segresp_LMICenc'
        neurogrouplabels = {'V1 IC-enc.', 'V1 seg.-resp.', 'LM IC-enc.'};

        if ismember(ICtrialtypes(typi), [105 106 506 1105 1201])
            neuICenc = ICsigall.(whichblock).ICencoder1;
            neusegresp = ICsigall.(whichblock).indin1 | ICsigall.(whichblock).indin3;
        elseif ismember(ICtrialtypes(typi), [109 111 511 1109 1299])
            neuICenc = ICsigall.(whichblock).ICencoder2;
            neusegresp = ICsigall.(whichblock).indin2 | ICsigall.(whichblock).indin4;
        else
            error('designate appropraite RFindclassic for image in IC image set')
        end
        
        neurongroups = zeros(size(neulocall));
        for g = 1:numel(neurogrouplabels)
            switch neurogrouplabels{g}
                case 'V1 IC-enc.'
                    tempneu = neuRS & neuV1 & neuICenc;
                case 'V1 seg.-resp.'
                    tempneu = neuRS & neuV1 & neusegresp;
                case 'LM IC-enc.'
                    tempneu = neuRS & neuLM & neuICenc;
            end
            neurongroups(tempneu)=g;
        end
        neuingroup = neurongroups>0;
end
neuoi = validT(typi,:)' & neuingroup;
Tpool = [Tpool Tmat(typi,neuoi)];
labelpool = [labelpool neurongroups(neuoi)'];
end

[p,tbl,stats]=kruskalwallis(Tpool, labelpool, 'off');
C = multcompare(stats,'display','off');
indsigcomp = find(C(:,end)<0.05);

figure('Position', [100+300*(t-1) 100 fw fh])
hold all
b=boxchart(labelpool, Tpool, 'notch' , 'on', 'linewidth', 2);
% b.JitterOutliers = 'on';
% b.MarkerStyle = '.';
b.MarkerStyle = 'none';
b.BoxFaceColor = 0.5*[1 1 1];
b.BoxEdgeColor = 0.5*[1 1 1];
for c = 1:numel(indsigcomp)
    h = yl(2)+0.02*range(yl)*(numel(indsigcomp)-c);
    plot([C(indsigcomp(c),1) C(indsigcomp(c),2)], [h h], 'k-')
    plot( mean([C(indsigcomp(c),1) C(indsigcomp(c),2)]), h+0.01*range(yl), 'k*')
end
text(numel(neurogrouplabels)+0.5, yl(1), sprintf('p=%.4f',p), 'FontSize', fs, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
set(gca, 'XTick', 1:numel(neurogrouplabels), 'XTickLabel', neurogrouplabels, 'XTickLabelRotation', xtr, 'FontSize', fs)
ylim([yl(1) yl(2)+0.02*range(yl)*numel(indsigcomp)])
xlim([0.5 numel(neurogrouplabels)+0.5])
ylabel(ylab, 'FontSize', fs)
switch poolopt
    case 'poolIC'
title(sprintf('I_C Trials p=%.4f', p), 'FontSize', fs, 'FontWeight', 'normal')
    case 'poolTRE'
title(sprintf('T_R_E Trials p=%.4f', p), 'FontSize', fs, 'FontWeight', 'normal')
    case 'poolXRE'
title(sprintf('X_R_E Trials p=%.4f', p), 'FontSize', fs, 'FontWeight', 'normal')
end
switch poolopt
    case 'poolIC'
title('I_C Trials')
    case 'poolTRE'
title('T_R_E Trials')
    case 'poolXRE'
title('X_R_E Trials')
end

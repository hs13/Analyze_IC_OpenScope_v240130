addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v230821')

datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

gazedistthresh = 20;

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

%% print genotype
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
genotypes = cell(size(nwbsessions));
tic
for ises = 1:numel(nwbsessions)
nwbfiles = cat(1, dir([datadir nwbsessions{ises} '/*.nwb']), dir([datadir nwbsessions{ises} '/*/*.nwb']));

% take filename  with shortest length or filename that does not contain probe
[~, fileind] = min(cellfun(@length, {nwbfiles.name}));
nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);

nwb = nwbRead(nwbspikefile);

genotypes{ises} = nwb.general_subject.genotype;
disp(nwb.general_subject.genotype)
end
toc

% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Pvalb-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt
% Sst-IRES-Cre/wt;Ai32(RCL-ChR2(H134R)_EYFP)/wt

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=14;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

unit_amplitude_agg = cell(1,Nsessions);
unit_isi_violations_agg = cell(1,Nsessions);
unit_wfdur_agg = cell(1,Nsessions);
unit_amplitude_cutoff_agg = cell(1,Nsessions);
unit_presence_ratio_agg = cell(1,Nsessions);

neuoindagg = cell(numel(probes),Nsessions);
probeneuagg = cell(numel(probes),Nsessions);
neulocagg = cell(numel(probes),Nsessions);
neupeakchagg = cell(numel(probes),Nsessions);
neuctxagg = cell(numel(probes),Nsessions);
sesneuagg = cell(numel(probes),Nsessions);

meanFRvecagg = cell(size(probes));
sponFRvecagg = cell(size(probes));

for b = 1:numel(ICblocks)
    ICsigagg.(ICblocks{b}) = cell(numel(probes), Nsessions);
end
RFCIagg = cell(numel(probes), Nsessions);
RFCIspinagg = cell(numel(probes), Nsessions);
sizeCIagg = cell(numel(probes), Nsessions);
oriparamsagg = cell(numel(probes), Nsessions);
ori4paramsagg = cell(numel(probes), Nsessions);


% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
allfields = fieldnames(ICsig_fixedgaze.ICwcfg1_presentations);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),2);
end
ICsigfields = allfields(validfields);
ICsigfieldsize2 = size2fields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
ICsig_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsig_fixedgazeagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) = [];
        end
    end
end

allfields = fieldnames(RFCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCI_fixedgaze.(allfields{f}),2);
end
RFCIfields = allfields(validfields);
RFCIfieldsize2 = size2fields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCI_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIfields)
        RFCI_fixedgazeagg.(probes{iprobe}).(RFCIfields{f}) = [];
    end
end

allfields = fieldnames(RFCIspin_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCIspin_fixedgaze.(allfields{f}),2);
end
RFCIspinfields = allfields(validfields);
RFCIspinfieldsize2 = size2fields(validfields);
RFCIspin_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIspinfields)
        RFCIspin_fixedgazeagg.(probes{iprobe}).(RFCIspinfields{f}) = [];
    end
end

allfields = fieldnames(sizeCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(sizeCI_fixedgaze.(allfields{f}),2);
end
sizeCIfields = allfields(validfields);
sizeCIfieldsize2 = size2fields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCI_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(sizeCIfields)
        sizeCI_fixedgazeagg.(probes{iprobe}).(sizeCIfields{f}) = [];
    end
end

allfields = fieldnames(oriparams_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(oriparams_fixedgaze.(allfields{f}),2);
end
oriparamsfields = allfields(validfields);
oriparamsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
oriparams_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(oriparamsfields)
        oriparams_fixedgazeagg.(probes{iprobe}).(oriparamsfields{f}) = [];
    end
end

allfields = fieldnames(ori4params_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ori4params_fixedgaze.(allfields{f}),2);
end
ori4paramsfields = allfields(validfields);
ori4paramsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
ori4params_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(ori4paramsfields)
        ori4params_fixedgazeagg.(probes{iprobe}).(ori4paramsfields{f}) = [];
    end
end

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    if ~exist([pathpp 'trackmouseeye.mat'], 'file')
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        continue
    end
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        tic
%         load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'


        if exist([pathpp, 'probes.mat'], 'file')
            probelist = load([pathpp, 'probes.mat']);
            warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        end
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        
        
%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        if nnz(floor(unit_peakch/1000)==probeind-1)==0
            fprintf('Probe %s Area %s: NO UNITS!!!\n', probes{iprobe}, visareas{iprobe} )
            continue
        end
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
%         disp(unique(probelocs)')
        
        load(sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        probeneuagg{iprobe} = cat(1, probeneuagg{iprobe}, neuoind);
        neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
        neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
        neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
        
        sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
        sesneuctx_fixedgazeagg{iprobe} = cat(1, sesneuctx_fixedgazeagg{iprobe}, ises*ones(nnz(neuctx),1));
        
%         % aggregate cortex neurons
%         neuoi = neuctx;
        % aggregate all neurons
        neuoi = true(size(neuoind));
        
        meanFRvecagg{iprobe} = cat(1, meanFRvecagg{iprobe}, meanFRvec(neuoi)');
        sponFRvecagg{iprobe} = cat(1, sponFRvecagg{iprobe}, sponFRvec(neuoi)');
        
        for b = 1:numel(ICblocks)
            for f= 1:numel(ICsigfields)
                if ~isfield(ICsig_fixedgaze.(ICblocks{b}), ICsigfields{f})
                    tempmat = NaN( nnz(neuoi), ICsigfieldsize2(f) );
                else
                    tempmat = ICsig_fixedgaze.(ICblocks{b}).(ICsigfields{f})(neuoi,:);
                end
                ICsig_fixedgazeagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) = ...
                    cat(1, ICsig_fixedgazeagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}), tempmat );
            end
        end
        
        % RFCIfields = fieldnames(RFCI);
        for f= 1:numel(RFCIfields)
            if ~isfield(RFCI_fixedgaze, RFCIfields{f})
                tempmat = NaN( nnz(neuoi), RFCIfieldsize2(f) );
            else
                tempmat = RFCI_fixedgaze.(RFCIfields{f})(neuoi,:);
            end
            RFCI_fixedgazeagg.(probes{iprobe}).(RFCIfields{f}) = cat(1, ...
                RFCI_fixedgazeagg.(probes{iprobe}).(RFCIfields{f}), tempmat );
        end
        
        for f= 1:numel(RFCIspinfields)
            if ~isfield(RFCIspin_fixedgaze, RFCIspinfields{f})
                tempmat = NaN( nnz(neuoi), RFCIspinfieldsize2(f) );
            else
                tempmat = RFCIspin_fixedgaze.(RFCIspinfields{f})(neuoi,:);
            end
            RFCIspin_fixedgazeagg.(probes{iprobe}).(RFCIspinfields{f}) = cat(1, ...
                RFCIspin_fixedgazeagg.(probes{iprobe}).(RFCIspinfields{f}), tempmat );
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        for f= 1:numel(sizeCIfields)
            if ~isfield(sizeCI_fixedgaze, sizeCIfields{f})
                tempmat = NaN( nnz(neuoi), sizeCIfieldsize2(f) );
            else
                tempmat = sizeCI_fixedgaze.(sizeCIfields{f})(neuoi,:);
            end
            sizeCI_fixedgazeagg.(probes{iprobe}).(sizeCIfields{f}) = cat(1, ...
                sizeCI_fixedgazeagg.(probes{iprobe}).(sizeCIfields{f}), tempmat );
        end
        
        for f= 1:numel(oriparamsfields)
            if ~isfield(oriparams_fixedgaze, oriparamsfields{f})
                tempmat = NaN( nnz(neuoi), oriparamsfieldsize2(f) );
            else
                tempmat = oriparams_fixedgaze.(oriparamsfields{f})(neuoi,:);
            end
            oriparams_fixedgazeagg.(probes{iprobe}).(oriparamsfields{f}) = cat(1, ...
                oriparams_fixedgazeagg.(probes{iprobe}).(oriparamsfields{f}), tempmat );
        end
        
        for f= 1:numel(ori4paramsfields)
            if ~isfield(ori4params_fixedgaze, ori4paramsfields{f})
                tempmat = NaN( nnz(neuoi), ori4paramsfieldsize2(f) );
            else
                tempmat = ori4params_fixedgaze.(ori4paramsfields{f})(neuoi,:);
            end
            ori4params_fixedgazeagg.(probes{iprobe}).(ori4paramsfields{f}) = cat(1, ...
                ori4params_fixedgazeagg.(probes{iprobe}).(ori4paramsfields{f}), tempmat );
        end
        
        toc
    end
    
end

Nrfs = RFCIfieldsize2(strcmp(RFCIfields, 'Rrfclassic'));
Nszs = sizeCIfieldsize2(strcmp(sizeCIfields, 'Rsizeclassic'));
Ndirs = oriparamsfieldsize2(strcmp(oriparamsfields, 'Rori'));
Noris = Ndirs/2;

ises=4;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

save(['/Users/hyeyoung/Documents/DATA/openscope_popavg_fixedgazeagg.mat'], ...
    'gazedistthresh', 'probes', 'visblocks', 'ICblocks', ...
    'nwbsessions', 'probeneuagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctx_fixedgazeagg', ...
    'meanFRvecagg', 'sponFRvecagg', 'vis', ...
    'ICsigfields', 'ICsigfieldsize2', 'ICsig_fixedgazeagg', ...
    'RFCIfields', 'RFCIfieldsize2', 'RFCI_fixedgazeagg', ...
    'RFCIspinfields', 'RFCIspinfieldsize2', 'RFCIspin_fixedgazeagg', ...
    'sizeCIfields', 'sizeCIfieldsize2', 'sizeCI_fixedgazeagg', ...
    'oriparamsfields', 'oriparamsfieldsize2', 'oriparams_fixedgazeagg', ...
    'ori4paramsfields', 'ori4paramsfieldsize2', 'ori4params_fixedgazeagg', '-v7.3')
copyfile('/Users/hyeyoung/Documents/DATA/openscope_popavg_fixedgazeagg.mat', '/Users/hyeyoung/Documents/DATA/ICexpts_submission22/')

%% aggregate Ron Roff and psthagg
aggpsth = true;
if ~aggpsth
probesR = {'C'};
else
    probesR = probes;
end

vistrialtypes_fixedgazeagg = struct();
vistrialrep_fixedgazeagg = struct();
vistrialorder_fixedgazeagg = struct();
vistrial_fixedgazeagg = struct();
if aggpsth
psthavg_fixedgazeagg = struct();
end
Ronavg_fixedgazeagg = struct();
Roffavg_fixedgazeagg = struct();
for b = 1:numel(visblocks)
%     vistrialtypes.(visblocks{b})=[];
%     vistrialrep.(visblocks{b})=[];
%     vistrialorder.(visblocks{b})=[];
%     vistrialfixedgaze.(visblocks{b})=[];
    for iprobe = 1:numel(probesR)
        if aggpsth
            psthavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe})=[];
        end
        Ronavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe})=[];
        Roffavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe})=[];
    end
end

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
try
    load([pathpp 'trackmouseeye.mat'])
catch
    fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
    continue
end
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probesR)
        tic
        if exist([pathpp, 'probes.mat'], 'file')
            probelist = load([pathpp, 'probes.mat']);
            warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        end
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        if ~ismember(probesR{iprobe}, probelist.probes)
            warning('%s: Probe%s does not exist', nwbsessions{ises}, probesR{iprobe})
            continue
        end

        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probesR{iprobe}))
        % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        %     load(sprintf('%svisresponses_probe%s.mat', pathpp, probesR{iprobe}))
        %     % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
                
        if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
            error('check neuoind')
        end
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        
        fprintf('Probe %s Area %s: %d/%d\n', probesR{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
        disp(unique(probelocs)')
        
%         probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
%         neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
%         neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
%         neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
%         
%         sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
%         sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
        
%         % aggregate cortex neurons
%         neuoi = neuctx;
        % aggregate all neurons
        neuoi = true(size(neuoind));
        
        
        % visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        %     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
        %ICblocks: each stim is 0.4s, inter-trial interval is 0.4s, static images
        %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
        for b = 1:numel(visblocks)
            vistrialtypes_fixedgazeagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
            vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            if strcmp(visblocks{b}, 'RFCI_presentations')
                vistrialtypes_fixedgazeagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder(1:4:end));
                vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder(1:4:end);
            end
            vistrialorder_fixedgazeagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            
            if ismember(visblocks{b}, ICblocks)
                tlon = psthtli>0 & psthtli<=400;
                tloff = psthtli>400 & psthtli<=800;
                temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
            elseif strcmp(visblocks{b}, 'RFCI_presentations')
                tlon = psthtli>0 & psthtli<=1000;
                tloff = psthtli>1000 & psthtli<=1000;
                tempgazedist = trialmaxdistmodecom.RFCI_presentations(1:4:end);
                templikelyblink = triallikelyblink.RFCI_presentations(1:4:end);
                temptrialsfixedgaze = tempgazedist<gazedistthresh & ~templikelyblink;
                temptrialsfixedgaze = reshape(repmat(temptrialsfixedgaze', 4,1),[],1);
            elseif strcmp(visblocks{b}, 'sizeCI_presentations')
                tlon = psthtli>0 & psthtli<=250;
                tloff = psthtli>250 & psthtli<=750;
                temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
            else
                error('vis block not recognized')
            end
            
            vistrial_fixedgazeagg(ises).(visblocks{b}) = temptrialsfixedgaze;
            Ntt = numel(vistrialtypes_fixedgazeagg(ises).(visblocks{b}));
            vistrialrep_fixedgazeagg(ises).(visblocks{b}) = zeros(Ntt,1);
            temppsth = zeros(length(psthtli), Ntt, nnz(neuoi));
            tempRonavg = zeros(Ntt, nnz(neuoi));
            tempRoffavg = zeros(Ntt, nnz(neuoi));
            for ii = 1:Ntt
                trialsoi = temptrialsfixedgaze & vis.(visblocks{b}).trialorder==vistrialtypes_fixedgazeagg(ises).(visblocks{b})(ii);
                vistrialrep_fixedgazeagg(ises).(visblocks{b})(ii) = nnz(trialsoi);
                temppsth(:,ii,:) = mean(1000*psth.(visblocks{b})(:,trialsoi,neuoi), 2);
            tempRonavg(ii,:) = squeeze(mean(1000*psth.(visblocks{b})(tlon,trialsoi,neuoi), [1 2]));
            tempRoffavg(ii,:) = squeeze(mean(1000*psth.(visblocks{b})(tloff,trialsoi,neuoi), [1 2]));
            end
            if aggpsth
                psthavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}) = cat(3, psthavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}), temppsth );
            end
            Ronavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}) = cat(2, Ronavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}), tempRonavg );
            Roffavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}) = cat(2, Roffavg_fixedgazeagg.(visblocks{b}).(probesR{iprobe}), tempRoffavg );
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        toc
    end
    
end

if aggpsth
    save([datadir 'postprocessed\openscope_psthavg_fixedgazeagg.mat'], 'probes', 'visareas', 'visind', 'nwbsessions', ...
    	'probeneuagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctx_fixedgazeagg', ...
        'vistrialtypes_fixedgazeagg', 'vistrialrep_fixedgazeagg', 'vistrialorder_fixedgazeagg', 'vistrial_fixedgazeagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavg_fixedgazeagg', 'Ronavg_fixedgazeagg', 'Roffavg_fixedgazeagg', '-v7.3')
    save('/Users/hyeyoung/Documents/DATA/ICexpts_submission22/openscope_psthavg_fixedgazeagg.mat', 'probes', 'visareas', 'visind', 'nwbsessions', ...
    	'probeneuagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctx_fixedgazeagg', ...
        'vistrialtypes_fixedgazeagg', 'vistrialrep_fixedgazeagg', 'vistrialorder_fixedgazeagg', 'vistrial_fixedgazeagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavg_fixedgazeagg', 'Ronavg_fixedgazeagg', 'Roffavg_fixedgazeagg', '-v7.3')
end

%% MATCH BETWEEN CRF POSITION FOR ALL TRIALS VS FIXED GAZE TRIALS: ~63% neurons show match
% load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'])

whichprobe = 'C';
iprobe = find(strcmp(probes, whichprobe));
neuinarea = strcmp(neulocagg{iprobe}, 'VISp2/3'); % 86% neurons show match
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp'); % 77% neurons show match
% neuinarea = neuctx_fixedgazeagg{iprobe}==1; % 63% neurons show match
% neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

all2fix = ismember(sesneuagg{iprobe}, sesneuagg{iprobe});
RFindclassic = RFCI_fixedgazeagg.(whichprobe).RFindclassic(all2fix);

%neuoi = neuinarea;
neuoi = neuinarea & RFCI_fixedgazeagg.(whichprobe).Pkw_rfclassic(all2fix)<0.05;
%neuoi = neuinarea & RFCIagg.(whichprobe).RFexclsigclassic(all2fix)==1;
figure; histogram2(RFindclassic(neuoi), RFCI_fixedgazeagg.(whichprobe).RFindclassic(neuoi), 'displaystyle', 'tile')

disp('match between CRF position for all trials vs fixed gaze trials')
disp(mean(RFindclassic(neuoi)==RFCI_fixedgazeagg.(whichprobe).RFindclassic(neuoi)))


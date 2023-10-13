% IMPORTANT VERSION UPDATE ON 231010 : aggregate to match the original order
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v230821')

datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

%% 
probevisareas = cell(numel(probes), Nsessions);
for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
                
        % if exist([pathpp, 'probes.mat'], 'file')
        %     probelist = load([pathpp, 'probes.mat']);
        %     warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        % else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        % end
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );

%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
%         disp(probes{iprobe})
%         disp(unique(neuloc(contains(neuloc, 'VIS')))')

        % probevisareas{iprobe, ises} = sprintf('%s ',unique(neuloc(contains(neuloc, 'VIS'))));
        tempcell = unique(neuloc(contains(neuloc, 'VIS')));
        probevisareas{iprobe, ises} = sprintf('%s ',tempcell{:});
    end
end

% sesvisareas = cell(1,Nsessions);
% for ises = 1:Nsessions
%     sesvisareas{ises} = sprintf('%s; ', probevisareas{:,ises});
% end

disp('Sessions with wrong CCF labels (all probes labeled as VISl)')
disp(nwbsessions( all(contains(probevisareas, 'VISl'),1) )')

disp('Sessions with no layer information')
disp(nwbsessions( ~all(contains(probevisareas, 'VISl'),1) & ~any(contains(probevisareas, '2'),1) )')

disp('Sessions with seemingly correct CCF labels & layer info')
disp(nwbsessions( any(contains(probevisareas, '2'),1) )')

% strcmp(probevisareas, 'VISl')

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=1;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_probeC.mat', pathpp))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
allfields = fieldnames(ICsig.ICwcfg1_presentations);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
end
ICsigfields = allfields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};

allfields = fieldnames(RFCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI.(allfields{f}),1)==Nneurons;
end
RFCIfields = allfields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};

allfields = fieldnames(RFCIspin);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin.(allfields{f}),1)==Nneurons;
end
RFCIspinfields = allfields(validfields);

% size vector [0, 4, 8, 16, 32, 64 ]
allfields = fieldnames(sizeCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI.(allfields{f}),1)==Nneurons;
end
sizeCIfields = allfields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};

allfields = fieldnames(oriparams);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams.(allfields{f}),1)==Nneurons;
end
oriparamsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};

allfields = fieldnames(ori4params);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params.(allfields{f}),1)==Nneurons;
end
ori4paramsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};



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
% sesneuctxagg = cell(numel(probes),Nsessions);

meanFRvecagg = cell(numel(probes),Nsessions);
sponFRvecagg = cell(numel(probes),Nsessions);

for b = 1:numel(ICblocks)
    ICsigagg.(ICblocks{b}) = cell(numel(probes), Nsessions);
end
RFCIagg = cell(numel(probes), Nsessions);
RFCIspinagg = cell(numel(probes), Nsessions);
sizeCIagg = cell(numel(probes), Nsessions);
oriparamsagg = cell(numel(probes), Nsessions);
ori4paramsagg = cell(numel(probes), Nsessions);


for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    load([pathpp 'qc_units.mat']) %'unit_amplitude', 'unit_isi_violations', 'unit_wfdur', 'unit_amplitude_cutoff', 'unit_presence_ratio'

unit_amplitude_agg{ises} = unit_amplitude;
unit_isi_violations_agg{ises} = unit_isi_violations;
unit_wfdur_agg{ises} = unit_wfdur;
unit_amplitude_cutoff_agg{ises} = unit_amplitude_cutoff;
unit_presence_ratio_agg{ises} = unit_presence_ratio;

    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        
        probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        
        if nnz(floor(unit_peakch/1000)==probeind-1)==0
            fprintf('Probe %s Area %s: NO UNITS!!!\n', probes{iprobe}, visareas{iprobe} )
            continue
        end
%         tic
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}), 'neuoind')
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
            error('check neuoind')
        end
        % neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
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
        
        neuoindagg{iprobe, ises} = neuoind;
        probeneuagg{iprobe, ises} = iprobe*ones(length(neuoind),1);
        neulocagg{iprobe, ises} = neuloc;
        neupeakchagg{iprobe, ises} = unit_peakch(neuoind);
        neuctxagg{iprobe, ises} = neuctx;
        
        sesneuagg{iprobe, ises} = ises*ones(length(neuoind),1);
        % sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
        
        meanFRvecagg{iprobe, ises} = meanFRvec';
        sponFRvecagg{iprobe, ises} = sponFRvec';
        
        for b = 1:numel(ICblocks)
            ICsigagg.(ICblocks{b}){iprobe,ises} = ICsig.(ICblocks{b});
        end
        
        RFCIagg{iprobe,ises} = RFCI;
        RFCIspinagg{iprobe,ises} = RFCIspin;
        sizeCIagg{iprobe,ises} = sizeCI;
        oriparamsagg{iprobe,ises} = oriparams;
        ori4paramsagg{iprobe,ises} = ori4params;

    end
    
end

Nrfs = size(RFCI.Rrfclassic, 2);
Nszs = size(sizeCI.Rsizeclassic, 2);
Ndirs = size(oriparams.Rori, 2);
Noris = Ndirs/2;

ises=1;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

%%
unit_amplitude_all = cat(1,unit_amplitude_agg{:});
unit_isi_violations_all = cat(1,unit_isi_violations_agg{:});
unit_wfdur_all = cat(1,unit_wfdur_agg{:});
unit_amplitude_cutoff_all = cat(1,unit_amplitude_cutoff_agg{:});
unit_presence_ratio_all = cat(1,unit_presence_ratio_agg{:});

neuoindall = cat(1,neuoindagg{:});
probeneuall = cat(1,probeneuagg{:});
neulocall = cat(1,neulocagg{:});
neupeakchall = cat(1,neupeakchagg{:});
neuctxall = cat(1,neuctxagg{:});
sesneuall = cat(1,sesneuagg{:});
meanFRvecall = cat(1,meanFRvecagg{:});
sponFRvecall = cat(1,sponFRvecagg{:});
% figure; plot(sesneuall*10+probeneuall, '.')

ICsigall = struct();
for b = 1:numel(ICblocks)
    temp = cat(1,ICsigagg.(ICblocks{b}){:});
    for f= 1:numel(ICsigfields)
        ICsigall.(ICblocks{b}).(ICsigfields{f}) = cat(1,temp.(ICsigfields{f}));
    end
end

RFCIall = struct();
for f= 1:numel(RFCIfields)
    temp = cat(1,RFCIagg{:});
    RFCIall.(RFCIfields{f}) = cat(1,temp.(RFCIfields{f}));
end

RFCIspinall = struct();
for f= 1:numel(RFCIspinfields)
    temp = cat(1,RFCIspinagg{:});
    RFCIspinall.(RFCIspinfields{f}) = cat(1,temp.(RFCIspinfields{f}));
end

sizeCIall = struct();
for f= 1:numel(sizeCIfields)
    temp = cat(1,sizeCIagg{:});
    sizeCIall.(sizeCIfields{f}) = cat(1,temp.(sizeCIfields{f}));
end

oriparamsall = struct();
for f= 1:numel(oriparamsfields)
    temp = cat(1,oriparamsagg{:});
    oriparamsall.(oriparamsfields{f}) = cat(1,temp.(oriparamsfields{f}));
end

ori4paramsall = struct();
for f= 1:numel(ori4paramsfields)
    temp = cat(1,ori4paramsagg{:});
    ori4paramsall.(ori4paramsfields{f}) = cat(1,temp.(ori4paramsfields{f}));
end

clearvars temp

%% sanity check: neuoind is ordered according to the probe order
for ises = 1:Nsessions
    if ~isequal( cat(1,neuoindagg{:,ises}), (1:length(cat(1,neuoindagg{:,ises})))' )
        error('neuoind is not ordered according to probe indices')
    end
end

if ~isequal(floor(neupeakchall/1000), probeneuall-1)
    error('check probeneuall')
end

%%
save(['S:\OpenScopeData\00248_v230821\postprocessed\openscope_popavg_all.mat'], ...
    'nwbsessions', 'neuoindall', 'probeneuall', 'neulocall', 'neupeakchall', 'sesneuall', 'neuctxall', ...
    'meanFRvecall', 'sponFRvecall', 'vis', ...
    'unit_amplitude_all', 'unit_isi_violations_all', 'unit_wfdur_all', 'unit_amplitude_cutoff_all', 'unit_presence_ratio_all', ...
    'ICsigfields', 'ICsigall', 'RFCIfields', 'RFCIall', ...
    'RFCIspinfields', 'RFCIspinall', 'sizeCIfields', 'sizeCIall', ...
    'oriparamsfields', 'oriparamsall', 'ori4paramsfields', 'ori4paramsall', '-v7.3')

%% aggregate Ron Roff and psthagg
aggpsth = true;
if ~aggpsth
probesR = {'C'};
else
    probesR = probes;
end

load([pathpp 'postprocessed_probeC.mat'], 'psthtli', 'vis')
Nneuronsall = numel(sesneuall);

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
vistrialtypesagg = struct();
vistrialrepagg = struct();
vistrialorderagg = struct();
if aggpsth
psthavgall = struct();
psthsemall = struct();
end
Ronavgall = struct();
Roffavgall = struct();
Ronstdall = struct();
Roffstdall = struct();
for b = 1:numel(visblocks)
    temptt = unique(vis.(visblocks{b}).trialorder);
    if strcmp(visblocks{b}, 'RFCI_presentations')
        temptt = unique(vis.(visblocks{b}).trialorder(1:4:end));
    end
    Ntt = numel(temptt);
    for iprobe = 1:numel(probesR)
        if aggpsth
            psthavgall.(visblocks{b}) = NaN(length(psthtli), Ntt, Nneuronsall);
            psthsemall.(visblocks{b}) = NaN(length(psthtli), Ntt, Nneuronsall);
        end
        Ronavgall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Roffavgall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Ronstdall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
        Roffstdall.(visblocks{b}) = NaN(Nneuronsall,Ntt);
    end
end

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probesR)
        tic
        probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        probeind = find( strcmp(probesR{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        if ~ismember(probesR{iprobe}, probelist.probes)
            warning('%s: Probe%s does not exist', nwbsessions{ises}, probesR{iprobe})
            continue
        end
           
        tempneusesprobe = sesneuall==ises & probeneuall==probeind;

        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probesR{iprobe}))
        % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        %     load(sprintf('%svisresponses_probe%s.mat', pathpp, probesR{iprobe}))
        %     % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
            error('check neuoind')
        end
        if numel(neuoind) ~= nnz(tempneusesprobe)
            error('check neuoind/sesneuall/probeneuall')
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
        
        
        % visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        %     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
        %ICblocks: each stim is 0.4s, inter-trial interval is 0.4s, static images
        %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
        for b = 1:numel(visblocks)
            vistrialtypesagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
            vistrialorderagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            if strcmp(visblocks{b}, 'RFCI_presentations')
                vistrialtypesagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder(1:4:end));
                vistrialorderagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder(1:4:end);
            end
            
            if ismember(visblocks{b}, ICblocks)
                tlon = psthtli>0 & psthtli<=400;
                tloff = psthtli>400 & psthtli<=800;
            elseif strcmp(visblocks{b}, 'RFCI_presentations')
                tlon = psthtli>0 & psthtli<=1000;
                tloff = psthtli>1000 & psthtli<=1000;
            elseif strcmp(visblocks{b}, 'sizeCI_presentations')
                tlon = psthtli>0 & psthtli<=250;
                tloff = psthtli>250 & psthtli<=750;
            else
                error('vis block not recognized')
            end
            
            Ntt = numel(vistrialtypesagg(ises).(visblocks{b}));
            vistrialrepagg(ises).(visblocks{b}) = zeros(Ntt,1);

            temppsth = NaN(length(psthtli), Ntt, nnz(tempneusesprobe));
            tempRonavg = NaN(nnz(tempneusesprobe), Ntt);
            tempRoffavg = NaN(nnz(tempneusesprobe), Ntt);

            temppsthsem = NaN(length(psthtli), Ntt, nnz(tempneusesprobe));
            tempRonstd = NaN(nnz(tempneusesprobe), Ntt);
            tempRoffstd = NaN(nnz(tempneusesprobe), Ntt);
            for ii = 1:Ntt
                trialsoi = vis.(visblocks{b}).trialorder==vistrialtypesagg(ises).(visblocks{b})(ii);
                vistrialrepagg(ises).(visblocks{b})(ii) = nnz(trialsoi);

                temppsth(:,ii,:) = mean(1000*psth.(visblocks{b})(:,trialsoi,:), 2);
                tempRonavg(:,ii) = squeeze(mean(1000*psth.(visblocks{b})(tlon,trialsoi,:), [1 2]));
                tempRoffavg(:,ii) = squeeze(mean(1000*psth.(visblocks{b})(tloff,trialsoi,:), [1 2]));

                temppsthsem(:,ii,:) = std(1000*psth.(visblocks{b})(:,trialsoi,:), 0,2)/sqrt(nnz(trialsoi));
                tempRonstd(:,ii) = squeeze( std( mean(1000*psth.(visblocks{b})(tlon,trialsoi,:),1), 0,2) );
                tempRoffstd(:,ii) = squeeze( std( mean(1000*psth.(visblocks{b})(tloff,trialsoi,:),1), 0,2) );
            end
            if aggpsth
                psthavgall.(visblocks{b})(:,:,tempneusesprobe) = temppsth;
                psthsemall.(visblocks{b})(:,:,tempneusesprobe) = temppsthsem;
            end
            Ronavgall.(visblocks{b})(tempneusesprobe,:) = tempRonavg ;
            Roffavgall.(visblocks{b})(tempneusesprobe,:) = tempRoffavg ;
            Ronstdall.(visblocks{b})(tempneusesprobe,:) = tempRonstd ;
            Roffstdall.(visblocks{b})(tempneusesprobe,:) = tempRoffstd ;
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        toc
    end
    
end

%%%%%%%%%% TODO: ADD SANITY CHECK THAT THERE ARE NO NAN VALUES %%%%%%%%%% 
for b = 5:numel(visblocks)
    tic
    n1 = nnz(isnan(psthavgall.(visblocks{b})));
    n2 = nnz(isnan(psthsemall.(visblocks{b})));
    n3 = nnz(isnan(Ronavgall.(visblocks{b})));
    n4 = nnz(isnan(Roffavgall.(visblocks{b})));
    n5 = nnz(isnan(Ronstdall.(visblocks{b})));
    n6 = nnz(isnan(Roffstdall.(visblocks{b})));
    if strcmp(visblocks{b}, 'RFCI_presentations')
        if ~all([n1 n2 n3 n5]==0)
            disp([n1 n2 n3 n5])
            error('%s NaN value found -- check', visblocks{b})
        end
    else
        if ~all([n1 n2 n3 n4 n5 n6]==0)
            disp([n1 n2 n3 n4 n5 n6])
            error('%s NaN value found -- check', visblocks{b})
        end
    end
    toc
end

if aggpsth
    save('S:\OpenScopeData\00248_v230821\postprocessed\openscope_psthavgall.mat', ...
        'probes', 'visareas', 'visind', 'nwbsessions', ...
        'neuoindall', 'probeneuall', 'neulocall', 'neupeakchall', 'sesneuall', 'neuctxall', ...
        'unit_amplitude_all', 'unit_isi_violations_all', 'unit_wfdur_all', 'unit_amplitude_cutoff_all', 'unit_presence_ratio_all', ...
        'vistrialtypesagg', 'vistrialrepagg', 'vistrialorderagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavgall', 'Ronavgall', 'Roffavgall', ...
        'psthsemall', 'Ronstdall', 'Roffstdall', '-v7.3')
end


%% report number of units in each area/session/probe
aggneuloc = cat(1,neulocagg{:});
aggsesneu = cat(1,sesneuagg{:});

[v,c]=uniquecnt(aggneuloc);
disp([v c])

areaunitsperses = zeros(length(v), Nsessions);
for ii = 1:length(v)
    neuoi = strcmp(aggneuloc, v(ii));
    hc=histcounts(aggsesneu(neuoi), 0.5:1:Nsessions+0.5);
    if ~( nnz(neuoi)==c(ii) && sum(hc)==c(ii) )
        error('mismatch between uniquecnt and strcmp -- check')
    end
    areaunitsperses(ii,:) = hc;
end

areaunitstab = table(v,c,areaunitsperses);
open areaunitstab


% % 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes,
% % 10-100s: which RFcenter, 1s: which direction
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
%szvec = [0, 4, 8, 16, 32, 64];

% 'gkcellind', 'kSLtunedneurons', 'kSLprefori', 'kSGtunedneurons', 'kSGprefori'

% ori actually means dir
function [sizeCI, oriparams, ori4params] = analyzesizeCI(R, trialorder)

Nneurons = size(R,1);

% determine tuning with fullscreen grating
% fulltrials = floor(trialorder/1000)==11;
% to increase trial repeats, also
% TODO: ADD Pkw_ori
Noris = 8;
Rori = zeros(Nneurons, Noris);
validoriparams = true;
for iori = 1:Noris
    if nnz(trialorder==11000+iori)==0
        validoriparams = false;
    end
    Rori(:,iori) = mean(R(:,trialorder==11000+iori),2);
end

if validoriparams
    [Rpref, prefiori] = max(Rori,[],2);
    if ~isequal(Rpref, Rori(sub2ind(size(Rori), (1:Nneurons)', prefiori)))
        error('check Rpref/prefiori')
    end
    
    orthiori = zeros(Nneurons,1);
    Rorth = NaN(Nneurons,1);
    for iori = 1:Noris
        iorth0 = mod(iori+2 -1, Noris)+1;
        iorth1 = mod(iori-2 -1, Noris)+1;
        neuoi0 = prefiori==iori & Rori(:,iorth0)<=Rori(:,iorth1);
        orthiori(neuoi0) = iorth0;
        Rorth(neuoi0) = Rori(neuoi0, iorth0);
        neuoi1 = prefiori==iori & Rori(:,iorth0)>Rori(:,iorth1);
        orthiori(neuoi1) = iorth1;
        Rorth(neuoi1) = Rori(neuoi1, iorth1);
    end
    if nnz(isnan(Rorth)) || nnz(orthiori==0)
        error('check orthiori and Rorth')
    end
    if ~isequal(Rorth, Rori(sub2ind(size(Rori), (1:Nneurons)', orthiori)))
        error('check Rorth/orthiori')
    end
    
    OSI = (Rpref-Rorth)./(Rpref+Rorth);
    
    % %%
    trialiori = mod(trialorder,100);
    
    OP = NaN(Nneurons,1);
    Pmww_OP = NaN(Nneurons, 1);
    Pkw_ori = NaN(Nneurons, 1);
    for ci = 1:Nneurons
        fgtrials = floor(trialorder/1000)==11;
        preforitrials = fgtrials & trialiori==prefiori(ci);
        orthoritrials = fgtrials & trialiori==orthiori(ci);
        
        scores = [R(ci, orthoritrials) R(ci, preforitrials)];
        labels = [zeros(1,nnz(orthoritrials)) ones(1,nnz(preforitrials))];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        OP(ci) = AUC;
        
        Pmww_OP(ci) = ranksum(R(ci, orthoritrials), R(ci, preforitrials));
        
        [p,tbl,stats] = kruskalwallis(R(ci,fgtrials), trialorder(fgtrials), 'off');
        Pkw_ori(ci) = p;
    end
end

Rori678 = zeros(Nneurons, Noris);
validori678params = true;
for iori = 1:Noris
    oritrials = ismember(trialorder, [6000+iori 11000+iori 12000+iori]);
    if nnz(oritrials)==0
        validori678params = false;
    end
    Rori678(:,iori) = mean(R(:,oritrials),2);
end

if validori678params
    [Rpref678, prefiori678] = max(Rori678,[],2);
    if ~isequal(Rpref678, Rori678(sub2ind(size(Rori678), (1:Nneurons)', prefiori678)))
        error('check Rpref/prefiori')
    end
    
    orthiori678 = zeros(Nneurons,1);
    Rorth678 = NaN(Nneurons,1);
    for iori = 1:Noris
        iorth0 = mod(iori+2 -1, Noris)+1;
        iorth1 = mod(iori-2 -1, Noris)+1;
        neuoi0 = prefiori678==iori & Rori678(:,iorth0)<=Rori678(:,iorth1);
        orthiori678(neuoi0) = iorth0;
        Rorth678(neuoi0) = Rori678(neuoi0, iorth0);
        neuoi1 = prefiori678==iori & Rori678(:,iorth0)>Rori678(:,iorth1);
        orthiori678(neuoi1) = iorth1;
        Rorth678(neuoi1) = Rori678(neuoi1, iorth1);
    end
    if nnz(isnan(Rorth678)) || nnz(orthiori678==0)
        error('check orthiori and Rorth')
    end
    if ~isequal(Rorth678, Rori678(sub2ind(size(Rori678), (1:Nneurons)', orthiori678)))
        error('check Rorth/orthiori')
    end
    
    OSI678 = (Rpref678-Rorth678)./(Rpref678+Rorth678);
    
    % %%
    trialiori = mod(trialorder,100);
    
    OP678 = NaN(Nneurons,1);
    Pmww_OP678 = NaN(Nneurons, 1);
    Pkw_ori678 = NaN(Nneurons, 1);
    for ci = 1:Nneurons
        trials678 = ismember(floor(trialorder/1000), [6 11 12]);
        preforitrials = trials678 & trialiori==prefiori678(ci);
        orthoritrials = trials678 & trialiori==orthiori678(ci);
        
        scores = [R(ci, orthoritrials) R(ci, preforitrials)];
        labels = [zeros(1,nnz(orthoritrials)) ones(1,nnz(preforitrials))];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        OP678(ci) = AUC;
        
        Pmww_OP678(ci) = ranksum(R(ci, orthoritrials), R(ci, preforitrials));
        
        [p,tbl,stats] = kruskalwallis(R(ci,trials678), trialorder(trials678), 'off');
        Pkw_ori678(ci) = p;
    end
end

oriparams = struct();
if validoriparams
    orifields = {'Rori', 'prefiori', 'orthiori', 'OSI', 'OP', 'Pmww_OP', 'Pkw_ori'};
    for f = 1:numel(orifields)
        oriparams.(orifields{f}) = eval(orifields{f});
    end
else
    warning('full grating oriparams skipped because not all orientations were measured')
end

if validori678params
    orifields = {'Rori678', 'prefiori678', 'orthiori678', 'OSI678', 'OP678', 'Pmww_OP678', 'Pkw_ori678'};
    for f = 1:numel(orifields)
        oriparams.(orifields{f}) = eval(orifields{f});
    end
else
    warning('CG64, IG0, IG4 oriparams skipped because not all orientations were measured')
end

%% ori4params
Nori4 = 4;
Rori4 = zeros(Nneurons, Nori4);
validori4params = true;
for iori = 1:Nori4
    oritrials = ismember(trialorder, [11000+iori 11004+iori]);
    if nnz(oritrials)==0
        validori4params = false;
    end
    Rori4(:,iori) = mean(R(:,oritrials),2);
end

if validori4params
    [Rpref4, prefiori4] = max(Rori4,[],2);
    if ~isequal(Rpref4, Rori4(sub2ind(size(Rori4), (1:Nneurons)', prefiori4)))
        error('check Rpref4/prefiori4')
    end
    
    orthiori4 = zeros(Nneurons,1);
    Rorth4 = NaN(Nneurons,1);
    for iori = 1:Nori4
        iorth = mod(iori+2 -1, Nori4)+1;
        neuoi = prefiori4==iori;
        orthiori4(neuoi) = iorth;
        Rorth4(neuoi) = Rori4(neuoi, iorth);
    end
    if nnz(isnan(Rorth4)) || nnz(orthiori4==0)
        error('check orthiori and Rorth')
    end
    if ~isequal(Rorth4, Rori4(sub2ind(size(Rori4), (1:Nneurons)', orthiori4)))
        error('check Rorth/orthiori')
    end
    
    OSI4 = (Rpref4-Rorth4)./(Rpref4+Rorth4);
    
    % %%
    trialiori4 = mod( mod(trialorder,100)-1, Nori4 )+1;
    
    OP4 = NaN(Nneurons,1);
    Pmww_OP4 = NaN(Nneurons, 1);
    Pkw_ori4 = NaN(Nneurons, 1);
    for ci = 1:Nneurons
        fgtrials = floor(trialorder/1000)==11;
        preforitrials = fgtrials & trialiori4==prefiori4(ci);
        orthoritrials = fgtrials & trialiori4==orthiori4(ci);
        
        scores = [R(ci, orthoritrials) R(ci, preforitrials)];
        labels = [zeros(1,nnz(orthoritrials)) ones(1,nnz(preforitrials))];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        OP4(ci) = AUC;
        
        Pmww_OP4(ci) = ranksum(R(ci, orthoritrials), R(ci, preforitrials));
        
        [p,tbl,stats] = kruskalwallis(R(ci,fgtrials), trialorder(fgtrials), 'off');
        Pkw_ori4(ci) = p;
    end
end

Rori4678 = zeros(Nneurons, Nori4);
validori4678params = true;
for iori = 1:Nori4
    oritrials = ismember(trialorder, [6000+iori 11000+iori 12000+iori 6004+iori 11004+iori 12004+iori]);
    if nnz(oritrials)==0
        validori4678params = false;
    end
    Rori4678(:,iori) = mean(R(:,oritrials),2);
end

if validori4678params
    [Rpref4678, prefiori4678] = max(Rori4678,[],2);
    if ~isequal(Rpref4678, Rori4678(sub2ind(size(Rori4678), (1:Nneurons)', prefiori4678)))
        error('check Rpref/prefiori')
    end
    
    orthiori4678 = zeros(Nneurons,1);
    Rorth4678 = NaN(Nneurons,1);
    for iori = 1:Nori4
        iorth = mod(iori+2 -1, Nori4)+1;
        neuoi = prefiori4678==iori;
        orthiori4678(neuoi) = iorth;
        Rorth4678(neuoi) = Rori4678(neuoi, iorth);
    end
    if nnz(isnan(Rorth4678)) || nnz(orthiori4678==0)
        error('check orthiori and Rorth')
    end
    if ~isequal(Rorth4678, Rori4678(sub2ind(size(Rori4678), (1:Nneurons)', orthiori4678)))
        error('check Rorth/orthiori')
    end
    
    OSI4678 = (Rpref4678-Rorth4678)./(Rpref4678+Rorth4678);
    
    % %%
    trialiori4 = mod( mod(trialorder,100)-1, Nori4 )+1;
    
    OP4678 = NaN(Nneurons,1);
    Pmww_OP4678 = NaN(Nneurons, 1);
    Pkw_ori4678 = NaN(Nneurons, 1);
    for ci = 1:Nneurons
        trials4678 = ismember(floor(trialorder/1000), [6 11 12]);
        preforitrials = trials4678 & trialiori4==prefiori4678(ci);
        orthoritrials = trials4678 & trialiori4==orthiori4678(ci);
        
        scores = [R(ci, orthoritrials) R(ci, preforitrials)];
        labels = [zeros(1,nnz(orthoritrials)) ones(1,nnz(preforitrials))];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        OP4678(ci) = AUC;
        
        Pmww_OP4678(ci) = ranksum(R(ci, orthoritrials), R(ci, preforitrials));
        
        [p,tbl,stats] = kruskalwallis(R(ci,trials4678), trialorder(trials4678), 'off');
        Pkw_ori4678(ci) = p;
    end
end

ori4params = struct();
if validori4params
    ori4fields = {'Rori4', 'prefiori4', 'orthiori4', 'OSI4', 'OP4', 'Pmww_OP4', 'Pkw_ori4'};
    for f = 1:numel(ori4fields)
        ori4params.(ori4fields{f}) = eval(ori4fields{f});
    end
else
    warning('full grating ori4params skipped because not all orientations were measured')
end

if validori4678params
    ori4fields = {'Rori4678', 'prefiori4678', 'orthiori4678', 'OSI4678', 'OP4678', 'Pmww_OP4678', 'Pkw_ori4678'};
    for f = 1:numel(ori4fields)
        ori4params.(ori4fields{f}) = eval(ori4fields{f});
    end
else
    warning('CG64, IG0, IG4 ori4params skipped because not all orientations were measured')
end

%%
Nsz = 6;
Rsizeclassic = NaN(Nneurons, Nsz);
Rsizeinverse = NaN(Nneurons, Nsz);
validsizeclassic = true;
validsizeinverse = true;
% 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
for szi = 1:Nsz
    csztrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 10000) / 1000) == szi ;
    if nnz(csztrials)==0
        validsizeclassic = false;
    end
    Rsizeclassic(:,szi) = mean(R(:, csztrials),2);
    
    isztrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 10000) / 1000) == szi ;
    if nnz(csztrials)==0
        validsizeinverse = false;
    end
    Rsizeinverse(:,szi) = mean(R(:, isztrials),2);
end

[~, sizeindclassic] = max(Rsizeclassic, [], 2);
[~, sizeindinverse] = max(Rsizeinverse, [], 2);

indRFCIsize = floor(mod(trialorder, 10000)/1000);
indRFCIloc = floor(mod(trialorder, 1000)/10);

% for each neuron, choose the preferred size
Pkw_sizeclassic = zeros(Nneurons, 1);
Pkw_sizeinverse = zeros(Nneurons, 1);
for ci = 1:Nneurons
    csztrials =  trialorder~=0 & floor(trialorder/10000) == 0;
    [p,tbl,stats] = kruskalwallis(R(ci,csztrials), indRFCIsize(csztrials), 'off');
    Pkw_sizeclassic(ci) = p;
    
    isztrials =  trialorder~=0 & floor(trialorder/10000) == 1;
    [p,tbl,stats] = kruskalwallis(R(ci,isztrials), indRFCIsize(isztrials), 'off');
    Pkw_sizeinverse(ci) = p;
end

sizeCIfields = {'validsizeclassic', 'validsizeinverse', 'Rsizeclassic', 'Rsizeinverse', 'sizeindclassic', 'sizeindinverse', ...
    'Pkw_sizeclassic', 'Pkw_sizeinverse'};%, 'pRsizeclassic', 'pRsizeinverse', 'psizeclassic', 'psizeinverse'};
sizeCI = struct();
for f = 1:numel(sizeCIfields)
    sizeCI.(sizeCIfields{f}) = eval(sizeCIfields{f});
end


end
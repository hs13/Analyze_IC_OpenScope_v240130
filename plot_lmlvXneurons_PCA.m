if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
else
    drivepath = 'G:/My Drive/';
end
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])

neuopt = 'V1';

switch neuopt
    case 'V1'
        spkcntICttagg = spkcntIChiV1agg;
    case 'visctx'
        spkcntICttagg = spkcntIChivisctxagg;
end

Nsessions = numel(spkcntICttagg);
Ntt = numel(spkcntICttagg{1});
spkcntICttacc = cat(2, spkcntICttagg{:})';

csz = cellfun(@size, spkcntICttacc(:,1), 'UniformOutput', false);
csz = cat(1,csz{:,1});
Nneuronsvec = csz(:,2);
Nneurons = sum(csz(:,2));

spkcntICttaccavg = cellfun(@mean, spkcntICttacc, 'uniformoutput', false);
spkcntICttaccvar = cellfun(@var, spkcntICttacc, 'uniformoutput', false);

spkcntICttavg = NaN(Nneurons, Ntt);
spkcntICttvar = NaN(Nneurons, Ntt);
for ii = 1:Ntt
    spkcntICttavg(:,ii) = cat(2,spkcntICttaccavg{:,ii});
    spkcntICttvar(:,ii) = cat(2,spkcntICttaccvar{:,ii});
end

clearvars spkcntICttacc spkcntICttaccavg spkcntICttaccvar

%% example session
ises = 21;
try
    tempspkcnt = cat(3,spkcntnsttagg{ises}{:});
    Nrep = size(tempspkcnt,1);
    Nneu = size(tempspkcnt,2);
catch
    % trial repetitions were not the same across natural scene trial types
    csz = cellfun(@size, spkcntnsttagg{ises}, 'UniformOutput', false);
    csz = cat(1,csz{:});
    Nrep = min(csz(:,1));
    if all(csz(:,2)==csz(1,2))
        Nneu = csz(1,2);
    else
        error('check number of neurons in session %d', ises)
    end
    tempspkcnt = NaN(Nrep, Nneu, Ntt);
    for n = 1:Ntt
        tempspkcnt(:,:,n) = spkcntnsttagg{ises}{n}(1:Nrep,:);
    end
end
spkcntnstt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu

%%
repsimfn = [drivepath 'RESEARCH/logmean_logvar/naturalscenes_representationsimilarity_', neuopt, '.mat'];
load(repsimfn, 'lmlvslope', 'lmlvyintercept')
disperses = -1:0.5:2;
movetoline = false;

EVdisperses = NaN(numel(disperses), Nneu);
EVttdisperses = NaN(Ntt, numel(disperses), Nneu);
for ii = 1:numel(disperses)
    tic
    spkcntres = spkcntnstt - mean(spkcntnstt,1); % Nrep * Nnstt * Nneu
    spkcntmu = mean(spkcntnstt,1); % 1XNimg X Nneurons
    spkcntvar = var(spkcntnstt,0,1); % 1XNimg X Nneurons
    temp = spkcntvar; temp(spkcntvar==0)=NaN;
    totvar = nanmean(temp,2);

    tempx = log10(spkcntmu);
    tempx(spkcntmu==0) = NaN;
    meanx = squeeze(nanmean(tempx,3)); % average across neurons: 1XNimg

    % when original log(mean) vs log(var) relationship is y=ax+b,
    % and new relationship is y=cx+d
    % d = (a-c)*mean(x)+b —> this keeps mean(y) constant
    % original variance: v0
    % new variance: v1
    % v1 = 10^(c/a * (log10(v0)-b) +d);
    % rescaled residual for every image and every neuron.
    % noise correlation does not change
    Avec = lmlvslope(ises,:);
    Bvec = lmlvyintercept(ises,:);
    Cvec = disperses(ii)*ones(1,Ntt);
    Dvec = (Avec-Cvec).*meanx + Bvec;

    newspkcntvar = 10.^( (Cvec./Avec).*(log10(spkcntvar)-Bvec) + Dvec);
    newspkcntres = spkcntres .* sqrt(newspkcntvar./spkcntvar);
    newspkcntnstt = mean(spkcntnstt,1)+newspkcntres;

    if movetoline
        spkcntvarlin = 10.^(lmlvslope(ises,:).*log10(spkcntmu)+ lmlvyintercept(ises,:) );
        spkcntreslin = spkcntres .* sqrt(spkcntvarlin./spkcntvar);
        spkcntnstt = spkcntmu+spkcntreslin;

        newspkcntvarlin = 10.^( (Cvec./Avec).*(log10(spkcntvarlin)-Bvec) + Dvec);
        newspkcntreslin = spkcntreslin .* sqrt(newspkcntvarlin./spkcntvarlin);
        newspkcntnstt = spkcntmu+newspkcntreslin;
    end

    [newcoeff, newscore, newlatent, newtsquared, newexplained, newmu] =pca(reshape(newspkcntnstt,Nrep*Ntt,Nneu));
    %figure; plot(log10(1:length(newexplained)), log10(newexplained), '.')
    EVdisperses(ii,1:length(newexplained)) = newexplained;

    for itt = 1:Ntt
    [newcoeff, newscore, newlatent, newtsquared, newexplained, newmu] =pca(squeeze(newspkcntnstt(:,itt,:)));
    EVttdisperses(itt,ii,1:length(newexplained)) = newexplained;
    end
end

figure;
plot(log10(1:Nneu), log10(EVdisperses))
colormap copper
legend

figure
plot(disperses, log10(EVdisperses(:,1))-log10(EVdisperses(:,10)), '.-')

xl = [disperses(1) disperses(end)];
figure
hold all
plot(disperses, squeeze(log10(EVttdisperses(:,:,1))-log10(EVttdisperses(:,:,10))), '.-')
plot(xl,[1 1], 'k-')
ylabel('PCA alpha')
xlabel('log(mean) log(var) slope')
title('Allen Brain Natural Scenes: PCA on each trial type')

% it appears that within each trial type, PCA alpha is equal to log(mean)
% log(var) slope!!!

figure; plot(squeeze(log10(mean(newspkcntnstt,1)))', squeeze(log10(var(newspkcntnstt,0,1)))', '.')
figure; plot(squeeze(log10(mean(spkcntnstt,1)))', squeeze(log10(var(spkcntnstt,0,1)))', '.')

%% fractal dimensions: box-counting
tempspkcntall = reshape(spkcntnstt, Nrep*Ntt, Nneu);

spkcntbins = 1:max(tempspkcntall(:));
Nvoxelsbins = zeros(1, numel(spkcntbins));
rankbins = zeros(1, numel(spkcntbins));
maxNneighborsbins= zeros(1, numel(spkcntbins));
for ibin = 1:numel(spkcntbins)
    tempspkbinned = floor(tempspkcntall/(spkcntbins(ibin)));
    uniqvecs = unique(tempspkbinned, 'rows');
    rankbins(ibin) = rank(tempspkbinned);
    % exclude the ones that are surrounded by occupied voxels (pick out
    % only the boudnary). max number of surrounding voxels is 2*Nneu
    numvoxelneighbors = zeros(size(uniqvecs,1),1);
    for ci = 1:Nneu
        uniqvecsplus1 = uniqvecs;
        uniqvecsplus1(:,ci) = uniqvecs(:,ci)+1;
        A = ismember(uniqvecsplus1, uniqvecs, 'rows');
        numvoxelneighbors = numvoxelneighbors+A;
        
        uniqvecsminus1 = uniqvecs;
        uniqvecsminus1(:,ci) = uniqvecs(:,ci)-1;
        B = ismember(uniqvecsminus1, uniqvecs, 'rows');
        numvoxelneighbors = numvoxelneighbors+B;
    end    
    maxNneighborsbins(ibin) = max(numvoxelneighbors);
    Nvoxelsbins(ibin) = size(uniqvecs,1);
end

figure;hold all
plot(rankbins, log10(Nvoxelsbins)./log10(spkcntbins))
plot(rankbins, rankbins, 'r-')

figure;plot(spkcntbins, maxNneighborsbins)
% % sanity check
% plus1neighbors = 0;
% for irow = 1:size(uniqvecsplus1)
%     plus1neighbors = plus1neighbors + ismember(uniqvecsplus1(irow,:), uniqvecs, 'rows');
% end

tempnewspkcntall = reshape(newspkcntnstt, Nrep*Ntt, Nneu);
tempnewspkcntall(tempnewspkcntall<0) = 0; 
% get rid of trials with NaN values
invaltrials = any(isnan(tempnewspkcntall),2);
tempnewspkcntall(invaltrials,:) = [];
newspkcntbins = 1:max(tempnewspkcntall(:));
newNvoxelsbins = zeros(size(newspkcntbins));
for ibin = 1:numel(newspkcntbins)
    tempnewspkbinned = floor(tempnewspkcntall/(newspkcntbins(ibin)));
    uniqvecs = unique(tempnewspkbinned, 'rows');
    newNvoxelsbins(ibin) = size(uniqvecs,1);
end

%figure; plot(spkcntbins, Nvoxelsbins, '.-')
figure; hold all
plot(log10(spkcntbins), log10(Nvoxelsbins), 'b.-')
plot(log10(newspkcntbins), log10(newNvoxelsbins), 'r.-')
% newslopes = diff(log10(newNvoxelsbins))./diff(log10(newspkcntbins));
% plot(log10(newspkcntbins(2:end)), newslopes, 'r-')

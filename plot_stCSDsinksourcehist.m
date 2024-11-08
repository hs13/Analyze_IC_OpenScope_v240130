

%% pick out the location of the strongest sink and strongest source in each area


%% spike triggered CSD and TFR
ctxelec = contains(lfpelecvec.location, 'VIS');
ctxelectop = find(ctxelec, 1, 'last');
ctxelecbottom = find(ctxelec, 1, 'first');
ctxcsdelectinds = find(ismember(csdelectinds, ctxelecbottom:ctxelectop));
eord = (1:numel(ctxcsdelectinds))';
ctxcsdblock = csdblock(ctxcsdelectinds,:);

trange = -50:0;
spiketraintrimmed = spiketrain;
spiketraintrimmed(1:-trange(1),:)=false;
spiketraintrimmed(end-trange(end)+1:end,:)=false;
stCSDsinkhist = zeros(Nneurons, numel(ctxcsdelectinds));
stCSDsourcehist = zeros(Nneurons, numel(ctxcsdelectinds));
% ~2000 units per session, 7 seconds per neuron, ~1 day per session. too long. choose neurons to do this with.
% ~200 V1 units per session, 7 seconds per neuron, 140 min per session.
tic
for ii = 1:Nneurons
    % tic
    tempspkinds = find(spiketraintrimmed(:,ii));
    tempNspikes = numel(tempspkinds);

    Ncutspks = 10000;
    if numel(tempspkinds)>Ncutspks % randomly pick Ncutspks so that code does not run out of memory
        cutspks = randperm(tempNspikes, Ncutspks);
    else
        cutspks = randi(tempNspikes, [1,Ncutspks]);
    end
    tempspkinds = tempspkinds(cutspks);

    tempstinds = tempspkinds+trange;

    tempstcsd = NaN( Ncutspks, numel(ctxcsdelectinds) );
    for e = 1:numel(ctxcsdelectinds)
        tempcsd = ctxcsdblock(e,:);
        tempstcsd(:,e) = mean(tempcsd(tempstinds),2);
    end

    [mv,mi]=max(tempstcsd,[],2);
    [v,c]=uniquecnt(mi);
    % if ~isequal(eord(ismember(eord,v)), v)
    %     error('expected uniquecnt to return values in order')
    % end
    stCSDsourcehist(ii,ismember(eord,v)) = c/Ncutspks;

    [mv,mi]=min(tempstcsd,[],2);
    [v,c]=uniquecnt(mi);
    stCSDsinkhist(ii,ismember(eord,v)) = c/Ncutspks;
    % toc
end
toc

%%
segall = ICsigall.ICwcfg1_presentations.indin1 | ICsigall.ICwcfg1_presentations.indin2 | ...
    ICsigall.ICwcfg1_presentations.indin3 | ICsigall.ICwcfg1_presentations.indin4;
V1segresp = find(segall(neuindV1));
V1ICenc = find(ICsigall.ICwcfg1_presentations.ICencoder(neuindV1));

figure; 
subplot(2,2,1)
imagesc(stCSDsourcehist(V1ICenc,:))
colorbar
subplot(2,2,2)
imagesc(stCSDsinkhist(V1ICenc,:))
colorbar
subplot(2,2,3)
imagesc(stCSDsourcehist(V1segresp,:))
colorbar
subplot(2,2,4)
imagesc(stCSDsinkhist(V1segresp,:))
colorbar

figure
subplot(1,2,1); hold all
plot(ctxcsdelectinds, stCSDsourcehist(V1segresp,:), 'Color', [0.7 0 0.7])
plot(ctxcsdelectinds, stCSDsourcehist(V1ICenc,:), 'Color', [0 0.7 0])
set(gca, 'XTick', ctxcsdelectinds, 'XTickLabel', lfpelecvec.location(ctxcsdelectinds))
title('source')

subplot(1,2,2); hold all
plot(ctxcsdelectinds, stCSDsinkhist(V1segresp,:), 'Color', [0.7 0 0.7])
plot(ctxcsdelectinds, stCSDsinkhist(V1ICenc,:), 'Color', [0 0.7 0])
set(gca, 'XTick', ctxcsdelectinds, 'XTickLabel', lfpelecvec.location(ctxcsdelectinds))
title('sink')

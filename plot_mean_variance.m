load('G:\My Drive\RESEARCH\ICexpts_revision23\openscope_psthavgall.mat')

% Siegle et al's single unit filter criteria
% isi_violations < 0.5 & amplitude_cutoff < 0.1 & presence_ratio > 0.9
neuRSall = unit_wfdur_all>0.4;
neufiltall = (unit_isi_violations_all<0.5 & unit_amplitude_cutoff_all<0.1 & unit_presence_ratio_all>0.9);


neuV1RSfiltind = find( contains(neulocall, 'VISp') & ~contains(neulocall, 'VISpm') & neuRSall & neufiltall );

blocks2plt = {'ICkcfg0_presentations', 'ICkcfg1_presentations', 'ICwcfg0_presentations', 'ICwcfg1_presentations'};
% blockcols = lines(numel(blocks2plt));
blockcols = [0 0 1; 0 0.7 0; 1 0.7 0; 1 0 0];
neu2plt = neuV1RSfiltind(randperm(numel(neuV1RSfiltind), 40));
figure
for ii = 1:numel(neu2plt)
    ci = neu2plt(ii);
    subplot(5,8,ii)
    hold all
    for b= 4 %1:numel(blocks2plt)
        highreptt = vistrialrepagg(1).(blocks2plt{b})==max(vistrialrepagg(1).(blocks2plt{b}));
        xvec = log2(0.4*Ronavgall.(blocks2plt{b})(ci,highreptt));
        yvec = log2( (0.4*Ronstdall.(blocks2plt{b})(ci,highreptt)).^2 );
        scatter( xvec, yvec, 'o', 'MarkerEdgeColor', blockcols(b,:) )
        scatter( xvec(1), yvec(1), '*', 'MarkerEdgeColor', blockcols(b,:) )
        % disp([min(xvec) max(xvec) yvec(1)-xvec(1)] )
        xl = [min(xvec(~isinf(xvec))) max(xvec(~isinf(xvec)))];
        plot(xl, xl+(yvec(1)-xvec(1)), '-', 'Color', blockcols(b,:) )
        % xmaxr = max(Ronavgall.(blocks2plt{b})(ci,:))/Ronavgall.(blocks2plt{b})(ci,1);
        % plot(xmaxr*[0 Ronavgall.(blocks2plt{b})(ci,1)], xmaxr*[0 Ronstdall.(blocks2plt{b})(ci,1).^2], '-', 'Color', blockcols(b,:))
    end
    xl = xlim;
    plot(xl,xl,'k-')
    title(neulocall{ci})
end

b=3;
[sv,si] = sort(mean(Ronavgall.(blocks2plt{b})(neuV1RSfiltind,:),2), 'descend');
neu2plt = neuV1RSfiltind(si(1:40));
figure
for ii = 1:numel(neu2plt)
    ci = neu2plt(ii);
    subplot(5,8,ii)
    hold all
    scatter(0.4*Ronavgall.(blocks2plt{b})(ci,:), (0.4*Ronstdall.(blocks2plt{b})(ci,:)).^2, 'o', 'MarkerEdgeColor', blockcols(b,:) )
    scatter(0.4*Ronavgall.(blocks2plt{b})(ci,1), (0.4*Ronstdall.(blocks2plt{b})(ci,1)).^2, '*', 'MarkerEdgeColor', blockcols(b,:) )
    xmaxr = max(Ronavgall.(blocks2plt{b})(ci,:))/Ronavgall.(blocks2plt{b})(ci,1);
    plot(xmaxr*[0 0.4*Ronavgall.(blocks2plt{b})(ci,1)], xmaxr*[0 (0.4*Ronstdall.(blocks2plt{b})(ci,1)).^2], '-', 'Color', blockcols(b,:))
    xl = xlim;
    plot(xl,xl,'k-')
    title(neulocall{ci})
end


b=1;
spkcnt = 0.4*Ronavgall.(blocks2plt{b});
FF = ( (0.4*Ronstdall.(blocks2plt{b})).^2)./(0.4*Ronavgall.(blocks2plt{b}));
FFV1 = FF(neuV1RSfiltind,:);
FFV1( any(isnan(FFV1),2), :)=[];
[p,tbl,stats]=friedman(FFV1);
figure;multcompare(stats)

spkcntevoked = spkcnt-spkcnt(:,1);

% blank trial firing rate vs fano factor
figure; hold all
plot(0.4*Ronavgall.(blocks2plt{b})(neuV1RSfiltind,1), FF(neuV1RSfiltind,1),'o')
corr(0.4*Ronavgall.(blocks2plt{b})(neuV1RSfiltind,1), FF(neuV1RSfiltind,1), 'type', 'spearman', 'rows','complete')

typi = 2;
figure; hold all
plot(spkcntevoked(neuV1RSfiltind,typi),FF(neuV1RSfiltind,typi)-FF(neuV1RSfiltind,1),'o')
xl =xlim; yl = ylim;
plot(xl,[0 0],'k-')
plot([0 0],yl,'k-')
xlabel('evoked activity')
ylabel('evoked-blank Fano Factor')

typi = 2;
figure; hold all
histogram2(spkcntevoked(neuV1RSfiltind,typi),FF(neuV1RSfiltind,typi)-FF(neuV1RSfiltind,1),'displaystyle', 'tile')
xl =xlim; yl = ylim;
plot(xl,[0 0],'w-', 'linewidth',1)
plot([0 0],yl,'w-', 'linewidth',1)
xlabel('evoked activity')
ylabel('evoked-blank Fano Factor')

% [sv,si]=sort(spkcnt,2);
[mv,mi]=max(spkcnt,[],2);
matind = sub2ind( size(spkcnt), (1:size(spkcnt,1))',mi );
if ~isequal(mv, spkcnt(matind))
    error('check matind')
end
FFpref = FF(matind);
spkcntpref = spkcnt(matind);
spkcntprefevoked = spkcntevoked(matind);
if ~isequal(spkcntprefevoked, spkcntpref-spkcnt(:,1))
    error('check')
end


figure; hold all
plot(FF(neuV1RSfiltind,1), FFpref(neuV1RSfiltind), 'o')
xl = xlim;
plot(xl,xl,'k-')
xlabel('blank FF')
ylabel('pref stim FF')
signrank(FF(neuV1RSfiltind,1), FFpref(neuV1RSfiltind)) % signifcant increase on pref stim
signrank(FF(neuV1RSfiltind,1), FFpref(neuV1RSfiltind), 'tail', 'left') 

[rstim, istim]=max(mean(Ronavgall.(blocks2plt{b})(neuV1RSfiltind,:),1));
signrank(FF(neuV1RSfiltind,1), FF(neuV1RSfiltind,istim)) % signifcantly higher on blank trials
signrank(FF(neuV1RSfiltind,1), FF(neuV1RSfiltind,istim), 'tail', 'right')

figure; hold all
histogram(FFpref(neuV1RSfiltind)-FF(neuV1RSfiltind,1), 'binwidth', 0.5) 
histogram(FF(neuV1RSfiltind,istim)-FF(neuV1RSfiltind,1), 'binwidth', 0.5) 

figure; hold all
xvec = spkcntprefevoked(neuV1RSfiltind);
yvec = FFpref(neuV1RSfiltind)-FF(neuV1RSfiltind,1);
h = histogram2(xvec, yvec, 'displaystyle', 'tile');
xbincenter = ( h.XBinEdges(1:end-1)+h.XBinEdges(2:end) )/2;
xbinmean = NaN(1,length(xbincenter));
xbinsem = NaN(1,length(xbincenter));
xbinmedian = NaN(1,length(xbincenter));
for ix = 1:length(xbincenter)
    tempxoi = xvec>=h.XBinEdges(ix) & xvec<h.XBinEdges(ix+1);
    xbinmean(ix) = nanmean(yvec(tempxoi));
    xbinsem(ix) = nanstd(yvec(tempxoi))/sqrt(nnz(~isnan(yvec(tempxoi))));
    xbinmedian(ix) = nanmedian(yvec(tempxoi));
end
errorbar(xbincenter, xbinmean, xbinsem, 'ro-', 'linewidth', 2)
plot(xbincenter, xbinmedian, 'g-', 'linewidth', 2)
length(h.XBinEdges)-1;
xl =xlim; yl = ylim;
plot(xl,[0 0],'w-', 'linewidth',1)
plot([0 0],yl,'w-', 'linewidth',1)
xlabel('evoked activity')
ylabel('pref-blank Fano Factor')

figure; hold all
xvec = Ronavgall.(blocks2plt{b})(neuV1RSfiltind,1);
yvec = FFpref(neuV1RSfiltind)-FF(neuV1RSfiltind,1);
h = histogram2(xvec, yvec, 'displaystyle', 'tile');
xbincenter = ( h.XBinEdges(1:end-1)+h.XBinEdges(2:end) )/2;
xbinmean = NaN(1,length(xbincenter));
xbinsem = NaN(1,length(xbincenter));
xbinmedian = NaN(1,length(xbincenter));
for ix = 1:length(xbincenter)
    tempxoi = xvec>=h.XBinEdges(ix) & xvec<h.XBinEdges(ix+1);
    xbinmean(ix) = nanmean(yvec(tempxoi));
    xbinsem(ix) = nanstd(yvec(tempxoi))/sqrt(nnz(~isnan(yvec(tempxoi))));
    xbinmedian(ix) = nanmedian(yvec(tempxoi));
end
errorbar(xbincenter, xbinmean, xbinsem, 'ro-', 'linewidth', 2)
plot(xbincenter, xbinmedian, 'g-', 'linewidth', 2)
length(h.XBinEdges)-1;
xl =xlim; yl = ylim;
plot(xl,[0 0],'w-', 'linewidth',1)
plot([0 0],yl,'w-', 'linewidth',1)
xlabel('blank response')
ylabel('pref-blank Fano Factor')

figure
hold all
errorbar(xbincenter, xbinmean, xbinsem, 'ro-', 'linewidth', 2)
plot(xbincenter, xbinmedian, 'g-', 'linewidth', 2)
xl = xlim;
plot(xl,[0 0],'k-', 'linewidth',1)
xlabel('blank response')
ylabel('pref-blank Fano Factor')

figure; hold all
xvec = Ronavgall.(blocks2plt{b})(neuV1RSfiltind,2)-Ronavgall.(blocks2plt{b})(neuV1RSfiltind,1);
yvec = FF(neuV1RSfiltind,2)-FF(neuV1RSfiltind,1);
h = histogram2(xvec, yvec, 'displaystyle', 'tile');
xbincenter = ( h.XBinEdges(1:end-1)+h.XBinEdges(2:end) )/2;
xbinmean = NaN(1,length(xbincenter));
xbinsem = NaN(1,length(xbincenter));
xbinmedian = NaN(1,length(xbincenter));
for ix = 1:length(xbincenter)
    tempxoi = xvec>=h.XBinEdges(ix) & xvec<h.XBinEdges(ix+1);
    xbinmean(ix) = nanmean(yvec(tempxoi));
    xbinsem(ix) = nanstd(yvec(tempxoi))/sqrt(nnz(~isnan(yvec(tempxoi))));
    xbinmedian(ix) = nanmedian(yvec(tempxoi));
end
errorbar(xbincenter, xbinmean, xbinsem, 'ro-', 'linewidth', 2)
plot(xbincenter, xbinmedian, 'g-', 'linewidth', 2)
length(h.XBinEdges)-1;
xl =xlim; yl = ylim;
plot(xl,[0 0],'w-', 'linewidth',1)
plot([0 0],yl,'w-', 'linewidth',1)
xlabel('evoked response')
ylabel('stim-blank Fano Factor')

figure
hold all
errorbar(xbincenter, xbinmean, xbinsem, 'ro-', 'linewidth', 2)
plot(xbincenter, xbinmedian, 'g-', 'linewidth', 2)
xl = xlim;
plot(xl,[0 0],'k-', 'linewidth',1)
xlabel('evoked response')
ylabel('stim-blank Fano Factor')

figure
plot(Ronavgall.(blocks2plt{b})(neuV1RSfiltind,:), FF(neuV1RSfiltind,:), 'o')

figure
histogram2( reshape(spkcnt(neuV1RSfiltind,:),[],1), reshape(FF(neuV1RSfiltind,:),[],1), 'displaystyle', 'tile')

figure
plot(Ronavgall.(blocks2plt{b})(neuV1RSfiltind,:), Ronstdall.(blocks2plt{b})(neuV1RSfiltind,:).^2, 'o')

figure
histogram2( reshape(Ronavgall.(blocks2plt{b})(neuV1RSfiltind,:),[],1), reshape(Ronstdall.(blocks2plt{b})(neuV1RSfiltind,:),[],1), 'displaystyle', 'tile')

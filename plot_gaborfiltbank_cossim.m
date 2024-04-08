phase36 = load('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130\visICtxiwcfg1\ICwcfg1_gaborfiltbank_cossim.mat');
phase2 = load('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130\visICtxiwcfg1\ICwcfg1_gaborfiltbank_cossim2.mat');
% phase1 = load('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130\visICtxiwcfg1\ICwcfg1_gaborfiltbank_cossim1.mat');


tt = {"0106", "0107", "0110", "0111", "1105", "1109"};
ttcol = [0 .4 0; .5 0.25 0; 1 0.5 0; 0 1 0; 0 1 1; 0 1 1];
figure; 
for ii = 1:2
    switch ii
        case 1
            tempphase = phase36;
        case 2
            tempphase = phase2;
    end
cossimmat = tempphase.cossimmat;
imlabel = tempphase.imlabel;

    subplot(1,2,ii)
imagesc(cossimmat); set(gca, 'XTick', 1:size(cossimmat,2), 'XTickLabel', imlabel, 'YTick', 1:size(cossimmat,1), 'YTickLabel', imlabel)
hold on
for typi = 1:numel(tt)
indim = find(strcmp(string(imlabel), tt{typi}));
plot([indim indim], [0.5 size(cossimmat,2)+0.5], '-', 'Color', ttcol(typi,:))
plot([0.5 size(cossimmat,1)+0.5], [indim indim], '-', 'Color', ttcol(typi,:))
end
colorbar
axis square
end

figure; plot(phase36.cossimmat, phase2.cossimmat, '.')
xl=xlim; hold on; plot(xl,xl,'k-')
axis square
xlabel('36 phases')
ylabel('2 phases')

corr(phase36.cossimmat(:), phase2.cossimmat(:), 'type', 'spearman')

% SHOULD GAZE BE DEFINED BASED ON EYETRACKING POSITION OR PUPILTRACKING POSITION?
% NOTE CURRENTLY I DEFINED GAZE BASED ON PUPIL POSITION
% Siegle et al Neuropixels platform paper:
% Across 50 mice with processed eye-tracking videos, we used 
% the gaze_mapping module of the AllenSDK to translate pupil position into 
% screen coordinates (in units of degrees). On average, 95% of gaze locations 
% fell within 6.4 ± 2.1° of the mean, with a maximum of 13.6°.


% <4 vis deg (stricter criterion) for fixed gaze and replicate Fig1 results (R2C1.1)
% 
% eye position on different trial types (esp. IC vs LC) (R1C1)
% pupil area on IC vs LC vs RE trials (R1C2)
% 
% Perhaps show that receptive field position is not different when using all trials vs fixed-gaze trials
% Alternatively, show that exclusively center responsive neurons defined with all trials do not respond to grating patches in other RF positions …

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );
Nsessions = numel(nwbsessions);
fs=12;
%% eye position on IC vs LC vs IRE trials
% 20 pix = 8 visual degrees
pixperdeg = 20/8;
eyecamframerate = 60;
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
whichblock = 'ICwcfg1_presentations';
tempsplt = strsplit(whichblock,'_');
whichICblock = tempsplt{1};

ises=2; pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'trialpupil.mat'], 'trackeyetli')

validfgsessions = true(Nsessions,1);
ttpupilposx = cell(numel(ICtrialtypes),Nsessions);
ttpupilposy = cell(numel(ICtrialtypes),Nsessions);
ttpupilareaz = cell(numel(ICtrialtypes),Nsessions);
trialpupilareazavg = NaN(length(trackeyetli), numel(ICtrialtypes), Nsessions);
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    if ~exist([pathpp 'trackmouseeye.mat'], 'file')
        validfgsessions(ises) = false;
        continue
    end
    load([pathpp 'postprocessed.mat'], 'vis')
    load([pathpp 'trackmouseeye.mat'], 'modecom', 'trialdistmodecom', 'pupiltracking')
    load([pathpp 'trialpupil.mat']) %'eyecamframerate', 'trackeyetli', 'trialpupildata', 'trialpupilarea'

    % trackeyetli = trialdistmodecom.(whichblock).trackeyetli;
    tloi = trackeyetli>=0 & trackeyetli<0.4*eyecamframerate;

    pupareamean = nanmean(pupiltracking.area);
    pupareastd = nanstd(pupiltracking.area);

    for t = 1:numel(ICtrialtypes)
        trialsoi = vis.(whichblock).trialorder==t-1;
        tempx = (trialpupildata.(whichblock).x(trialsoi,tloi) - modecom(1)) / pixperdeg;
        tempy = (trialpupildata.(whichblock).y(trialsoi,tloi) - modecom(2)) / pixperdeg;

        ttpupilposx{t,ises} = tempx;
        ttpupilposy{t,ises} = tempy;

        ttpupilareaz{t,ises} = (trialpupilarea.(whichblock)(trialsoi,tloi)-pupareamean)/pupareastd;
        trialpupilareazavg(:,t,ises) = nanmean(trialpupilarea.(whichblock)(trialsoi,:),1);
    end
end
valfgses = find(validfgsessions);

save(['G:\My Drive\RESEARCH\ICexpts_revision23\openscope_pupilarea_', whichICblock, '.mat'], ...
    'pixperdeg', 'eyecamframerate', 'ICtrialtypes', 'valfgses', 'trackeyetli', 'tloi', ...
    'ttpupilposx', 'ttpupilposy', 'ttpupilareaz', 'trialpupilareazavg')

%% 'ttpupilposx', 'ttpupilposy', 'ttpupilareaz', 'trialpupilareazavg' all have visual degrees as units
drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/RESEARCH/ICexpts_revision23/';
codepath = '/Users/hyeyoung/Documents/CODE/';
addpath(genpath([codepath 'helperfunctions']))


load([drivepath 'openscope_pupilarea_ICwcfg1.mat'])
load([drivepath 'openscope_popavg_all.mat'], 'nwbsessions')
Nsessions = numel(nwbsessions);
blankwcfg0 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg0/100000.tif']));
blankwcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110000.tif']));
IC1wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110106.tif']));
LC1wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110107.tif']));
LC2wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110110.tif']));
IC2wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110111.tif']));
IRE1wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110506.tif']));
IRE2wcfg1 = double(imread([codepath 'Display_IC/visOpenScope/visICtxiwcfg1/110511.tif']));

halfhorz = blankwcfg0(round(size(blankwcfg0,1)/2),:);
cumhalfhorz = cumsum(halfhorz);
npix16deg = nnz(cumhalfhorz==cumhalfhorz(round(length(cumhalfhorz)/2)))+1;
nICimpixperdeg = npix16deg/16;

fs=18;

%% 1d histogram of gaze distribution for each trial type
tt2p = [106 107 110 111];
ttdesc = {'I_C_1', 'L_C_1', 'L_C_2', 'I_C_2'};
hbe = 0:0.5:16;
hbc = ( hbe(1:end-1)+hbe(2:end) )/2;
histtrialeachses = cell(size(tt2p));
histtrialallses = cell(size(tt2p));
histtrialmaxeachses = cell(size(tt2p));
histtrialmaxallses = cell(size(tt2p));
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    
    histtrialeachses{s} = NaN(length(hbe)-1, Nsessions);
    histtrialmaxeachses{s} = NaN(length(hbe)-1, Nsessions);
    for ises = 1:Nsessions
        if isempty(ttpupilposx{1,ises})
            continue
        end
        tempx = ttpupilposx{typi,ises};
        tempy = ttpupilposy{typi,ises};
        tempdeg = sqrt(tempx.^2+tempy.^2);
        [P,EDGES]=histcounts(tempdeg(:), hbe, 'normalization', 'cdf');
        histtrialeachses{s}(:,ises) = P;
        
        [P,EDGES]=histcounts(max(tempdeg,[],2), hbe, 'normalization', 'cdf');
        histtrialmaxeachses{s}(:,ises) = P;
    end
    
    tempx = cat(1,ttpupilposx{typi,:});
    tempy = cat(1,ttpupilposy{typi,:});
    tempdeg = sqrt(tempx.^2+tempy.^2);
    [P,EDGES]=histcounts(tempdeg(:), hbe, 'normalization', 'cdf');
    histtrialallses{s} = P;

    [P,EDGES]=histcounts(max(tempdeg,[],2), hbe, 'normalization', 'cdf');
    histtrialmaxallses{s} = P;
end


tempx = cat(1,ttpupilposx{ismember(ICtrialtypes, tt2p),:});
tempy = cat(1,ttpupilposy{ismember(ICtrialtypes, tt2p),:});
tempdeg = sqrt(tempx.^2+tempy.^2);
mean(tempdeg(:)<8)

%figure('Position', [400 100 630 120])
figure('Position', [400 100 1260 240])
%annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    subplot(1,4,s)
    hold all
    plot(hbe(2:end), histtrialeachses{s})
    plot(hbe(2:end), histtrialallses{s}, 'k-', 'LineWidth', 2)
    plot([8 8], [0 1], 'k--', 'LineWidth', 0.5)
    axis([0 16 0 1])
    set(gca, 'XTick', 0:8:16, 'FontSize', fs)
    xlabel('Visual Degrees', 'FontSize', fs)
    ylabel('Gaze Deviation CDF', 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
end

%figure('Position', [400 100 630 120])
figure('Position', [400 100 1260 240])
%annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    subplot(1,4,s)
    hold all
    plot(hbe(2:end), histtrialmaxeachses{s})
    plot(hbe(2:end), histtrialmaxallses{s}, 'k-', 'LineWidth', 2)
    plot([8 8], [0 1], 'k--', 'LineWidth', 0.5)
    axis([0 16 0 1])
    set(gca, 'XTick', 0:8:16, 'FontSize', fs)
    xlabel('Visual Degrees', 'FontSize', fs)
    ylabel('Gaze Max Deviation CDF', 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
end

%% 2d histogram each session each trialtype
tt2p = [106 107 110 111 506 511];
ttdesc = {'I_C_1', 'L_C_1', 'L_C_2', 'I_C_2', 'I_R_E_1', 'I_R_E_2'};
figure
sescnt=0;
for ises = 1:Nsessions
    if isempty(ttpupilposx{1,ises})
        continue
    else
    sescnt = sescnt+1;
    end
    for s = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(s);
        tempx = ttpupilposx{typi,ises};
        tempy = ttpupilposy{typi,ises};

        subplot(numel(tt2p),10,10*(s-1)+sescnt)
        histogram2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'probability', 'displaystyle', 'tile')
        axis equal
        axis([-10 10 -10 10])
        caxis([0 0.08])
        colormap redblue

        % xl = prctile(trialpupildata.(whichblock).x(:),[0.5 99.5]);
        % yl = prctile(trialpupildata.(whichblock).y(:),[0.5 99.5]);
        % [N,XEDGES,YEDGES] = histcounts2(tempx(:),tempy(:));
        % xctrs = ( XEDGES(1:end-1)+XEDGES(2:end) )/2;
        % yctrs = ( YEDGES(1:end-1)+YEDGES(2:end) )/2;
        % hold on
        % [M,c] = contour(repmat(xctrs,length(yctrs),1), repmat(yctrs',1,length(xctrs)), N', 'linewidth', 3);
        %
        % disp(c.LevelList)
        % if ~ismember(200, c.LevelList)
        %     error('200 is not one of the levels..')
        % end
        % Mcell = cell(size(c.LevelList));
        % cnt = 1;
        % for ilev = 1:length(c.LevelList)
        %     Mcell{ilev} = M(:,cnt+1:cnt+M(2,cnt));
        %     cnt = cnt+M(2,cnt)+1;
        % end
        % plot(Mcell{c.LevelList==250}(1,:), Mcell{c.LevelList==250}(2,:), 'r:', 'linewidth', 2)
        % axis([xl yl])

        title(sprintf('%d %s %d', ises, nwbsessions{ises}, tt2p(s)))
        colorbar
    end
end

%% contour plot of gaze position across sessions
contourlevel = 0.01;
contourl = cell(numel(tt2p),Nsessions);
contourlcell = cell(numel(tt2p),Nsessions);

figure
for s = 1:numel(tt2p)
    sescnt=0;
    subplot(2,3,s)
    hold all
    for ises = 1:Nsessions
        typi = ICtrialtypes==tt2p(s);
        if isempty(ttpupilposx{typi,ises})
            continue
        else
            sescnt = sescnt+1;
        end
        % subplot(numel(tt2p),10,10*(s-1)+sescnt)

        tempx = ttpupilposx{typi,ises};
        tempy = ttpupilposy{typi,ises};

        % histogram2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'pdf', 'displaystyle', 'tile')
        % axis equal
        % axis([-10 10 -10 10])
        % % clim([0 0.08])
        % colorbar
        % colormap redblue

        % xl = prctile(trialpupildata.(whichblock).x(:),[0.5 99.5]);
        % yl = prctile(trialpupildata.(whichblock).y(:),[0.5 99.5]);
        [N,XEDGES,YEDGES] = histcounts2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'pdf');
        xctrs = ( XEDGES(1:end-1)+XEDGES(2:end) )/2;
        yctrs = ( YEDGES(1:end-1)+YEDGES(2:end) )/2;
        hold on
        [M,c] = contour(repmat(xctrs,length(yctrs),1), repmat(yctrs',1,length(xctrs)), N', [contourlevel contourlevel], 'linewidth', 1);
        contourl{s,ises} = M;

        axis([-10 10 -10 10])

        % disp(c.LevelList)
        % if ~ismember(0.04, c.LevelList)
        %     error('200 is not one of the levels..')
        % end
        Mcell = cell(size(c.LevelList));
        cnt = 1;
        ilev = 1;
        while cnt<size(M,2)
            Mcell{ilev} = M(:,cnt+1:cnt+M(2,cnt));
            cnt = cnt+M(2,cnt)+1;
            ilev = ilev+1;
        end
        contourlcell{s,ises} = Mcell;
        % plot(Mcell{c.LevelList==0.04}(1,:), Mcell{c.LevelList==0.04}(2,:), 'r:', 'linewidth', 2)
        % axis([xl yl])

        title(sprintf('%d %s %d', ises, nwbsessions{ises}, tt2p(s)))
    end
end

% %% across sessions 1 percentile contour plot
linecols = lines(Nsessions);
figure('Position', [400 100 630 180])
annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    sescnt=0;
    subplot(1,4,s)
    hold all
    for ises = 1:Nsessions
        if ~isempty(contourlcell{s,ises})
            for ii = 1:numel(contourlcell{s,ises})
                if ~isempty(contourlcell{s,ises}{ii})
                plot(contourlcell{s,ises}{ii}(1,:), contourlcell{s,ises}{ii}(2,:), 'Color', linecols(ises,:))
                end
            end
        end
    end
    axis([-10 10 -10 10])
    axis square
    set(gca, 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
    % xlabel('Pupil Position (vis. deg.)')
    % if s==1
    %     ylabel('Pupil Position (vis. deg.)', 'FontSize', fs)
    % end
end

linecols = lines(Nsessions);
%figure('Position', [400 100 630 120])
figure('Position', [400 100 1260 240])
%annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    switch tt2p(s)
        case 106
            tempim = IC1wcfg1;
        case 107
            tempim = LC1wcfg1;
        case 110
            tempim = LC2wcfg1;
        case 111
            tempim = IC2wcfg1;
        case 506
            tempim = IRE1wcfg1;
        case 511
            tempim = IRE2wcfg1;
        otherwise
            error('designate tempim')
    end
    sescnt=0;
    subplot(1,4,s)
    xax = -floor(size(blankwcfg1,2)/2)+(0:size(blankwcfg1,2)-1);
    xax = xax/nICimpixperdeg;
    yax = -floor(size(blankwcfg1,1)/2)+(0:size(blankwcfg1,1)-1);
    yax = yax/nICimpixperdeg;
    imagesc(xax, yax, repmat(tempim,1,1,3), 'AlphaData', 0.1)
    hold on
    for ises = 1:Nsessions
        if ~isempty(contourlcell{s,ises})
            for ii = 1:numel(contourlcell{s,ises})
                if ~isempty(contourlcell{s,ises}{ii})
                plot(contourlcell{s,ises}{ii}(1,:), contourlcell{s,ises}{ii}(2,:), 'Color', linecols(ises,:), 'linewidth', 2)
                end
            end
        end
    end
    axis equal
    axis([xax(1) xax(end) yax(1) yax(end)])
    axis off
    %set(gca, 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
    % xlabel('Pupil Position (vis. deg.)')
    % if s==1
    %     ylabel('Pupil Position (vis. deg.)', 'FontSize', fs)
    % end
end

%% example session
ises = 2;
pltcb=false; 
cl=[0 0.08];
figure('Position', [400 100 630 180])
annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    tempx = ttpupilposx{typi,ises};
    tempy = ttpupilposy{typi,ises};

    subplot(1,4,s)
    h = histogram2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'pdf', 'displaystyle', 'tile', 'edgecolor', 'none');
    caxis(cl)
    set(gca, 'FontSize', fs)
    axis([-10 10 -10 10])
    if pltcb && s==4
        set(gca,'XTick',[],'YTick',[])
    cb = colorbar;
    cb.Ticks = [cl];
    cb.FontSize = fs;
    else
    title(ttdesc{s}, 'FontSize', fs)
    end
    axis square
    colormap redblue
    % xlabel('Pupil Position (vis. deg.)')
    % if s==1
    %     ylabel('Pupil Position (vis. deg.)', 'FontSize', fs)
    % end
end


%figure('Position', [400 100 630 120])
figure('Position', [400 100 1260 240])
%annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    tempx = ttpupilposx{typi,ises};
    tempy = ttpupilposy{typi,ises};
    
    switch tt2p(s)
        case 106
            tempim = IC1wcfg1;
        case 107
            tempim = LC1wcfg1;
        case 110
            tempim = LC2wcfg1;
        case 111
            tempim = IC2wcfg1;
        case 506
            tempim = IRE1wcfg1;
        case 511
            tempim = IRE2wcfg1;
        otherwise
            error('designate tempim')
    end
    sescnt=0;
    subplot(1,4,s)
    xax = -floor(size(blankwcfg1,2)/2)+(0:size(blankwcfg1,2)-1);
    xax = xax/nICimpixperdeg;
    yax = -floor(size(blankwcfg1,1)/2)+(0:size(blankwcfg1,1)-1);
    yax = yax/nICimpixperdeg;
    imagesc(xax, yax, repmat(tempim,1,1,3), 'AlphaData', 0.1)
    hold on
    h = histogram2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'pdf', 'displaystyle', 'tile', 'edgecolor', 'none');
    caxis(cl)
    if pltcb && s==4
    cb = colorbar;
    cb.Ticks = [cl];
    cb.FontSize = fs;
    else
    title(ttdesc{s}, 'FontSize', fs)
    end
    axis equal
    axis([xax(1) xax(end) yax(1) yax(end)])
    axis off
    %set(gca, 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
    % xlabel('Pupil Position (vis. deg.)')
    % if s==1
    %     ylabel('Pupil Position (vis. deg.)', 'FontSize', fs)
    % end
end
colormap redblue

%%
% density plot with transparent scatterplot
%figure('Position', [400 100 630 120])
figure('Position', [400 100 1260 240])
%annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    tempx = ttpupilposx{typi,ises};
    tempy = ttpupilposy{typi,ises};
    
    switch tt2p(s)
        case 106
            tempim = IC1wcfg1;
        case 107
            tempim = LC1wcfg1;
        case 110
            tempim = LC2wcfg1;
        case 111
            tempim = IC2wcfg1;
        case 506
            tempim = IRE1wcfg1;
        case 511
            tempim = IRE2wcfg1;
        otherwise
            error('designate tempim')
    end
    sescnt=0;
    subplot(1,4,s)
    %figure
    xax = -floor(size(blankwcfg1,2)/2)+(0:size(blankwcfg1,2)-1);
    xax = xax/nICimpixperdeg;
    yax = -floor(size(blankwcfg1,1)/2)+(0:size(blankwcfg1,1)-1);
    yax = yax/nICimpixperdeg;
    imagesc(xax, yax, repmat(tempim,1,1,3), 'AlphaData', 0.1)
    hold on
    scatter(tempx(:),tempy(:), 1, 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.1);
    caxis(cl)
    if pltcb && s==4
    cb = colorbar;
    cb.Ticks = [cl];
    cb.FontSize = fs;
    else
    title(ttdesc{s}, 'FontSize', fs)
    end
    axis equal
    axis([xax(1) xax(end) yax(1) yax(end)])
    axis off
    %set(gca, 'FontSize', fs)
    title(ttdesc{s}, 'FontSize', fs)
    % xlabel('Pupil Position (vis. deg.)')
    % if s==1
    %     ylabel('Pupil Position (vis. deg.)', 'FontSize', fs)
    % end
end

%% pupil area on IC vs LC vs RE trials (R1C2)
trackeyetl = trackeyetli/eyecamframerate;
tt2p = [106 107 110 111 506 511];
ttdesc = {'I_C_1', 'L_C_1', 'L_C_2', 'I_C_2', 'I_R_E_1', 'I_R_E_2'};
ttcol = [0 .4 0; .5 0.25 0; 1 0.5 0; 0 1 0; 0 0 0.4; 0 0 1];

tt2pinds = zeros(size(tt2p));
for s = 1:numel(tt2p)
    tt2pinds(s) = find(ICtrialtypes==tt2p(s));
end

figure; hold all
for s = 1:numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    temparea = squeeze(nanmean(cat(3,ttpupilareaz{typi,:}),[1,2]));
    histogram(temparea(:), -3:0.1:3, 'normalization', 'probability')
end

% figure('Position', [800 300 240 200])
figure('Position', [800 300 300 240])
hold all
for s = 1:numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    temparea = squeeze(nanmean(cat(3,ttpupilareaz{typi,:}),[1,2]));
b = boxchart(s*ones(numel(temparea),1), temparea(:), 'BoxFaceColor', ttcol(s,:), 'MarkerStyle', 'none', 'linewidth', 1);%, 'Notch' , 'on'); %, 'FontName', 'Arial')
%b.MarkerColor = ttcol(s,:);
b.BoxFaceColor = ttcol(s,:);
b.BoxEdgeColor = ttcol(s,:);
b.BoxFaceAlpha = 0.5;
end
xlim([0.5 numel(tt2p)+0.5])
ylim([-1 1])
set(gca, 'FontSize', fs, 'XTick', 1:numel(tt2p), 'XTickLabel', ttdesc, 'XTickLabelRotation', 0)
ylabel('z-Pupil Area', 'FontSize', fs)
title(' ', 'FontSize', fs)
xlabel(' ', 'FontSize', fs)

figure
hold all
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
plot(trackeyetl, squeeze(trialpupilareazavg(:,typi,:)), 'Color', ttcol(s,:))
end
xlabel('Time (s)')
ylabel('z-Pupil Area')
title('Each Session')

figure
hold all
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
plot(trackeyetl, nanmean(trialpupilareazavg(:,typi,:),3), 'Color', ttcol(s,:), 'LineWidth', 1)
end
xlabel('Time (s)')
ylabel('z-Pupil Area')
title('Average across Sessions')


tempmat = squeeze(mean(trialpupilareazavg(tloi, tt2pinds,:),1))';
tempmat(any(isnan(tempmat),2), :)=[];
[p,tbl,stats] = friedman(tempmat);
figure; multcompare(stats)
disp(p)

ttpupilareazavg = NaN(numel(valfgses), numel(tt2p));
for s = 1:numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    ttpupilareazavg(:,s) = squeeze(nanmean(cat(3,ttpupilareaz{typi,:}),[1,2]));
end
[p,tbl,stats] = friedman(ttpupilareazavg);
figure; multcompare(stats)
disp(p) % 0.7220

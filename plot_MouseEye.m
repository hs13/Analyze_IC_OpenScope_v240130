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
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

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

ises=2; pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'trialpupil.mat'], 'trackeyetli')

ttpupilposx = cell(numel(ICtrialtypes),Nsessions);
ttpupilposy = cell(numel(ICtrialtypes),Nsessions);
ttpupilareaz = cell(numel(ICtrialtypes),Nsessions);
trialpupilareazavg = NaN(length(trackeyetli), numel(ICtrialtypes), Nsessions);
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    if ~exist([pathpp 'trackmouseeye.mat'], 'file')
        continue
    end
    load([pathpp 'postprocessed.mat'], 'vis')
    load([pathpp 'trackmouseeye.mat'], 'modecom', 'trialdistmodecom', 'pupiltracking')
    load([pathpp 'trialpupil.mat']) %'eyecamframerate', 'trackeyetli', 'trialpupildata', 'trialpupilarea'

    % trackeyetli = trialdistmodecom.(whichblock).trackeyetli;
    tloi = trackeyetli>0 & trackeyetli<=0.4*eyecamframerate;

    pupareamean = nanmean(pupiltracking.area);
    pupareastd = nanstd(pupiltracking.area);

    for t = 1:numel(ICtrialtypes)
        trialsoi = vis.(whichblock).trialorder==t-1;
        tempx = (trialpupildata.(whichblock).x(trialsoi,tloi) - modecom(1)) / pixperdeg;
        tempy = (trialpupildata.(whichblock).y(trialsoi,tloi) - modecom(2)) / pixperdeg;

        ttpupilposx{t,ises} = tempx;
        ttpupilposy{t,ises} = tempy;

        ttpupilareaz{t,ises} = (trialpupilarea.(whichblock)(trialsoi,tloi)-pupareamean)/pupareastd;
        trialpupilareazavg(:,t,ises) = mean(trialpupilarea.(whichblock)(trialsoi,:),1);
    end
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
        clim([0 0.08])
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
contour01 = cell(numel(tt2p),Nsessions);
contour01cell = cell(numel(tt2p),Nsessions);

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
        [M,c] = contour(repmat(xctrs,length(yctrs),1), repmat(yctrs',1,length(xctrs)), N', [0.01 0.01], 'linewidth', 1);
        contour01{s,ises} = M;

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
        contour01cell{s,ises} = Mcell;
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
        if ~isempty(contour01cell{s,ises})
            for ii = 1:numel(contour01cell{s,ises})
                plot(contour01cell{s,ises}{ii}(1,:), contour01cell{s,ises}{ii}(2,:), 'Color', linecols(ises,:))
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

%% example session
ises = 2;
pltcb=false; cl=[0 0.08];
figure('Position', [400 100 630 180])
annotation('textbox',[0 0.05 1 0.1], 'string', 'Pupil Position (vis. deg.)', 'FontSize', fs, 'edgecolor', 'none', 'HorizontalAlignment','center')
for s = 1:4%numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    tempx = ttpupilposx{typi,ises};
    tempy = ttpupilposy{typi,ises};

    subplot(1,4,s)
    h = histogram2(tempx(:),tempy(:), 'binwidth', 1, 'normalization', 'pdf', 'displaystyle', 'tile', 'edgecolor', 'none');
    set(gca, 'FontSize', fs)
    clim(cl)
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

%% pupil area on IC vs LC vs RE trials (R1C2)
trackeyetl = trackeyetli/eyecamframerate;
tt2p = [106 107 110 111 506 511];
ttdesc = {'I_C_1', 'L_C_1', 'L_C_2', 'I_C_2', 'I_R_E_1', 'I_R_E_2'};
ttcol = [0 .4 0; .5 0.25 0; 1 0.5 0; 0 1 0; 0 0 0.4; 0 0 1];

figure; hold all
for s = 1:numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
% plot(trackeyetl, squeeze(trialpupilareazavg(:,typi,:)) )
plot(trackeyetl, nanmean(trialpupilareazavg(:,typi,:),3), '-', 'color', ttcol(s,:), 'linewidth', 2)
end

tt2p = [106 107 110 111 506 511];
figure; hold all
for s = 1:numel(tt2p)
    typi = ICtrialtypes==tt2p(s);
    temparea = cat(3,ttpupilareaz{typi,:});
    histogram(temparea(:), -3:0.1:3, 'normalization', 'probability')
end

% figure('Position', [800 300 240 200])
figure('Position', [800 300 300 240])
hold all
for s = 1:numel(tt2p)
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

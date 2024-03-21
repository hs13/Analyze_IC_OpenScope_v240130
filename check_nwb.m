% requested by Ahad: Could you look at sub-630506_ses-1192952692_ogen.nwb just to confirm everything looks correct?
if ismac
    addpath(genpath('/Users/hyeyoung/Documents/CODE/Analyze_OpenScope/matnwb'))
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
    addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v230821')
    addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))
    drivepath = 'G:\My Drive\';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath(genpath([codepath 'Analyze_IC_OpenScope_v230821']))
addpath(genpath([codepath 'helperfunctions']))

% nwbspikefile = "G:\My Drive\RESEARCH\IllusionOpenScope\sub-625554_ses-1181330601-acq-FINAL_ogen.nwb";
% nwbspikefile = "G:\My Drive\RESEARCH\IllusionOpenScope\sub-625554_ses-1181330601.nwb";
% nwbspikefile = "G:\My Drive\RESEARCH\IllusionOpenScope\sub-625554_ses-1181330601_ogen.nwb";
nwbspikefile = "G:\My Drive\RESEARCH\IllusionOpenScope\sub-619296_ses-1187930705_ogen.nwb";
nwb = nwbRead(nwbspikefile);

%% check opto
optopsth_v230821 = load('S:\OpenScopeData\00248_v230821\postprocessed\sub-625554\psth_opto_probeC.mat');
optostim_v230821 = optopsth_v230821.opto.actualcond(optopsth_v230821.opto.optotrials)';

optocond = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('condition').data.load();
optostim = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('stimulus_name').data.load();
optodur = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('duration').data.load();
optolevel = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('level').data.load();

% "A Single
% unique(optocond)
%     {'1 second square pulse: continuously on for 1s'}
%     {'10 ms pulses at 1 Hz'                         }
%     {'2 ms pulses at 1 Hz'                          }
%     {'A single 30hz pulse'                          }
%     {'A single 40hz pulse'                          }
%     {'A single 50hz pulse'                          }
%     {'a single 10hz pulse'                          }
%     {'a single 20hz pulse'                          }
%     {'a single 5hz pulse'                           }
%     {'a single 60hz pulse'                          }
%     {'a single 80hz pulse'                          }
%     {'cosine pulse'                                 }


temp = [optostim, optostim_v230821];
sprintf('')
tmp = sortrows(temp,1);
unique(tmp, 'rows')

%% check vis
% changes in visblock keys
% in IC blocks, "frame" -> "Image"
% in RFCI and sizeCI block, viskeys "orientation" -> "Ori"
% in RFCI block, "x_position" -> "Pos_x", "y_position" -> "Pos_y"
% request adding fields: trialorder, trialtypedescription,
% MaskDiaVisDeg (RFCI_presentations and sizeCI_presentations),
% RFcentersVisDeg (RFCI_presentations)
visblocks = nwb.intervals.keys;
for b = 1:numel(visblocks)
    disp(visblocks{b})
    viskeys =  nwb.intervals.get(visblocks{b}).vectordata.keys;
    disp(viskeys)
end

visblocks = nwb.intervals.keys;
vis = struct();
for b = 1:numel(visblocks)
    disp(visblocks{b})
    vis.(visblocks{b}).start_time = nwb.intervals.get(visblocks{b}).start_time.data.load();
    vis.(visblocks{b}).stop_time = nwb.intervals.get(visblocks{b}).stop_time.data.load();
    viskeys =  nwb.intervals.get(visblocks{b}).vectordata.keys;
    for k = 1:numel(viskeys)
        vis.(visblocks{b}).(viskeys{k}) = nwb.intervals.get(visblocks{b}).vectordata.get(viskeys{k}).data.load();
    end
    
    % IC blocks
    if ismember('frame', viskeys)
        % expect 61, 31, 31, 31 for frame_firsttrial
        frame_firsttrial = find(vis.(visblocks{b}).frame~=0, 1, 'first');
        if mod(frame_firsttrial, 10) ~=1
            warning('first trial was blank')
            frame_firsttrial = 10*floor(frame_firsttrial/10)+1;
        end
        frame_lasttrial = find(vis.(visblocks{b}).frame~=0, 1, 'last');
        if mod(frame_lasttrial, 10) ~=9
            warning('last trial was blank')
            frame_lasttrial = 10*floor(frame_lasttrial/10)+9;
        end
        
        trialframeinds = frame_firsttrial:2:frame_lasttrial;
        vis.(visblocks{b}).trialtrialorderinds = trialframeinds;
        vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time(trialframeinds);
        vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time(trialframeinds);
        vis.(visblocks{b}).numtrials = length(trialframeinds);
        vis.(visblocks{b}).trialorder = vis.(visblocks{b}).frame(trialframeinds);
        if contains(visblocks{b}, 'cfg1')
            vis.(visblocks{b}).trialtypedescription = {'Blank', 'X', 'TC1', 'IC1', 'LC1', 'TC2', 'LC2', 'IC2', ...
                'IRE1', 'IRE2', 'TRE1', 'TRE2', 'XRE1', 'XRE2', ...
                'InBR', 'InBL', 'InTL', 'InTR', 'OutBR', 'OutBL', 'OutTL', 'OutTL'};
        elseif contains(visblocks{b}, 'cfg0')
            vis.(visblocks{b}).trialtypedescription = {'Blank', 'X', 'TC1', 'IC1', 'LC1', 'TC2', 'LC2', 'IC2', ...
                'IRE1', 'IRE2', 'TRE1', 'TRE2', 'XRE1', 'XRE2', ...
                'InR', 'InB', 'InL', 'InT', 'OutR', 'OutB', 'OutL', 'OutT'};
        else
            error('unrecognized configuration')
        end
        vis.(visblocks{b}).ICtrialtypes = [0 101 105 106 107 109 110 111 ...
            506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        
        disp([frame_firsttrial, frame_lasttrial vis.(visblocks{b}).numtrials])
        if contains(visblocks{b}, 'ICwcfg1')
            expectedNtrials = 5300; % 12*400+10*50
        else
            expectedNtrials = 22*30;
        end
        if ~( vis.(visblocks{b}).numtrials==expectedNtrials )
            error('check numtrials')
        end
    end
    
    % RFCIblocks
    % (+right,+up))
    % rfpos = [(0,0), (0,-203.3786), (203.3786/2**0.5,-203.3786/2**0.5), (203.3786,0), \
    %             (203.3786/2**0.5,203.3786/2**0.5), (0,203.3786), (-203.3786/2**0.5,203.3786/2**0.5), \
    %             (-203.3786,0), (-203.3786/2**0.5,-203.3786/2**0.5)]
    % '10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
    if ismember('y_position', viskeys)
        vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time;
        vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time;
        
        yx_position = [vis.(visblocks{b}).y_position vis.(visblocks{b}).x_position];
        vis.(visblocks{b}).sizepix = 203.3786;
        vis.(visblocks{b}).MaskDiaVisDeg = 16;
        % note, made the order of RFcentersrel match that in matlab (start
        % down, then move counterclockwise)
        vis.(visblocks{b}).RFcentersrel = [-1 0
            -1/sqrt(2) 1/sqrt(2)
            0 1
            1/sqrt(2) 1/sqrt(2)
            1 0
            1/sqrt(2) -1/sqrt(2)
            0 -1
            -1/sqrt(2) -1/sqrt(2)];
        vis.(visblocks{b}).RFcenters = vis.(visblocks{b}).sizepix * vis.(visblocks{b}).RFcentersrel;
        vis.(visblocks{b}).RFcentersVisDeg = vis.(visblocks{b}).MaskDiaVisDeg * vis.(visblocks{b}).RFcentersrel;
        if ~all(ismember(vis.(visblocks{b}).RFcenters, unique(yx_position, 'rows'), 'rows'))
            error('check RFcenters')
        end
        
        vis.(visblocks{b}).directions = unique(vis.(visblocks{b}).orientation);
        vis.(visblocks{b}).MaskList = unique(vis.(visblocks{b}).Mask);
        disp(vis.(visblocks{b}).MaskList)
        
        % vertical is zero, then rotates clockwise (45 is 1.5o'clock)
        vis.(visblocks{b}).numtrials = length(vis.(visblocks{b}).orientation);
        vis.(visblocks{b}).trialorder = zeros(vis.(visblocks{b}).numtrials, 1);
        for typi = 1:numel(vis.(visblocks{b}).directions)
            trialsoi = vis.(visblocks{b}).orientation==vis.(visblocks{b}).directions(typi);
            vis.(visblocks{b}).trialorder(trialsoi) = typi + vis.(visblocks{b}).trialorder(trialsoi);
        end
        for typi = 1:size(vis.(visblocks{b}).RFcenters,1)
            trialsoi = ismember(yx_position, vis.(visblocks{b}).RFcenters(typi,:), 'rows');
            vis.(visblocks{b}).trialorder(trialsoi) = 10*typi + vis.(visblocks{b}).trialorder(trialsoi);
        end
        for typi = 1:numel(vis.(visblocks{b}).MaskList)
            trialsoi = strcmp(vis.(visblocks{b}).Mask, vis.(visblocks{b}).MaskList(typi));
            tempss = strsplit(vis.(visblocks{b}).MaskList{typi}, '\');
            tempss = strsplit(tempss{end}, '.tif');
            maskno = str2num(tempss{1});
            vis.(visblocks{b}).trialorder(trialsoi) = maskno + vis.(visblocks{b}).trialorder(trialsoi);
        end
        vis.(visblocks{b}).trialtypedescription = ['10000s: classic 0 vs inverse 1,', ...
            ' 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'];
    end
    
    if contains(visblocks{b}, 'sizeCI')
        vis.(visblocks{b}).trialstart = vis.(visblocks{b}).start_time;
        vis.(visblocks{b}).trialend = vis.(visblocks{b}).stop_time;
        
        vis.(visblocks{b}).directions = unique(vis.(visblocks{b}).orientation);
        vis.(visblocks{b}).MaskList = unique(vis.(visblocks{b}).Mask);
        disp(vis.(visblocks{b}).MaskList)
        
        vis.(visblocks{b}).numtrials = length(vis.(visblocks{b}).orientation);
        vis.(visblocks{b}).trialorder = zeros(vis.(visblocks{b}).numtrials, 1);
        for typi = 1:numel(vis.(visblocks{b}).directions)
            trialsoi = vis.(visblocks{b}).orientation==vis.(visblocks{b}).directions(typi);
            vis.(visblocks{b}).trialorder(trialsoi) = typi + vis.(visblocks{b}).trialorder(trialsoi);
        end
        for typi = 1:numel(vis.(visblocks{b}).MaskList)
            trialsoi = strcmp(vis.(visblocks{b}).Mask, vis.(visblocks{b}).MaskList(typi));
            tempss = strsplit(vis.(visblocks{b}).MaskList{typi}, '\');
            tempss = strsplit(tempss{end}, '.tif');
            maskno = str2num(tempss{1});
            vis.(visblocks{b}).trialorder(trialsoi) = maskno + vis.(visblocks{b}).trialorder(trialsoi);
        end
        vis.(visblocks{b}).MaskDiaVisDeg = [0, 4, 8, 16, 32, 64];
        vis.(visblocks{b}).trialtypedescription = ['10000s: classic 0 vs inverse 1,', ...
            ' 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'];
    end
    
    % disp(unique(vis.(visblocks{b}).trialorder)')
end

%% check that stimulus_templates order is consistent with image file name order in all blocks
% This is especially a problem in the four IC blocks (ICwcfg1, ICwcfg0,
% ICkcfg1, ICkcfg0) because the "Image" (used to be "frame") field numbers
% do not match up with the stimulus_templates order...
% I really think the best way to avoid confusion would be to call them by their actual file names (110000, 110101, etc) instead of "ICwcfg1_presentations0", "ICwcfg1_presentations1", etc.

Nblocks2plot = numel(visblocks)-1;
figure
for b = 1:Nblocks2plot;%numel(visblocks)-1
    disp(visblocks{b})
    if contains(visblocks{b}, 'IC')
        trialtypes = vis.(visblocks{b}).ICtrialtypes;
    else
        trialtypes = unique(vis.(visblocks{b}).trialorder);
    end
    imkeys = nwb.stimulus_templates.get(visblocks{b}).image.keys;
    for ii = 1:numel(imkeys)
        %tempim=nwb.stimulus_templates.get(visblocks{b}).image.get(imkeys{ii}).data.load();
        tempimkey = [visblocks{b} num2str(ii-1)];
        tempim=nwb.stimulus_templates.get(visblocks{b}).image.get(tempimkey).data.load();
        if size(tempim,1)==3
            tempim = permute(tempim, [3 2 1]);
        else
            tempim = permute(tempim, [2 1]);
        end
        subplot(Nblocks2plot, 22, 22*(b-1)+ii)
        imshow(tempim)
        if contains(visblocks{b}, 'RFCI')
            hold on
            plot(size(tempim,1)/2*[1 1], [0.5 size(tempim,2)+0.5], 'r-')
            plot([0.5 size(tempim,1)+0.5], size(tempim,2)/2*[1 1], 'r-')
        end
        % title([num2str(trialtypes(ii)) ': ' imkeys{ii}], 'FontSize', 8)
        title(trialtypes(ii), 'FontSize', 8)
    end
end


visICblocks = {'visICtxiwcfg1', 'visICtxiwcfg0', 'visICtxikcfg1', 'visICtxikcfg0'};
blankno = [110000 100000 010000 000000];
figure
for b = 1:numel(visICblocks)
    for ii = 1:numel(ICtrialtypes)
        tempim = imread(['G:\My Drive\RESEARCH\IllusionOpenScope\IllusionOpenScope_220314\' visICblocks{b} '\' sprintf('%06d', blankno(b)+ICtrialtypes(ii)), '.tif']);
        subplot(numel(visICblocks),numel(ICtrialtypes), numel(ICtrialtypes)*(b-1)+ii)
        imshow(tempim)
        title(ICtrialtypes(ii))
    end
end

for b = 1:numel(visICblocks)
    figure
    for ii = 1:numel(ICtrialtypes)
        tempim = imread(['G:\My Drive\RESEARCH\IllusionOpenScope\IllusionOpenScope_220314\' visICblocks{b} '\' sprintf('%06d', blankno(b)+ICtrialtypes(ii)), '.tif']);
        subplot(4,6,ii)
        imshow(tempim)
        title(ICtrialtypes(ii))
    end
end


%% IC blocks figure
ICblocksimgs = zeros(1200, 1200*8);

for b = 1:4
    tempimkey = [visblocks{b} '0'];
    tempim = nwb.stimulus_templates.get(visblocks{b}).image.get(tempimkey).data.load();
end


%% RFCI and sizeCI example figure

%% all image figure
% figure out number of pixelsthat correspond to 16 degrees
tempim = nwb.stimulus_templates.get('ICwcfg0_presentations').image.get('ICwcfg0_presentations0').data.load();
tempim = tempim';
tempvec = tempim(size(tempim,1)/2,:);
tempcumvec = cumsum(double(tempvec));
visdeg16 = nnz(tempcumvec==tempcumvec(round(length(tempcumvec)/2)));

% figure; hold all
% plot(tempcumvec, 'linewidth', 1)
% plot(visdeg16*[-0.5 0.5]+round(length(tempcumvec)/2), tempcumvec(round(length(tempcumvec)/2))*[1 1], 'r-')
% figure; hold all
% plot(tempim(:,size(tempim,2)/2), 'k-', 'linewidth', 2)
% plot(visdeg16*[-0.5 0.5]+round(size(tempim,1)/2),[0 0], 'r-', 'linewidth', 1)

% images are 1920 by 1200
ICtrialtypedescription = {'Blank', 'X', 'T_C_1', 'I_C_1', 'L_C_1', 'T_C_2', 'L_C_2', 'I_C_2', ...
    'I_R_E_1', 'I_R_E_2', 'T_R_E_1', 'T_R_E_2', 'X_R_E_1', 'X_R_E_2', ...
    'In_B_R', 'In_B_L', 'In_T_L', 'In_T_R', 'Out_B_R', 'Out_B_L', 'Out_T_L', 'Out_T_L'};

% put all images into a 2 by 11 grid
whichblock = 'ICwcfg1_presentations';
ICwcfg1allimgs = zeros(1200*2, 1200*11);
imkeys = nwb.stimulus_templates.get(whichblock).image.keys;
for ii = 1:numel(imkeys)
    r = ceil(ii/11); c=mod(ii-1,11)+1;
    tempimkey = [whichblock num2str(ii-1)];
    tempim = nwb.stimulus_templates.get(whichblock).image.get(tempimkey).data.load();
    if size(tempim,1)==3
        tempim = permute(tempim, [3 2 1]);
    else
        tempim = permute(tempim, [2 1]);
    end
    tempim = tempim(:, size(tempim,2)/2-600+1:size(tempim,2)/2+600);
    ICwcfg1allimgs(1200*(r-1)+1:1200*r, 1200*(c-1)+1:1200*c) = tempim;
end
ICwcfg1allimgs(1200+[-2:3],:)=1;
for c = 1:11-1
    ICwcfg1allimgs(:,1200*c+[-2:3])=1;
end

% figure;
% imshow(ICwcfg1allimgs)
% hold on
% plot(30+[0 visdeg16], 1200-60+[0 0], 'c-', 'LineWidth', 2)
% text(30+visdeg16/2, 1200-60, '16°', 'Color', 'c', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% for ii = 1:numel(ICtrialtypedescription)
%     r = ceil(ii/11); c=mod(ii-1,11)+1;
% text(30+(c-1)*1200,30+(r-1)*1200, sprintf('(%d) %s', ii-1, ICtrialtypedescription{ii}), 'Color', 'c', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% end

% put all images into a 3 by 8 grid
rlist = [1*ones(1,8) 2*ones(1,6) 3*ones(1,8)];
clist = [1:8 1:6 1:8];
ICwcfg1_allimgs = zeros(1200*3, 1200*8);
imkeys = nwb.stimulus_templates.get(whichblock).image.keys;
for ii = 1:numel(imkeys)
    r = rlist(ii); c = clist(ii);
    tempimkey = [whichblock num2str(ii-1)];
    tempim = nwb.stimulus_templates.get(whichblock).image.get(tempimkey).data.load();
    if size(tempim,1)==3
        tempim = permute(tempim, [3 2 1]);
    else
        tempim = permute(tempim, [2 1]);
    end
    tempim = tempim(:, size(tempim,2)/2-600+1:size(tempim,2)/2+600);
    ICwcfg1_allimgs(1200*(r-1)+1:1200*r, 1200*(c-1)+1:1200*c) = tempim;
end
for r = 1:3
    if r<3
        ICwcfg1_allimgs(1200*r+[-4:5],:)=1;
    end
    for c = clist(rlist==r)
        ICwcfg1_allimgs(1200*(r-1)+1:1200*r, 1200*c+[-4:5])=1;
    end
end

fs=20;
figure;
imshow(ICwcfg1_allimgs)
hold on
plot(30+[0 visdeg16], 1200-45+[0 0], 'c-', 'LineWidth', 2)
text(30+visdeg16/2, 1200-45, '16°', 'Color', 'c', 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
for ii = 1:numel(ICtrialtypedescription)
    r = rlist(ii); c = clist(ii);
    text(30+(c-1)*1200,0+(r-1)*1200, sprintf('(%d) %s', ii-1, ICtrialtypedescription{ii}), 'Color', 'c', 'FontName', 'Arial Bold', 'FontSize', fs, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
end

%% resize RFCI and sizeCI (only the vertical orientation for now 
% need to figure out something else for other orientations
RFCIorigpath = [codepath 'Display_IC' filesep 'visRFpos' filesep];
RFCIsavepath = [drivepath 'RESEARCH/IllusionOpenScope/IllusionOpenScope_220314/visRFCI/'];

sizeCIorigpath = [codepath 'Display_IC' filesep 'vissize' filesep];
sizeCIsavepath = [drivepath 'RESEARCH/IllusionOpenScope/IllusionOpenScope_220314/vissizeCI/'];
tiflist=dir([sizeCIorigpath '*.tif']);
for ii = 1:numel(tiflist)
    tf = strsplit(tiflist(ii).name, '.');
    if mod(str2double(tf{1}),4)==1
    tempim = imread([tiflist(ii).folder filesep tiflist(ii).name]);
    newim = zeros(1200, 1920, class(tempim));
    vertstartind = (size(newim,1)-size(tempim,1))/2;
    vertinds = vertstartind+(1:size(tempim,1));
    newim(vertinds,:) = tempim;
    newim(1:vertinds(1)-1,:) = repmat(tempim(1,:), vertstartind,1);
    newim(vertinds(end)+1:end,:) = repmat(tempim(end,:), size(newim,1)-(size(tempim,1)+vertstartind),1);
    imwrite(newim, [sizeCIsavepath filesep tiflist(ii).name])
    end
end

% order of blocks: ICwcfg1, ICwcfg0, ICkcfg1, ICkcfg0, RFCI, sizeCI

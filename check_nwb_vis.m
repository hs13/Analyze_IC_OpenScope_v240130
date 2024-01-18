%% check vis
% changes in visblock keys 
% in IC blocks, "frame" -> "Image"
% in RFCI and sizeCI block, viskeys "orientation" -> "Ori"
% in RFCI block, "x_position" -> "Pos_x", "y_position" -> "Pos_y"
% request adding fields: trialorder, trialtypedescription,
% MaskDiaVisDeg (RFCI_presentations and sizeCI_presentations),
% RFcentersVisDeg (RFCI_presentations)

nwbspikefile = "G:\My Drive\RESEARCH\IllusionOpenScope\sub-625554_ses-1181330601_ogen.nwb";
nwb = nwbRead(nwbspikefile); 

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

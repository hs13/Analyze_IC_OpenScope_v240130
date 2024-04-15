datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

neuopt = 'RS';
svmdesc = 'trainREx';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
whichblock = [whichICblock '_presentations'];
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
probes = {'A', 'B', 'C', 'D', 'E', 'F'};

%%
for ises = 1:numel(nwbsessions)
    clearvars -except datadir nwbsessions ises svmdesc preproc whichSVMkernel whichICblock whichblock visareas probes 
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'SVM_' svmdesc '_selectareas' filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];
% load([pathpp 'spiketimes.mat']) %, 'neuctx')

    load([pathpp 'postprocessed.mat'])
    load([pathpp 'postprocessed_probeC.mat'], 'psthtli', 'vis')
    psthall = struct();
    psthall.(whichblock) = false(length(psthtli), length(vis.(whichblock).trialorder), length(neuallloc));
    for iprobe = 1:numel(probes)
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        psthall.(whichblock)(:,:,neuoind) = psth.(whichblock);
    end
    tloi = psthtli>0 & psthtli<=400;
    tempRall = 1000*squeeze(mean(psthall.(whichblock)(tloi,:,:),1));
    if ~isequal(Rall.(whichblock), tempRall)
        error('check psthall')
    end
    
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        fprintf('%d/%d %s %s\n', ises, numel(nwbsessions), svmdesc, whichvisarea)
        tic

        svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
        svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
        if ~exist(svmmdlfn, 'file')
            disp([svmmdlfn ' does not exist'])
            continue
        end
        load(svmfn)
        load(svmmdlfn)
        % load([pathpp 'postprocessed_probe' probes{iprobe} '.mat'])

        switch svmdesc
            case 'trainICRC'
                SVMout = SVMtrainICRC;
                SVM_models = SVMtrainICRC_models;
            case 'trainREx'
                SVMout = SVMtrainREx;
                SVM_models = SVMtrainREx_models;
        end

        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        traintrialtypes = SVMout.trialtypes;
        
        if length(SVMout.neu2anal) ~= size(psthall.(whichblock),3)
            error('check match between pathpp and loaded SVMout')
        end
        
        neu2anal = SVMout.neu2anal;
        Nclasses = length(traintrialtypes);
        Nsplits = numel(SVM_models.spkcnt);

        % SVM PSTH
        Twin = 5;
        psthbinTends = Twin:Twin:400;

        trialorder = ICtrialtypes(vis.(whichblock).trialorder + 1);
        temppsth = psthall.(whichblock)(:,:,neu2anal);
        Nbins = size(psthbinTends,2);
        cumpsthbin = NaN(Nbins, length(trialorder), nnz(neu2anal) );
        for ibin = 1:Nbins
            temptli = psthtli>0 & psthtli<=psthbinTends(ibin);
            cumpsthbin(ibin, :, :) = sum(temppsth(temptli, :, :), 1)*1000/400;
        end
        tempR = Rall.(whichblock)(:,neu2anal);

        % figure; hold all
        % plot(psthtli,1000*mean(temppsth,[2,3]))
        % plot(psthtli(psthbinTinds(13,:)),mean(psthbin,[2,3]), 'linewidth', 1)

        SVMcumpsth = struct();
        SVMcumpsth.Twin = Twin;
        SVMcumpsth.trialorder = trialorder;
        SVMcumpsth.traintrialinds = SVMout.spkcnt.traintrialinds;
        SVMcumpsth.testtrialinds = SVMout.spkcnt.testtrialinds;
        SVMcumpsth.traintrialtypes = SVMout.trialtypes;
        SVMcumpsth.Ylabs = SVMout.spkcnt.Ylabs;
        SVMcumpsth.psthtli = psthtli;
        SVMcumpsth.psthbinTends = psthbinTends;
        SVMcumpsth.cumpsthbin = cumpsthbin;
        SVMcumpsth.label = NaN(Nbins, length(trialorder), Nsplits);
        SVMcumpsth.score = NaN(Nbins, length(trialorder), Nclasses, Nsplits);
        for isplit = 1:Nsplits
            traintrialinds = SVMout.spkcnt.traintrialinds(:,isplit);
            testtrialinds = SVMout.spkcnt.testtrialinds(:,isplit);
            switch preproc
                case 'none'
                    normppsthbin = cumpsthbin;
                case 'zscore'
                    % Z-score
                    trainRmean = reshape( mean(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    trainRstd = reshape( std(tempR(traintrialinds,:),0,1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (cumpsthbin-trainRmean)./trainRstd;
                    normppsthbin(isnan(normppsthbin))=0;
                case 'minmax'
                    trainRmin = reshape( min(tempR(traintrialinds,:),[],1), 1,1,nnz(neu2anal));
                    trainRrange = reshape( range(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (cumpsthbin-trainRmin)./trainRrange;
                case 'meancenter'
                    trainRmean = reshape( mean(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (cumpsthbin-trainRmean);
            end
            
            tempSVMmodel = SVM_models.spkcnt{isplit};
            for ibin = 1:Nbins
                Xtemp = squeeze( normppsthbin(ibin,:,:) );
                [templabel,tempscore] = predict(tempSVMmodel,Xtemp); % Xtemp is Ntrials X Nneurons
                SVMcumpsth.label(ibin,:,isplit) = templabel;
                SVMcumpsth.score(ibin,:,:,isplit) = tempscore;
            end
        end
        %{
        isequal(tempR, squeeze(cumpsthbin(end,:,:)) )
        Xtemp = ( tempR-mean(tempR(traintrialinds,:),1) )./std(tempR(traintrialinds,:),0,1);
        Xtemp(isnan(Xtemp))=0;
        [tempRlabel,tempRscore] = predict(tempSVMmodel,Xtemp); % Xtemp is Ntrials X Nneurons
        isequal(squeeze(SVMout.spkcnt.all.label(:,isplit)), tempRlabel)
        isequal(squeeze(SVMout.spkcnt.all.score(:,:,isplit)), tempRscore)
        figure; plot( tempRscore(:), reshape(squeeze(SVMout.spkcnt.all.score(:,:,isplit)),[],1), '.')
        % SVMout score is roughly twice that of the new calculation... but
        % only in MAC!!! (in Windows, where it was originally calculated,
        % these two almost match, although not exactly...
        
        tempscore = squeeze(SVMcumpsth.score(end,:,:,:));
        isequal(2*tempscore, SVMout.spkcnt.all.score)
        max(abs(2*tempscore(:)-SVMout.spkcnt.all.score(:)))
        %}
        save([pathsvm, 'SVMcumpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, neuopt, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], ...
            'SVMcumpsth', '-v7.3')
        toc
    end
end

%%


SVMcumpsthall.V1 = SVMcumpsth;

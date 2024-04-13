datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

svmdesc = 'trainREx';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
whichblock = [whichICblock '_presentations'];
visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam', 'LGd', 'LP'};
probes = {'A', 'B', 'C', 'D', 'E', 'F'};

%%
for ises = 1:numel(nwbsessions)
    clearvars -except datadir nwbsessions ises svmdesc preproc whichSVMkernel whichICblock probes 
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
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

        svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
        svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
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
        Twin = 25; % bins 5, 25, 100, 400ms
        if Twin<=50
            psthbinTinds = (0:Twin-1)'+(find(psthtli==-100):5:find(psthtli==500)-Twin);
        else
            psthbinTinds = (0:Twin-1)'+(find(psthtli==-200):10:find(psthtli==600)-Twin);
        end

        trialorder = ICtrialtypes(vis.(whichblock).trialorder + 1);
        temppsth = psthall.(whichblock)(:,:,neu2anal);
        Nbins = size(psthbinTinds,2);
        psthbin = NaN(Nbins, length(trialorder), nnz(neu2anal) );
        for ibin = 1:Nbins
            temptli = psthbinTinds(:,ibin);
            psthbin(ibin, :, :) = 1000*mean(temppsth(temptli, :, :), 1);
        end
        tempR = Rall.(whichblock)(:,neu2anal);

        % figure; hold all
        % plot(psthtli,1000*mean(temppsth,[2,3]))
        % plot(psthtli(psthbinTinds(13,:)),mean(psthbin,[2,3]), 'linewidth', 1)

        SVMpsth = struct();
        SVMpsth.Twin = Twin;
        SVMpsth.trialorder = trialorder;
        SVMpsth.traintrialinds = SVMout.spkcnt.traintrialinds;
        SVMpsth.testtrialinds = SVMout.spkcnt.testtrialinds;
        SVMpsth.traintrialtypes = SVMout.trialtypes;
        SVMpsth.Ylabs = SVMout.spkcnt.Ylabs;
        SVMpsth.psthtli = psthtli;
        SVMpsth.psthbinTinds = psthbinTinds;
        SVMpsth.psthbin = psthbin;
        SVMpsth.label = NaN(Nbins, length(trialorder), Nsplits);
        SVMpsth.score = NaN(Nbins, length(trialorder), Nclasses, Nsplits);
        for isplit = 1:Nsplits
            traintrialinds = SVMout.spkcnt.traintrialinds(:,isplit);
            testtrialinds = SVMout.spkcnt.testtrialinds(:,isplit);
            switch preproc
                case 'none'
                    normppsthbin = psthbin;
                case 'zscore'
                    % Z-score
                    trainRmean = reshape( mean(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    trainRstd = reshape( std(tempR(traintrialinds,:),0,1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (psthbin-trainRmean)./trainRstd;
                    normppsthbin(isnan(normppsthbin))=0;
                case 'minmax'
                    trainRmin = reshape( min(tempR(traintrialinds,:),[],1), 1,1,nnz(neu2anal));
                    trainRrange = reshape( range(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (psthbin-trainRmin)./trainRrange;
                case 'meancenter'
                    trainRmean = reshape( mean(tempR(traintrialinds,:),1), 1,1,nnz(neu2anal));
                    
                    normppsthbin = (psthbin-trainRmean);
            end
            
            for ibin = 1:Nbins
                Xtemp = squeeze( normppsthbin(ibin,:,:) );
                tempSVMmodel = SVM_models.spkcnt{isplit};
                [templabel,tempscore] = predict(tempSVMmodel,Xtemp); % Xtemp is Ntrials X Nneurons
                SVMpsth.label(ibin,:,isplit) = templabel;
                SVMpsth.score(ibin,:,:,isplit) = tempscore;
            end
        end
        save([pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], ...
            'SVMpsth', '-v7.3')
        toc
    end
end

%% decoder weight SVMbeta
% weight_vector=c1.Beta;
% bais_vector=c1.Bias;
for ises = 1:numel(nwbsessions)
    clearvars -except datadir nwbsessions ises svmdesc preproc whichSVMkernel whichICblock probes 
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];
    load([pathpp 'spiketimes.mat'], 'neuctx')
    
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        fprintf('%d/%d %s %s\n', ises, numel(nwbsessions), svmdesc, whichvisarea)
        tic

        svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
        svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
        if ~exist(svmmdlfn, 'file')
            disp([svmmdlfn ' does not exist'])
            continue
        end
        load(svmfn)
        load(svmmdlfn)

        switch svmdesc
            case 'trainICRCtestRE'
                SVMout = SVMtrainICRC;
                SVM_models = SVMtrainICRC_models;
            case 'trainRExtestICRC'
                SVMout = SVMtrainREx;
                SVM_models = SVMtrainREx_models;
        end

        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        traintrialtypes = SVMout.trialtypes;
        c = nchoosek(1:length(traintrialtypes),2);
        traintrialcombos = traintrialtypes(c);

        neu2anal = SVMout.neu2anal;
        Nclasses = length(traintrialtypes);
        Nbinarylearners = size(traintrialcombos,1);
        Nsplits = numel(SVM_models.spkcnt);

        SVMcodingname = cell(Nsplits,1);
        Nlearners = zeros(Nsplits,1);
        betalearners = cell(Nsplits,1);
        for isplit = 1:Nsplits
            Ylabs = SVMout.spkcnt.Ylabs;
            splitcombo = sort(Ylabs(c),2);
            Nlearners(isplit) = numel(SVM_models.spkcnt{isplit}.BinaryLearners);
            SVMcodingname{isplit} = SVM_models.spkcnt{isplit}.CodingName;
            betalearners{isplit} = NaN(nnz(neu2anal), Nlearners(isplit));
            for ibl = 1:Nlearners(isplit)
                switch SVM_models.spkcnt{isplit}.CodingName
                    case 'onevsone'
                        indbl = find(ismember(splitcombo, traintrialcombos(ibl,:), 'rows' ));
                        if ~isequal(sort(Ylabs(c(indbl,:))), traintrialcombos(ibl,:))
                            error('sanity check did not pass, check indbl')
                        end
                        if isequal(Ylabs(c(indbl,:)), traintrialcombos(ibl,:)) % same order
                            betagain = 1;
                        else % flipped order
                            betagain = -1;
                        end
                    case 'onevsall'
                        if length(traintrialtypes)==2
                            indbl = 1;
                            if isequal(Ylabs, traintrialtypes) % same order
                                betagain = 1;
                            else % flipped order
                                betagain = -1;
                            end
                        else
                            indbl = find(Ylabs==traintrialtypes(ibl));
                            betagain = 1;
                        end
                    otherwise
                        error('unexpected CodingName -- need to make another exception case')
                end
                betalearners{isplit}(:,ibl) = betagain * SVM_models.spkcnt{isplit}.BinaryLearners{indbl}.Beta;
            end
        end
        save([pathsvm, 'SVMbeta_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], ...
            'ICtrialtypes', 'traintrialtypes', 'traintrialcombos', 'neuctxinprobe', 'SVMcodingname', 'Nlearners', 'betalearners', '-v7.3')

        toc
    end
end

%%

betabl = cat(3,betalearners{:});

for ibl = 1:Nbinarylearners
    corrbeta = corr(squeeze(betabl(:,ibl,Nlearners==Nbinarylearners)));
disp([ibl mean(corrbeta(triu(true(size(corrbeta)), 1))) mean(abs(corrbeta(triu(true(size(corrbeta)), 1))))])
end

figure; histogram( corrbeta(triu(true(size(corrbeta)))), -1:0.05:1 )
figure; imagesc(corrbeta); caxis([-1 1]); colorbar
figure; plot(betabl(:,1,1), betabl(:,1,6), '.')

Z = linkage(corrbeta,'complete','correlation');
figure
[H,T,outperm] = dendrogram(Z);

figure; imagesc(corrbeta(outperm,outperm)); caxis([-1 1]); colorbar

Z2 = linkage(squeeze(betabl(:,ibl,Nlearners=Nbinarylearners))','complete','correlation');
figure
[H2,T2,outperm2] = dendrogram(Z2);

figure; imagesc(corrbeta(outperm2,outperm2)); caxis([-1 1]); colorbar

% note, Z and Z2 are very similar

figure
hold all
for itt = 1:numel(traintrialtypes)
    typi = traintrialtypes(itt);
    trialsoi = trialorder==typi;
    temppsth = mean(SVMpsth.label(:,trialsoi,:)==typi,[2,3]);
    plot(psthtli(psthbinTinds(51,:)), temppsth)
end

%%
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])

load([drivepath 'RESEARCH/logmean_logvar/OpenScope_UMAP_V1RS_lmlvslopes.mat'])

trialorder = trialorderacc{1};
hireptt = unique(trialorder);
Nhireptt = numel(hireptt);
Ndims = 2;
Nsplits = 2;
Nsessions = numel(nwbsessions);
Nslopes = numel(lmlvslope_list);

trialorderind= trialorder;
for ii = 1:Nhireptt
    trialorderind(trialorder==hireptt(ii))=ii;
end
trialcol = zeros(Nhireptt,3);
trialcol(hireptt==0,:) = [0 0 0];
trialcol(hireptt==101,:) = [1 0 0];
trialcol(hireptt==105,:) = [0 1 1];
trialcol(hireptt==109,:) = [0 0.5 0.5];
trialcol(hireptt==106,:) = [0 1 0];
trialcol(hireptt==111,:) = [0 0.5 0];
trialcol(hireptt==107,:) = [1 0.5 0];
trialcol(hireptt==110,:) = [0.5 0.25 0];
trialcol(hireptt==1105,:) = [0 0 1];
trialcol(hireptt==1109,:) = [0 0 0.5];
trialcol(hireptt==1201,:) = [1 0 1];
trialcol(hireptt==1299,:) = [0.5 0 0.5];

fs = 14;

%% centroid distance on UMAP between pairs of trial types

umapcentroid = struct();
umapcentroid.unsup_origall = NaN(Nhireptt, Ndims, numel(nwbsessions));
umapcentroid.unsup_origtest = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.unsup_origtrain = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.semisup_origall = NaN(Nhireptt, Ndims, numel(nwbsessions));
umapcentroid.semisup_origtest = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
umapcentroid.semisup_origtrain = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);

umapcentroid.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umapcentroid.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umapcentroid.unsup_lmlvtrain = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
umapcentroid.semisup_lmlvtrain = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umapcentroid.unsup_lmlvall{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions));
    umapcentroid.unsup_lmlvtest{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.unsup_lmlvtrain{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.semisup_lmlvall{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions));
    umapcentroid.semisup_lmlvtest{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
    umapcentroid.semisup_lmlvtrain{islope} = NaN(Nhireptt, Ndims, numel(nwbsessions), Nsplits);
end

umapctrdist = struct();
umapctrdist.unsup_origall = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
umapctrdist.unsup_origtest = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.unsup_origtrain = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.semisup_origall = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
umapctrdist.semisup_origtest = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
umapctrdist.semisup_origtrain = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);

umapctrdist.unsup_lmlvall = cell(1,numel(lmlvslope_list));
umapctrdist.unsup_lmlvtest = cell(1,numel(lmlvslope_list));
umapctrdist.unsup_lmlvtrain = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvall = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvtest = cell(1,numel(lmlvslope_list));
umapctrdist.semisup_lmlvtrain = cell(1,numel(lmlvslope_list));
for islope = 1:numel(lmlvslope_list)
    umapctrdist.unsup_lmlvall{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
    umapctrdist.unsup_lmlvtest{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.unsup_lmlvtrain{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.semisup_lmlvall{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions));
    umapctrdist.semisup_lmlvtest{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
    umapctrdist.semisup_lmlvtrain{islope} = NaN(Nhireptt, Nhireptt, numel(nwbsessions), Nsplits);
end

for ises = 1:numel(nwbsessions)
    tic
    for ieuc = 1:2
        switch ieuc
            case 1
                eucfield = 'semisup_origall';
                tempembedding = UMAPorigallagg(ises).embeddings;
                temptrialorder = trialorderacc{ises};
            case 2
                eucfield = 'unsup_origall';
                tempembedding = UMAPorigall_unsupagg(ises).embeddings;
                temptrialorder = trialorderacc{ises};
        end
        tempcentroid = NaN(Nhireptt, Ndims);
        for ii = 1:Nhireptt
            trialsoi = temptrialorder==hireptt(ii);
            tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
        end
        umapcentroid.(eucfield)(:,:,ises) = tempcentroid;
        
        tempctrdist = zeros(Nhireptt,Nhireptt);
        for idim = 1:Ndims
            tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
        end
        tempctrdist = sqrt(tempctrdist);
        umapctrdist.(eucfield)(:,:,ises) = tempctrdist;
        
    end
    
    for ieuc = 1:4
        for isplit = 1:Nsplits
            switch ieuc
                case 1
                    eucfield = 'semisup_origtest';
                    tempembedding = squeeze(UMAPorigagg(ises).test_embeddings(isplit,:,:));
                    temptrialorder = UMAPorigagg(ises).test_truelabels(isplit,:);
                case 2
                    eucfield = 'semisup_origtrain';
                    tempembedding = squeeze(UMAPorigagg(ises).train_embeddings(isplit,:,:));
                    temptrialorder = trialorderacc{ises}( 1+UMAPorigagg(ises).train_index(isplit,:) );
                case 3
                    eucfield = 'unsup_origtest';
                    tempembedding = squeeze(UMAPorig_unsupagg(ises).test_embeddings(isplit,:,:));
                    temptrialorder = UMAPorig_unsupagg(ises).test_truelabels(isplit,:);
                case 4
                    eucfield = 'unsup_origtrain';
                    tempembedding = squeeze(UMAPorig_unsupagg(ises).train_embeddings(isplit,:,:));
                    temptrialorder = trialorderacc{ises}( 1+UMAPorig_unsupagg(ises).train_index(isplit,:) );
            end
            tempcentroid = NaN(Nhireptt, Ndims);
            for ii = 1:Nhireptt
                trialsoi = temptrialorder==hireptt(ii);
                tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
            end
            umapcentroid.(eucfield)(:,:,ises,isplit) = tempcentroid;
            
            tempctrdist = zeros(Nhireptt,Nhireptt);
            for idim = 1:Ndims
                tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
            end
            tempctrdist = sqrt(tempctrdist);
            umapctrdist.(eucfield)(:,:,ises,isplit) = tempctrdist;
        end
    end
    
    for islope = 1:numel(lmlvslope_list)
        for ieuc = 1:2
            switch ieuc
                case 1
                    eucfield = 'semisup_lmlvall';
                    tempembedding = UMAPlmlvallacc{ises,islope}.embeddings;
                    temptrialorder = trialorderacc{ises};
                case 2
                    eucfield = 'unsup_lmlvall';
                    tempembedding = UMAPlmlvall_unsupacc{ises,islope}.embeddings;
                    temptrialorder = trialorderacc{ises};
            end
            tempcentroid = NaN(Nhireptt, Ndims);
            for ii = 1:Nhireptt
                trialsoi = temptrialorder==hireptt(ii);
                tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
            end
            umapcentroid.(eucfield){islope}(:,:,ises) = tempcentroid;
            
            tempctrdist = zeros(Nhireptt,Nhireptt);
            for idim = 1:Ndims
                tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
            end
            tempctrdist = sqrt(tempctrdist);
            umapctrdist.(eucfield){islope}(:,:,ises) = tempctrdist;
            
        end
        
        for ieuc = 1:4
            for isplit = 1:Nsplits
                switch ieuc
                    case 1
                        eucfield = 'semisup_lmlvtest';
                        tempembedding = squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(isplit,:,:));
                        temptrialorder = UMAPlmlvacc{ises,islope}.test_truelabels(isplit,:);
                    case 2
                        eucfield = 'semisup_lmlvtrain';
                        tempembedding = squeeze(UMAPlmlvacc{ises,islope}.train_embeddings(isplit,:,:));
                        temptrialorder = trialorderacc{ises}( 1+UMAPlmlvacc{ises,islope}.train_index(isplit,:) );
                    case 3
                        eucfield = 'unsup_lmlvtest';
                        tempembedding = squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,:,:));
                        temptrialorder = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(isplit,:);
                    case 4
                        eucfield = 'unsup_lmlvtrain';
                        tempembedding = squeeze(UMAPlmlv_unsupacc{ises,islope}.train_embeddings(isplit,:,:));
                        temptrialorder = trialorderacc{ises}( 1+UMAPlmlv_unsupacc{ises,islope}.train_index(isplit,:) );
                end
                tempcentroid = NaN(Nhireptt, Ndims);
                for ii = 1:Nhireptt
                    trialsoi = temptrialorder==hireptt(ii);
                    tempcentroid(ii,:) = mean(tempembedding(trialsoi,:),1);
                end
                umapcentroid.(eucfield){islope}(:,:,ises,isplit) = tempcentroid;
                
                tempctrdist = zeros(Nhireptt,Nhireptt);
                for idim = 1:Ndims
                    tempctrdist = tempctrdist + ( tempcentroid(:,idim)-tempcentroid(:,idim)' ).^2;
                end
                tempctrdist = sqrt(tempctrdist);
                umapctrdist.(eucfield){islope}(:,:,ises,isplit) = tempctrdist;
                
            end
        end
    end
    toc
end

%% calculate overlap between UMAPs after procrustes transformation: between sessions, UMAP on all trials
% [d,Z,transform] = procrustes(X,Y) also returns the transformation that maps Y to Z.
% transform is a structure array with fields:
% c — Translation component
% T — Orthogonal rotation and reflection component
% b — Scale component
%
% when X and Y are Npoints X Ndim
% c = transform.c; % Npoints X Ndim
% T = transform.T; % Ndim X Ndim
% b = transform.b; % scalar
%
% Z = b*Y*T + c;

% [k,v] = boundary(x,y)

warning('off')

for supind = 1:2
    tic
    switch supind
        case 1
            UMAPopt = 'semisup';
            UMAPoption = 'semisup_origall';
        case 2
            UMAPopt = 'unsup';
            UMAPoption = 'unsup_origall';
    end
    
    UMAPareasall.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPtfareasall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt); % row is template index, col is transformed umap index
    UMAPoverlapall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPoverlapportionall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPoverlaptfportionall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPnormoverlapall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPnormoverlaptfall.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    
    % transform between semisupervised all trials sessions 2 and 9
    for ises = 1:Nsessions
        for jses = 1:Nsessions
            
            switch UMAPoption
                case 'semisup_origall'
                    umap1 = double( UMAPorigallagg(ises).embeddings );
                    umap2 = double( UMAPorigallagg(jses).embeddings );
                    
                    umapctr1 = umapcentroid.semisup_origall(:,:,ises);
                    umapctr2 = umapcentroid.semisup_origall(:,:,jses);
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
                case 'unsup_origall'
                    umap1 = double( UMAPorigall_unsupagg(ises).embeddings );
                    umap2 = double( UMAPorigall_unsupagg(jses).embeddings );
                    
                    umapctr1 = umapcentroid.unsup_origall(:,:,ises);
                    umapctr2 = umapcentroid.unsup_origall(:,:,jses);
                    
                    trialord1 = trialorderacc{ises};
                    trialord2 = trialorderacc{jses};
            end
            
            [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
            
            if any(isnan(UMAPareasall.(UMAPopt)(ises,:)))
                UMAPareasall.(UMAPopt)(ises,:) = templateareas;
            else
                if ~isequal(UMAPareasall.(UMAPopt)(ises,:)', templateareas(:))
                    error('UMAPareas unexpected mismatch: please check')
                end
            end
            UMAPtfareasall.(UMAPopt)(ises,jses,:) = transformedareas;
            UMAPoverlapall.(UMAPopt)(ises,jses,:) = diag(overlapareas);
            UMAPoverlapportionall.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(templateareas,[],1);
            UMAPoverlaptfportionall.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(transformedareas,[],1);
            
            UMAPnormoverlapall.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
            UMAPnormoverlaptfall.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
        end
    end
    toc
end

UMAPareas_lmlvall = struct();
UMAPtfareas_lmlvall = struct();
UMAPoverlap_lmlvall = struct();
UMAPoverlapportion_lmlvall = struct();
UMAPoverlaptfportion_lmlvall = struct();
UMAPnormoverlap_lmlvall = struct();
UMAPnormoverlaptf_lmlvall = struct();
for islope = 1:Nslopes
    tic
    for supind = 1:2
        switch supind
            case 1
                UMAPopt = 'semisup';
                UMAPoption = 'semisup_lmlvall';
            case 2
                UMAPopt = 'unsup';
                UMAPoption = 'unsup_lmlvall';
        end
        
        UMAPareas_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPtfareas_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt); % row is template index, col is transformed umap index
        UMAPoverlap_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPoverlapportion_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPoverlaptfportion_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPnormoverlap_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPnormoverlaptf_lmlvall(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        
        % transform between semisupervised all trials sessions 2 and 9
        for ises = 1:Nsessions
            for jses = 1:Nsessions
                
                switch UMAPoption
                    case 'semisup_lmlvall'
                        umap1 = double( UMAPlmlvallacc{ises, islope}.embeddings );
                        umap2 = double( UMAPlmlvallacc{jses, islope}.embeddings );
                        
                        umapctr1 = umapcentroid.semisup_lmlvall{islope}(:,:,ises);
                        umapctr2 = umapcentroid.semisup_lmlvall{islope}(:,:,jses);
                        
                        trialord1 = trialorderacc{ises};
                        trialord2 = trialorderacc{jses};
                    case 'unsup_lmlvall'
                        umap1 = double( UMAPlmlvall_unsupacc{ises, islope}.embeddings );
                        umap2 = double( UMAPlmlvall_unsupacc{jses, islope}.embeddings );
                        
                        umapctr1 = umapcentroid.unsup_lmlvall{islope}(:,:,ises);
                        umapctr2 = umapcentroid.unsup_lmlvall{islope}(:,:,jses);
                        
                        trialord1 = trialorderacc{ises};
                        trialord2 = trialorderacc{jses};
                end
                
                [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
                
                if any(isnan(UMAPareas_lmlvall(islope).(UMAPopt)(ises,:)))
                    UMAPareas_lmlvall(islope).(UMAPopt)(ises,:) = templateareas;
                else
                    if ~isequal(UMAPareas_lmlvall(islope).(UMAPopt)(ises,:)', templateareas(:))
                        error('UMAPareas unexpected mismatch: please check')
                    end
                end
                UMAPtfareas_lmlvall(islope).(UMAPopt)(ises,jses,:) = transformedareas;
                UMAPoverlap_lmlvall(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas);
                UMAPoverlapportion_lmlvall(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(templateareas,[],1);
                UMAPoverlaptfportion_lmlvall(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(transformedareas,[],1);
                UMAPnormoverlap_lmlvall(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
                UMAPnormoverlaptf_lmlvall(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
            end
        end
    end
    toc
end


% calculate overlap between UMAPs after procrustes transformation: between sessions, cross-validated UMAP test trials
isplit = 1;
for supind = 1:2
    tic
    switch supind
        case 1
            UMAPopt = 'semisup';
            UMAPoption = 'semisup_origtest';
        case 2
            UMAPopt = 'unsup';
            UMAPoption = 'unsup_origtest';
    end
    
    UMAPareastest.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPtfareastest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt); % row is template index, col is transformed umap index
    UMAPoverlaptest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPoverlapportiontest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPoverlaptfportiontest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPnormoverlaptest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    UMAPnormoverlaptftest.(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
    
    % transform between semisupervised all trials sessions 2 and 9
    for ises = 1:Nsessions
        for jses = 1:Nsessions
            
            switch UMAPoption
                case 'semisup_origtest'
                    umap1 = double( squeeze(UMAPorigagg(ises).test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPorigagg(jses).test_embeddings(isplit,:,:)) );
                    
                    umapctr1 = umapcentroid.semisup_origtest(:,:,ises,isplit);
                    umapctr2 = umapcentroid.semisup_origtest(:,:,jses,isplit);
                    
                    trialord1 = UMAPorigagg(ises).test_truelabels(isplit,:);
                    trialord2 = UMAPorigagg(jses).test_truelabels(isplit,:);
                case 'unsup_origtest'
                    umap1 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(isplit,:,:)) );
                    umap2 = double( squeeze(UMAPorig_unsupagg(jses).test_embeddings(isplit,:,:)) );
                    
                    umapctr1 = umapcentroid.unsup_origtest(:,:,ises,isplit);
                    umapctr2 = umapcentroid.unsup_origtest(:,:,jses,isplit);
                    
                    trialord1 = UMAPorig_unsupagg(ises).test_truelabels(isplit,:);
                    trialord2 = UMAPorig_unsupagg(jses).test_truelabels(isplit,:);
            end
            
            [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
            
            if any(isnan(UMAPareastest.(UMAPopt)(ises,:)))
                UMAPareastest.(UMAPopt)(ises,:) = templateareas;
            else
                if ~isequal(UMAPareastest.(UMAPopt)(ises,:)', templateareas(:))
                    error('UMAPareas unexpected mismatch: please check')
                end
            end
            UMAPtfareastest.(UMAPopt)(ises,jses,:) = transformedareas;
            UMAPoverlaptest.(UMAPopt)(ises,jses,:) = diag(overlapareas);
            UMAPoverlapportiontest.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(templateareas,[],1);
            UMAPoverlaptfportiontest.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(transformedareas,[],1);
            UMAPnormoverlaptest.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
            UMAPnormoverlaptftest.(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
        end
    end
    toc
end

UMAPareas_lmlvtest = struct();
UMAPtfareas_lmlvtest = struct();
UMAPoverlap_lmlvtest = struct();
UMAPoverlapportion_lmlvtest = struct();
UMAPoverlaptfportion_lmlvtest = struct();
UMAPnormoverlap_lmlvtest = struct();
UMAPnormoverlaptf_lmlvtest = struct();
for islope = 1:Nslopes
    for supind = 1:2
        tic
        switch supind
            case 1
                UMAPopt = 'semisup';
                UMAPoption = 'semisup_lmlvtest';
            case 2
                UMAPopt = 'unsup';
                UMAPoption = 'unsup_lmlvtest';
        end
        
        UMAPareas_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPtfareas_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt); % row is template index, col is transformed umap index
        UMAPoverlap_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPoverlapportion_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPoverlaptfportion_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPnormoverlap_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        UMAPnormoverlaptf_lmlvtest(islope).(UMAPopt) = NaN(Nsessions, Nsessions, Nhireptt);
        
        % transform between semisupervised all trials sessions 2 and 9
        for ises = 1:Nsessions
            for jses = 1:Nsessions
                
                switch UMAPoption
                    case 'semisup_lmlvtest'
                        umap1 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(isplit,:,:)) );
                        umap2 = double( squeeze(UMAPlmlvacc{jses,islope}.test_embeddings(isplit,:,:)) );
                        
                        umapctr1 = umapcentroid.semisup_lmlvtest{islope}(:,:,ises,isplit);
                        umapctr2 = umapcentroid.semisup_lmlvtest{islope}(:,:,jses,isplit);
                        
                        trialord1 = UMAPlmlvacc{ises,islope}.test_truelabels(isplit,:);
                        trialord2 = UMAPlmlvacc{jses,islope}.test_truelabels(isplit,:);
                    case 'unsup_lmlvtest'
                        umap1 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(isplit,:,:)) );
                        umap2 = double( squeeze(UMAPlmlv_unsupacc{jses,islope}.test_embeddings(isplit,:,:)) );
                        
                        umapctr1 = umapcentroid.unsup_lmlvtest{islope}(:,:,ises,isplit);
                        umapctr2 = umapcentroid.unsup_lmlvtest{islope}(:,:,jses,isplit);
                        
                        trialord1 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(isplit,:);
                        trialord2 = UMAPlmlv_unsupacc{jses,islope}.test_truelabels(isplit,:);
                end
                
                [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
                
                if any(isnan(UMAPareas_lmlvtest(islope).(UMAPopt)(ises,:)))
                    UMAPareas_lmlvtest(islope).(UMAPopt)(ises,:) = templateareas;
                else
                    if ~isequal(UMAPareas_lmlvtest(islope).(UMAPopt)(ises,:)', templateareas(:))
                        error('UMAPareas unexpected mismatch: please check')
                    end
                end
                UMAPtfareas_lmlvtest(islope).(UMAPopt)(ises,jses,:) = transformedareas;
                UMAPoverlap_lmlvtest(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas);
                UMAPoverlapportion_lmlvtest(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(templateareas,[],1);
                UMAPoverlaptfportion_lmlvtest(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(transformedareas,[],1);
                UMAPnormoverlap_lmlvtest(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
                UMAPnormoverlaptf_lmlvtest(islope).(UMAPopt)(ises,jses,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
            end
        end
        toc
    end
end


% calculate overlap between UMAPs after procrustes transformation: between test sets within session
for supind = 1:2
    tic
    switch supind
        case 1
            UMAPopt = 'semisup';
            UMAPoption = 'semisup_origcv';
        case 2
            UMAPopt = 'unsup';
            UMAPoption = 'unsup_origcv';
    end
    
    UMAPareascv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPtfareascv.(UMAPopt) = NaN(Nsessions, Nhireptt); % row is template index, col is transformed umap index
    UMAPoverlapcv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPoverlapportioncv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPoverlaptfportioncv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPnormoverlapcv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    UMAPnormoverlaptfcv.(UMAPopt) = NaN(Nsessions, Nhireptt);
    
    % transform between semisupervised all trials sessions 2 and 9
    for ises = 1:Nsessions
        
        switch UMAPoption
            case 'semisup_origcv'
                umap1 = double( squeeze(UMAPorigagg(ises).test_embeddings(1,:,:)) );
                umap2 = double( squeeze(UMAPorigagg(ises).test_embeddings(2,:,:)) );
                
                umapctr1 = umapcentroid.semisup_origtest(:,:,ises,1);
                umapctr2 = umapcentroid.semisup_origtest(:,:,ises,2);
                
                trialord1 = UMAPorigagg(ises).test_truelabels(1,:);
                trialord2 = UMAPorigagg(ises).test_truelabels(2,:);
            case 'unsup_origcv'
                umap1 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(1,:,:)) );
                umap2 = double( squeeze(UMAPorig_unsupagg(ises).test_embeddings(2,:,:)) );
                
                umapctr1 = umapcentroid.unsup_origtest(:,:,ises,1);
                umapctr2 = umapcentroid.unsup_origtest(:,:,ises,2);
                
                trialord1 = UMAPorig_unsupagg(ises).test_truelabels(1,:);
                trialord2 = UMAPorig_unsupagg(ises).test_truelabels(2,:);
        end
        
        [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
        
        if any(isnan(UMAPareascv.(UMAPopt)(ises,:)))
            UMAPareascv.(UMAPopt)(ises,:) = templateareas;
        else
            if ~isequal(UMAPareascv.(UMAPopt)(ises,:)', templateareas(:))
                error('UMAPareas unexpected mismatch: please check')
            end
        end
        UMAPtfareascv.(UMAPopt)(ises,:) = transformedareas;
        UMAPoverlapcv.(UMAPopt)(ises,:) = diag(overlapareas);
        UMAPoverlapportioncv.(UMAPopt)(ises,:) = diag(overlapareas)./reshape(templateareas,[],1);
        UMAPoverlaptfportioncv.(UMAPopt)(ises,:) = diag(overlapareas)./reshape(transformedareas,[],1);
        UMAPnormoverlapcv.(UMAPopt)(ises,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
        UMAPnormoverlaptfcv.(UMAPopt)(ises,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
    end
    toc
end

UMAPareas_lmlvcv = struct();
UMAPtfareas_lmlvcv = struct();
UMAPoverlap_lmlvcv = struct();
UMAPoverlapportion_lmlvcv = struct();
UMAPoverlaptfportion_lmlvcv = struct();
UMAPnormoverlap_lmlvcv = struct();
UMAPnormoverlaptf_lmlvcv = struct();
for islope = 1:Nslopes
    for supind = 1:2
        tic
        switch supind
            case 1
                UMAPopt = 'semisup';
                UMAPoption = 'semisup_lmlvcv';
            case 2
                UMAPopt = 'unsup';
                UMAPoption = 'unsup_lmlvcv';
        end
        
        UMAPareas_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPtfareas_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt); % row is template index, col is transformed umap index
        UMAPoverlap_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPoverlapportion_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPoverlaptfportion_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPnormoverlap_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        UMAPnormoverlaptf_lmlvcv(islope).(UMAPopt) = NaN(Nsessions, Nhireptt);
        
        % transform between semisupervised all trials sessions 2 and 9
        for ises = 1:Nsessions
            
            switch UMAPoption
                case 'semisup_lmlvcv'
                    umap1 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPlmlvacc{ises,islope}.test_embeddings(2,:,:)) );
                    
                    umapctr1 = umapcentroid.semisup_lmlvtest{islope}(:,:,ises,1);
                    umapctr2 = umapcentroid.semisup_lmlvtest{islope}(:,:,ises,2);
                    
                    trialord1 = UMAPlmlvacc{ises,islope}.test_truelabels(1,:);
                    trialord2 = UMAPlmlvacc{ises,islope}.test_truelabels(2,:);
                case 'unsup_lmlvcv'
                    umap1 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(1,:,:)) );
                    umap2 = double( squeeze(UMAPlmlv_unsupacc{ises,islope}.test_embeddings(2,:,:)) );
                    
                    umapctr1 = umapcentroid.unsup_lmlvtest{islope}(:,:,ises,1);
                    umapctr2 = umapcentroid.unsup_lmlvtest{islope}(:,:,ises,2);
                    
                    trialord1 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(1,:);
                    trialord2 = UMAPlmlv_unsupacc{ises,islope}.test_truelabels(2,:);
            end
            
            [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2);
            
            if any(isnan(UMAPareas_lmlvcv(islope).(UMAPopt)(ises,:)))
                UMAPareas_lmlvcv(islope).(UMAPopt)(ises,:) = templateareas;
            else
                if ~isequal(UMAPareas_lmlvcv(islope).(UMAPopt)(ises,:)', templateareas(:))
                    error('UMAPareas unexpected mismatch: please check')
                end
            end
            UMAPtfareas_lmlvcv(islope).(UMAPopt)(ises,:) = transformedareas;
            UMAPoverlap_lmlvcv(islope).(UMAPopt)(ises,:) = diag(overlapareas);
            UMAPoverlapportion_lmlvcv(islope).(UMAPopt)(ises,:) = diag(overlapareas)./reshape(templateareas,[],1);
            UMAPoverlaptfportion_lmlvcv(islope).(UMAPopt)(ises,:) = diag(overlapareas)./reshape(transformedareas,[],1);
            UMAPnormoverlap_lmlvcv(islope).(UMAPopt)(ises,:) = diag(overlapareas)./reshape(sum(overlapareas,2),[],1);
            UMAPnormoverlaptf_lmlvcv(islope).(UMAPopt)(ises,:) = diag(overlapareas)./reshape(sum(overlapareas,1),[],1);
        end
        toc
    end
end

save([drivepath 'RESEARCH/logmean_logvar/OpenScope_UMAPoverlap_V1RS_lmlvslopes.mat'], 'umapcentroid', 'umapctrdist', ...
    'UMAPareasall', 'UMAPtfareasall', 'UMAPoverlapall', 'UMAPoverlapportionall', 'UMAPoverlaptfportionall', 'UMAPnormoverlapall', 'UMAPnormoverlaptfall', ...
    'UMAPareas_lmlvall', 'UMAPtfareas_lmlvall', 'UMAPoverlap_lmlvall', 'UMAPoverlapportion_lmlvall', 'UMAPoverlaptfportion_lmlvall', 'UMAPnormoverlap_lmlvall', 'UMAPnormoverlaptf_lmlvall', ...
    'UMAPareastest', 'UMAPtfareastest', 'UMAPoverlaptest', 'UMAPoverlapportiontest', 'UMAPoverlaptfportiontest', 'UMAPnormoverlaptest', 'UMAPnormoverlaptftest', ...
    'UMAPareas_lmlvtest', 'UMAPtfareas_lmlvtest', 'UMAPoverlap_lmlvtest', 'UMAPoverlapportion_lmlvtest', 'UMAPoverlaptfportion_lmlvtest', 'UMAPnormoverlap_lmlvtest', 'UMAPnormoverlaptf_lmlvtest', ...
    'UMAPareascv', 'UMAPtfareascv', 'UMAPoverlapcv', 'UMAPoverlapportioncv', 'UMAPoverlaptfportioncv', 'UMAPnormoverlapcv', 'UMAPnormoverlaptfcv', ...
    'UMAPareas_lmlvcv', 'UMAPtfareas_lmlvcv', 'UMAPoverlap_lmlvcv', 'UMAPoverlapportion_lmlvcv', 'UMAPoverlaptfportion_lmlvcv', 'UMAPnormoverlap_lmlvcv', 'UMAPnormoverlaptf_lmlvcv')

%%
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_UMAPoverlap_V1RS_lmlvslopes.mat'])

%% plot overlap between UMAPs after procrustes transformation: between sessions, UMAP on all trials
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', 'UMAP all trials: Proscrutes transform between sessions', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPoverlapportionall.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlapportion_lmlvall(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPoverlaptfportionall.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlaptfportion_lmlvall(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/transformed area')
    title(UMAPopt)
end

xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', 'UMAP all trials: Proscrutes transform between sessions', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPnormoverlapall.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlap_lmlvall(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPnormoverlaptfall.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlaptf_lmlvall(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap transformed area')
    title(UMAPopt)
end

%% calculate overlap between UMAPs after procrustes transformation: between sessions, cross-validated UMAP test trials

xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', sprintf('%d/2-fold cross-validated UMAP test trials: Proscrutes transform between sessions', isplit), 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPoverlapportiontest.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlapportion_lmlvtest(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPoverlaptfportiontest.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlaptfportion_lmlvtest(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/transformed area')
    title(UMAPopt)
    
end

xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', sprintf('%d/2-fold cross-validated UMAP test trials: Proscrutes transform between sessions', isplit), 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPnormoverlaptest.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlap_lmlvtest(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPnormoverlaptftest.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlaptf_lmlvtest(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap transformed area')
    title(UMAPopt)
    
end

%% calculate overlap between UMAPs after procrustes transformation: between test sets within session
xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', '2-fold cross-validated UMAP: Proscrutes transform between test sets', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPoverlapportioncv.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlapportion_lmlvcv(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPoverlaptfportioncv.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPoverlaptfportion_lmlvcv(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('overlap/transformed area')
    title(UMAPopt)
end


xl = [lmlvslope_list(1) lmlvslope_list(end)];
figure
annotation('textbox',[0.1 0.91 0.9 0.1], 'string', '2-fold cross-validated UMAP: Proscrutes transform between test sets', 'edgecolor', 'none', 'fontsize', fs)
for supind = 1:2
    switch supind
        case 1
            UMAPopt = 'semisup';
        case 2
            UMAPopt = 'unsup';
    end
    
    tempmat = mean(UMAPnormoverlapcv.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlap_lmlvcv(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap template area')
    title(UMAPopt)
    
    tempmat = mean(UMAPnormoverlaptfcv.(UMAPopt),3);
    tempmat(eye(Nsessions)==1) = NaN;
    tempmat = tempmat(:);
    overlapavgorig = nanmean(tempmat);
    overlapsemorig = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    
    overlapavglmlvs = NaN(1,Nslopes);
    overlapportionsemlmlvs = NaN(1,Nslopes);
    for islope = 1:Nslopes
        tempmat = mean(UMAPnormoverlaptf_lmlvcv(islope).(UMAPopt),3);
        tempmat(eye(Nsessions)==1) = NaN;
        tempmat = tempmat(:);
        overlapavglmlvs(islope) = nanmean(tempmat);
        overlapportionsemlmlvs(islope) = nanstd(tempmat,0,1)/sqrt(nnz(~isnan(tempmat)));
    end
    subplot(2,2,2+supind)
    hold all
    errorbar(lmlvslope_list, overlapavglmlvs, overlapportionsemlmlvs, 'ko-', 'linewidth', 1)
    plot(xl, overlapavgorig*[1 1], 'c-')
    xlabel('LMLV slopes')
    ylabel('norm-overlap transformed area')
    title(UMAPopt)
end


%% define function UMAPProscrutesOverlap
function [templateareas, transformedareas, overlapareas]=UMAPProscrutesOverlap(umap1, umap2, umapctr1, umapctr2, trialord1, trialord2)
warning('off')
[d,Z,transform] = procrustes(umapctr1, umapctr2);
% figure; plot(Z, transform.b*umapctr2*transform.T + transform.c, '.')

trialtypes = unique([trialord1 trialord2]);
Ntt = numel(trialtypes);

umap2tf = NaN(size(umap2));
for itt = 1:Ntt
    trialsoj = find( trialord2==trialtypes(itt) );
    umap2tf(trialsoj,:)= transform.b*umap2(trialsoj,:)*transform.T + transform.c(itt,:);
end


templateareas = NaN(Ntt,1);
transformedareas = NaN(1,Ntt);
polys = struct();
for itt = 1:Ntt
    for ses2plt = 1:2
        switch ses2plt
            case 1
                trialsoi = find( trialord1==trialtypes(itt) );
                P = umap1(trialsoi,:);
            case 2
                trialsoi = find( trialord2==trialtypes(itt) );
                P = umap2tf(trialsoi,:);
        end
        distances = mahal(P, P);
        % Define a threshold for outlier detection (e.g., 99% confidence level)
        threshold = chi2inv(0.95, size(P, 2)); % Chi-square critical value
        outliers = distances > threshold;
        % Remove outliers
        P_cleaned = P(~outliers, :);
        [k,v]=boundary(P_cleaned);
        poly = polyshape( P_cleaned(k,1), P_cleaned(k,2) );
        if abs(area(poly)-v)>10^-10
            error('boundary and polyshape outputs mismatch: please check')
        end
        
        if isempty(fieldnames(polys))
            polys = poly;
        else
            polys(ses2plt,itt) = poly;
        end
        switch ses2plt
            case 1
                templateareas(itt) = v;
            case 2
                transformedareas(itt) = v;
        end
    end
end

overlapareas = NaN(Ntt, Ntt);
for itt = 1:Ntt
    for jtt = 1:Ntt
        overlapRegion = intersect(polys(1,itt), polys(2,jtt));
        overlapArea = area(overlapRegion);
        overlapareas(itt,jtt) = overlapArea;
    end
end

end

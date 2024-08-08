addpath('C:\Users\USER\GitHub\helperfunctions')

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

hireptt = [0	101	105	106	107	109	110	111	1105	1109	1201	1299];
Ntt = numel(hireptt);
spkcntIChiV1agg = cell(1, Nsessions);
spkcntIChivisctxagg = cell(1, Nsessions);
unitlocvisctx = cell(1, Nsessions);

for ises = 1:Nsessions
pathpp = [datadir 'postprocessed\' nwbsessions{ises} '\']; % session with highest cross validation accuracy
load([pathpp, 'postprocessed.mat'])
load([pathpp, 'qc_units.mat'])

whichblock = 'ICwcfg1_presentations';
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
blocktrialorder = ICtrialtypes(vis.(whichblock).trialorder+1);
spkcnt= Rall.(whichblock)*0.4;

[v,c]=uniquecnt(blocktrialorder);
if ~isequal(hireptt, v(c==400))
    error(['check ' nwbsessions{ises}])
end

neuV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
neuvis = contains(neuallloc, 'VIS');
neuRS = unit_wfdur>0.4;
neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);

neuV1RS = neuV1 & neuRS;
neuvisRS = neuvis & neuRS;

    spkcntttcvisctx = cell(Ntt,1);
    spkcntttcV1 = cell(Ntt,1);
    for cnt = 1:Ntt
        trialsoi = blocktrialorder==hireptt(cnt);
        spkcntttcvisctx{cnt} = spkcnt(trialsoi,neuvisRS);
        spkcntttcV1{cnt} = spkcnt(trialsoi,neuV1RS);
    end

    spkcntIChiV1agg{ises} = spkcntttcV1;
    spkcntIChivisctxagg{ises} = spkcntttcvisctx;
    unitlocvisctx{ises} = neuallloc(neuvisRS);

end
save('G:\My Drive\RESEARCH\logmean_logvar\OpenScope_spkcnt_ICwcfg1.mat', 'nwbsessions', 'hireptt', 'spkcntIChiV1agg', 'spkcntIChivisctxagg', 'unitlocvisctx')

%%
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
else
    drivepath = 'G:/My Drive/';
end
load([drivepath 'RESEARCH/logmean_logvar/OpenScope_spkcnt_ICwcfg1.mat'])

neuopt = 'V1';

switch neuopt
    case 'V1'
        spkcntICttagg = spkcntIChiV1agg;
    case 'visctx'
        spkcntICttagg = spkcntIChivisctxagg;
end

Nsessions = numel(spkcntICttagg);
Ntt = numel(spkcntICttagg{1});
spkcntICttacc = cat(2, spkcntICttagg{:})';

csz = cellfun(@size, spkcntICttacc(:,1), 'UniformOutput', false);
csz = cat(1,csz{:,1});
Nneurons = sum(csz(:,2));

spkcntICttaccavg = cellfun(@mean, spkcntICttacc, 'uniformoutput', false);
spkcntICttaccvar = cellfun(@var, spkcntICttacc, 'uniformoutput', false);

spkcntICttavg = NaN(Nneurons, Ntt);
spkcntICttvar = NaN(Nneurons, Ntt);
for ii = 1:Ntt
    spkcntICttavg(:,ii) = cat(2,spkcntICttaccavg{:,ii});
    spkcntICttvar(:,ii) = cat(2,spkcntICttaccvar{:,ii});
end

%% fit log mean log var relationship for every session, every image
Ntt = numel(spkcntICttagg{1});
Nsessions = numel(spkcntICttagg);
lmlvslope = NaN(Nsessions, Ntt);
lmlvyintercept = NaN(Nsessions, Ntt);
tic
for ises = 1:Nsessions
    try
        tempspkcnt = cat(3,spkcntICttagg{ises}{:});
        Nrep = size(tempspkcnt,1);
    catch
        % trial repetitions were not the same across natural scene trial types
        csz = cellfun(@size, spkcntICttagg{ises}, 'UniformOutput', false);
        csz = cat(1,csz{:});
        Nrep = min(csz(:,1));
        if all(csz(:,2)==csz(1,2))
            Nneu = csz(1,2);
        else
            error('check number of neurons in session %d', ises)
        end
        tempspkcnt = NaN(Nrep, Nneu, Ntt);
        for n = 1:Ntt
            tempspkcnt(:,:,n) = spkcntICttagg{ises}{n}(1:Nrep,:);
        end
    end
    
    % fit log(mean) vs log(var)
    spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu
    spkcntmu = squeeze(mean(spkcntICtt,1)); % Nimg X Nneurons
    spkcntvar = squeeze(var(spkcntICtt,0,1)); % Nimg X Nneurons
    temp = spkcntvar; temp(spkcntvar==0)=NaN;
    totvar = nanmean(temp,2);
    
    for typi = 1:Ntt
        tempx = log10(spkcntmu(typi,:));
        tempy = log10(spkcntvar(typi,:));
        valneu = spkcntmu(typi,:)>0;
        LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
        b=LM.Coefficients.Estimate(1);
        a=LM.Coefficients.Estimate(2);
        lmlvslope(ises, typi) = a;
        lmlvyintercept(ises, typi) = b;
    end
end
toc

%% original data: representation similarity across sessions
repsimpairtrials = zeros(Nsessions, 1);
repsimmeantrial = zeros(Nsessions, 1);
repsimallmeantrial = zeros(Nsessions, 1);

Ntt = numel(spkcntICttagg{1});
Nsessions = numel(spkcntICttagg);

cossimperses = NaN(Ntt, Ntt, Nsessions); % 14 seconds
cossimbsperses = NaN(Ntt-1, Ntt-1, Nsessions);
for ises = 1:Nsessions
    try
        tempspkcnt = cat(3,spkcntICttagg{ises}{:});
        Nrep = size(tempspkcnt,1);
        Nneu = size(tempspkcnt,2);
    catch
        % trial repetitions were not the same across natural scene trial types
        csz = cellfun(@size, spkcntICttagg{ises}, 'UniformOutput', false);
        csz = cat(1,csz{:});
        Nrep = min(csz(:,1));
        if all(csz(:,2)==csz(1,2))
            Nneu = csz(1,2);
        else
            error('check number of neurons in session %d', ises)
        end
        tempspkcnt = NaN(Nrep, Nneu, Ntt);
        for n = 1:Ntt
            tempspkcnt(:,:,n) = spkcntICttagg{ises}{n}(1:Nrep,:);
        end
    end

    spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu
    spkcntres = spkcntICtt - mean(spkcntICtt,1); % Nrep * Nnstt * Nneu
    spkcntmu = mean(spkcntICtt,1); % 1XNimg X Nneurons
    spkcntvar = var(spkcntICtt,0,1); % 1XNimg X Nneurons

    % skipped shuffling repetition order per trial type

    % cosine similarity base subtracted
    bsspkcntICtt = spkcntICtt(:,2:end,:) - mean(spkcntICtt(:,1,:), 1);

    % cosine similarity
    C = nchoosek((1:Nrep)',2);
    Ntrialpair = size(C,1);

    newmu = squeeze(mean(spkcntICtt,1));
    mucs = normr(newmu)*normr(newmu)';
    triuind = find(triu(true(Ntt),1));
    temprepsim = NaN(Ntrialpair,1);
    temprepsimmean = NaN(Ntrialpair,1);
    temprepsimallmean = NaN(Ntrialpair,1);

    tempcossimsum = zeros(Ntt, Ntt);
    tempcossimbssum = zeros(Ntt-1, Ntt-1);
    %tic
    for c = 1:Ntrialpair
        temp1 = squeeze(spkcntICtt(C(c,1),:,:)); % Nttns*Nneurons
        temp2 = squeeze(spkcntICtt(C(c,2),:,:)); % Nttns*Nneurons

        tempbs1 = squeeze(bsspkcntICtt(C(c,1),:,:));
        tempbs2 = squeeze(bsspkcntICtt(C(c,2),:,:));

        % normalized_M = normr(M) takes a single matrix or cell array of matrices, M, and returns the matrices with rows normalized to a length of one.
        % tempdot = temp1*temp2';
        % tempnorm = sqrt(diag(temp1*temp1')) * sqrt(diag(temp2*temp2'))';
        %figure; plot(normr(temp1)*normr(temp2)', tempdot./tempnorm, '.') % essentially equal
        tempcs = normr(temp1)*normr(temp2)';
        tempcossimsum = tempcossimsum + tempcs;

        temprepsimmean(c) = corr( mucs(triuind), tempcs(triuind) );
        temprepsimallmean(c) = corr( mucs(:), tempcs(:) );

        tempcs1 = normr(temp1)*normr(temp1)';
        tempcs2 = normr(temp2)*normr(temp2)';
        temprepsim(c) = corr( tempcs1(triuind), tempcs2(triuind) );

        tempcossimbssum = tempcossimbssum + normr(tempbs1)*normr(tempbs2)';
    end
    %toc
    cossimperses(:,:, ises) = tempcossimsum/Ntrialpair;
    cossimbsperses(:,:, ises) = tempcossimbssum/Ntrialpair;

    repsimpairtrials(ises) = mean(temprepsim);
    repsimmeantrial(ises) = mean(temprepsimmean);
    repsimallmeantrial(ises) = mean(temprepsimallmean);
end


cossimvecses = reshape(cossimperses,Ntt^2,Nsessions);
repsimmatall = corr(cossimvecses);
cossimtriuses = cossimvecses(find(triu(true(Ntt),1)),:);
repsimmat = corr(cossimtriuses);
repsimvec = repsimmat(triu(true(Nsessions),1));

cossimdiagses = mean( cossimvecses(find(eye(Ntt)),:), 1);
if numel(cossimdiagses)~=Nsessions
    error('check cossimdiagses')
end

fprintf('original data diag mean %.4f rep sim mean %.4f, median %.4f (Q1: %.4f, Q3: %.4f)\n', ...
    mean(cossimdiagses), mean(repsimvec), median(repsimvec), prctile(repsimvec,25), prctile(repsimvec,75))

save([drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_representationsimilarity_', neuopt, '.mat'], ...
    'lmlvslope', 'lmlvyintercept', 'cossimperses', 'cossimbsperses', ...
    'repsimpairtrials', 'repsimmeantrial', 'repsimallmeantrial')

%% vary dispersion and fano factor (keep total variance per image same): representation similarity across sessions
disperses = -1:0.1:2;
alltrialpairs = false;
Nshuffle = 1000; % 0 if we don't want to shuffle, 1000 otherwise (10000 takes 1400s per dispersion, roughly 30 min)

if Nshuffle>0
    cossimfn = [drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_keeptotalvariance_', neuopt, '_shuffle.mat'];
else
    cossimfn = [drivepath 'RESEARCH/logmean_logvar/OpenScopeIC_keeptotalvariance_', neuopt, '.mat'];
end

fffitperses = zeros(Ntt, length(disperses), Nsessions);
cossimtvperses = cell(size(disperses));
cossimbstvperses = cell(size(disperses));

repsimtvpairtrials = zeros(Nsessions, numel(disperses));
repsimtvmeantrial = zeros(Nsessions, numel(disperses));
repsimalltvmeantrial = zeros(Nsessions, numel(disperses));

Ntt = numel(spkcntICttagg{1});
Nsessions = numel(spkcntICttagg);

PTRE1asIC1vsLC1 = zeros(Nsessions, numel(disperses));
PTRE2asIC2vsLC2 = zeros(Nsessions, numel(disperses));
PIC1asXRE1vsXRE2 = zeros(Nsessions, numel(disperses));
PIC2asXRE2vsXRE1 = zeros(Nsessions, numel(disperses));


% 45sec per iteration
for ii = 1:numel(disperses)
    
    dc = tic;
    tempcossim = NaN(Ntt, Ntt, Nsessions); 
    tempcossimbs = NaN(Ntt-1, Ntt-1, Nsessions); 
    for ises = 1:Nsessions
        try
            tempspkcnt = cat(3,spkcntICttagg{ises}{:});
            Nrep = size(tempspkcnt,1);
            Nneu = size(tempspkcnt,2);
        catch
            % trial repetitions were not the same across natural scene trial types
            csz = cellfun(@size, spkcntICttagg{ises}, 'UniformOutput', false);
            csz = cat(1,csz{:});
            Nrep = min(csz(:,1));
            if all(csz(:,2)==csz(1,2))
                Nneu = csz(1,2);
            else
                error('check number of neurons in session %d', ises)
            end
            tempspkcnt = NaN(Nrep, Nneu, Ntt);
            for n = 1:Ntt
                tempspkcnt(:,:,n) = spkcntICttagg{ises}{n}(1:Nrep,:);
            end
        end
        
        % % destroy noise correlations: shuffle trial order for every neuron and every trial type
        % if Nshuffle>0
        %     tempspkcnt0 = tempspkcnt;
        %     for ci = 1:Nneu
        %         for n = 1:Nhireptt
        %             tempvec = squeeze(tempspkcnt0(:,ci,n));
        %             tempspkcnt(:,ci,n) = tempvec(randperm(numel(tempvec)));
        %         end
        %     end
        % end
        
        % fit log(mean) vs log(var)
        spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu
        spkcntres = spkcntICtt - mean(spkcntICtt,1); % Nrep * Nnstt * Nneu
        spkcntmu = mean(spkcntICtt,1); % 1XNimg X Nneurons
        spkcntvar = var(spkcntICtt,0,1); % 1XNimg X Nneurons
        temp = spkcntvar; temp(spkcntvar==0)=NaN;
        totvar = nanmean(temp,2);
        
        tempx = log10(spkcntmu);
        tempx(spkcntmu==0) = NaN;
        meanx = squeeze(nanmean(tempx,3)); % average across neurons: 1XNimg
        
        % when original log(mean) vs log(var) relationship is y=ax+b,
        % and new relationship is y=cx+d
        % d = (a-c)*mean(x)+b —> this keeps mean(y) constant
        % original variance: v0
        % new variance: v1
        % v1 = 10^(c/a * (log10(v0)-b) +d);
        % rescaled residual for every image and every neuron.
        % noise correlation does not change
        Avec = lmlvslope(ises,:);
        Bvec = lmlvyintercept(ises,:);
        Cvec = disperses(ii)*ones(1,Ntt);
        Dvec = (Avec-Cvec).*meanx + Bvec;
        newspkcntvar = 10.^( (Cvec./Avec).*(log10(spkcntvar)-Bvec) + Dvec);
        newspkcntres = spkcntres .* sqrt(newspkcntvar./spkcntvar);
        newspkcntICtt = mean(spkcntICtt,1)+newspkcntres;
        
        fffitperses(:,ii,ises) = Dvec;


        %{
% sanity check: fitting new newspkcntnstt
figure; hold all
plot( log10(squeeze(mean(spkcntICtt(:,typi,:),1))), log10(squeeze(var(spkcntICtt(:,typi,:),0,1))), 'k.')
plot( log10(squeeze(mean(newspkcntICtt(:,typi,:),1))), log10(squeeze(var(newspkcntICtt(:,typi,:),0,1))), 'r.')
xl = xlim; yl = ylim;
tempx = log10(squeeze(mean(spkcntICtt(:,typi,:),1)));
tempy = log10(squeeze(var(spkcntICtt(:,typi,:),0,1)));
valneu = squeeze(mean(spkcntICtt(:,typi,:),1))>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
plot(xl, a*xl+b, 'k-')
text(xl(2), yl(1)+0.1*range(yl), sprintf('original data: y=%.2fx+%.2f, fitted: y=%.2fx+%.2f', lmlvslope(ises,typi), lmlvyintercept(ises,typi), a, b), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom')
tempx = log10(squeeze(mean(newspkcntICtt(:,typi,:),1)));
tempy = log10(squeeze(var(newspkcntICtt(:,typi,:),0,1)));
valneu = squeeze(mean(newspkcntICtt(:,typi,:),1))>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
d=LM.Coefficients.Estimate(1);
c=LM.Coefficients.Estimate(2);
plot(xl, c*xl+d, 'r-')
text(xl(2), yl(1), sprintf('intended rescaling: y=%.2fx+%.2f, fitted: y=%.2fx+%.2f', Cvec(typi), Dvec(typi), c,d), 'HorizontalAlignment','right', 'VerticalAlignment', 'bottom')
axis([xl yl])
xlabel('log_1_0(mean)')
ylabel('log_1_0(var)')

% noise correlation remains constant 
figure; plot( corr(squeeze(spkcntICtt(:,typi,:))), corr(squeeze(newspkcntICtt(:,typi,:))), '.' )
figure; histogram(corr(squeeze(spkcntICtt(:,typi,:))))
   
% check that total variance remains constant, i.e., original mean(y) and new mean(y) are (approx.) equal across images
varorig = var(spkcntICtt,0,1);
varorig(spkcntmu==0) = NaN;
varnew = var(newspkcntICtt,0,1);
varnew(spkcntmu==0) = NaN;
figure; plot( squeeze(nanmean(log10(varorig) , 3)), ...
    squeeze(nanmean( log10(varnew), 3)), '.' );
hold on; xl = xlim; plot(xl,xl,'r-')
        %}
        


        % shuffle the trial order for each trial type
        newspkcntICttshuf = newspkcntICtt; % Nrep*Nttns*Nneurons
        for cnt = 1:Ntt
            temprepord = randperm(Nrep);
            tempmat = squeeze(newspkcntICtt(:,cnt,:));
            newspkcntICttshuf(:,cnt,:) = tempmat(temprepord,:);
        end
        if max(abs( squeeze(mean(tempspkcnt,1))' - squeeze(mean(newspkcntICttshuf,1))), [],'all')>2^-30
            error('check spkcntnsttshuf')
        end

        % cosine similarity base subtracted
        newbsspkcntICttshuf = newspkcntICttshuf(:,2:end,:) - mean(newspkcntICttshuf(:,1,:), 1);

        % cosine similarity
        if Nshuffle==0
            C = nchoosek((1:Nrep)',2);
            if alltrialpairs
                Ntrialpair = size(C,1);
            else
                Ntrialpair = min([size(C,1) 10000]);
            end
        else
            Ntrialpair = Nshuffle;
        end

        newmu = squeeze(mean(newspkcntICttshuf,1));
        mucs = normr(newmu)*normr(newmu)';
        triuind = find(triu(true(Ntt),1));
        temprepsim = NaN(Ntrialpair,1);
        temprepsimmean = NaN(Ntrialpair,1);
        temprepsimallmean = NaN(Ntrialpair,1);

        tempNTRE1asIC1vsLC1 = 0;
        tempNTRE2asIC2vsLC2 = 0;
        tempNIC1asXRE1vsXRE2 = 0;
        tempNIC2asXRE2vsXRE1 = 0;
        %tic
        tempcossimsum = zeros(Ntt, Ntt);
        tempcossimbssum = zeros(Ntt-1, Ntt-1);
        for c = 1:Ntrialpair
            if Nshuffle==0
                temp1 = squeeze(newspkcntICttshuf(C(c,1),:,:)); % Nttns*Nneurons
                temp2 = squeeze(newspkcntICttshuf(C(c,2),:,:)); % Nttns*Nneurons

            tempbs1 = squeeze(newbsspkcntICttshuf(C(c,1),:,:));
            tempbs2 = squeeze(newbsspkcntICttshuf(C(c,2),:,:));
            else
                % since newspkcntICttshuf is already shuffled along the
                % trial axis, we only need to shuffle along the neuron axis
                temp12 = newspkcntICttshuf;
                tempbs12 = newbsspkcntICttshuf;
                for ci = 1:Nneu
                    temprepord = randperm(Nrep);
                    tempmat = newspkcntICttshuf(:,:,ci);
                    temp12(:,:,ci) = tempmat(temprepord,:);

                    tempbsmat = newbsspkcntICttshuf(:,:,ci);
                    tempbs12(:,:,ci) = tempbsmat(temprepord,:);
                end
                temp1 = squeeze(temp12(1,:,:));
                temp2 = squeeze(temp12(2,:,:));

                tempbs1 = squeeze(tempbs12(1,:,:));
                tempbs2 = squeeze(tempbs12(2,:,:));
            end
            % normalized_M = normr(M) takes a single matrix or cell array of matrices, M, and returns the matrices with rows normalized to a length of one.
            % tempdot = temp1*temp2';
            % tempnorm = sqrt(diag(temp1*temp1')) * sqrt(diag(temp2*temp2'))';
            %figure; plot(normr(temp1)*normr(temp2)', tempdot./tempnorm, '.') % essentially equal
            tempcs = normr(temp1)*normr(temp2)';
            tempcossimsum = tempcossimsum + tempcs;
            
            temprepsimmean(c) = corr( mucs(triuind), tempcs(triuind) );
            temprepsimallmean(c) = corr( mucs(:), tempcs(:) );
            
            tempcs1 = normr(temp1)*normr(temp1)';
            tempcs2 = normr(temp2)*normr(temp2)';
            temprepsim(c) = corr( tempcs1(triuind), tempcs2(triuind) );

        tempNTRE1asIC1vsLC1 = tempNTRE1asIC1vsLC1 + nnz(tempcs(hireptt==1105, hireptt==106)>tempcs(hireptt==1105, hireptt==107));
        tempNTRE2asIC2vsLC2 = tempNTRE2asIC2vsLC2 + nnz(tempcs(hireptt==1109, hireptt==111)>tempcs(hireptt==1109, hireptt==110));
        tempNIC1asXRE1vsXRE2 = tempNIC1asXRE1vsXRE2 + nnz(tempcs(hireptt==106, hireptt==1201)>tempcs(hireptt==106, hireptt==1299));
        tempNIC2asXRE2vsXRE1 = tempNIC2asXRE2vsXRE1 + nnz(tempcs(hireptt==111, hireptt==1201)<tempcs(hireptt==111, hireptt==1299));

            tempcossimbssum = tempcossimbssum + normr(tempbs1)*normr(tempbs2)';
        end
        %toc
        tempcossim(:,:, ises) = tempcossimsum/Ntrialpair;
        tempcossimbs(:,:, ises) = tempcossimbssum/Ntrialpair;
        
        repsimtvpairtrials(ises,ii) = mean(temprepsim);
        repsimtvmeantrial(ises,ii) = mean(temprepsimmean);
        repsimalltvmeantrial(ises,ii) = mean(temprepsimallmean);

        PTRE1asIC1vsLC1(ises,ii) = tempNTRE1asIC1vsLC1/Ntrialpair;
        PTRE2asIC2vsLC2(ises,ii) = tempNTRE2asIC2vsLC2/Ntrialpair;
        PIC1asXRE1vsXRE2(ises,ii) = tempNIC1asXRE1vsXRE2/Ntrialpair;
        PIC2asXRE2vsXRE1(ises,ii) = tempNIC2asXRE2vsXRE1/Ntrialpair;
                
    end
    
    cossimtvperses{ii} = tempcossim;
    cossimbstvperses{ii} = tempcossimbs;
    
    
    toc(dc)
    
    
    cossimvecses = reshape(tempcossim,Ntt^2,Nsessions);
    repsimmatall = corr(cossimvecses);
    cossimtriuses = cossimvecses(find(triu(true(Ntt),1)),:);
    repsimmat = corr(cossimtriuses);
    repsimvec = repsimmat(triu(true(Nsessions),1));
    
    cossimdiagses = mean( cossimvecses(find(eye(Ntt)),:), 1);
    if numel(cossimdiagses)~=Nsessions
        error('check cossimdiagses')
    end

    fprintf('slope %.2f diag mean %.4f rep sim mean %.4f, median %.4f (Q1: %.4f, Q3: %.4f)\n', ...
        disperses(ii), mean(cossimdiagses), mean(repsimvec), median(repsimvec), prctile(repsimvec,25), prctile(repsimvec,75))
    
end
save(cossimfn, ...
    'lmlvslope', 'lmlvyintercept', 'fffitperses', 'disperses', 'cossimtvperses', 'cossimbstvperses', ...
    'PTRE1asIC1vsLC1', 'PTRE2asIC2vsLC2', 'PIC1asXRE1vsXRE2', 'PIC2asXRE2vsXRE1', ...
    'repsimtvpairtrials', 'repsimtvmeantrial', 'repsimalltvmeantrial')

%% cosine similarity and representational similarity
incldiagrepsim = true;
cossimavgagg = normr(spkcntICttavg')*normr(spkcntICttavg')';
cossimavgaggvec = reshape(cossimavgagg,Ntt^2,1);
cossimavgaggtriu = cossimavgaggvec(find(triu(true(Ntt),1)));


% original data
tempcossim = cossimperses ;

cossimvecses = reshape(tempcossim,Ntt^2,Nsessions);
repsimmatall = corr(cossimvecses);
cossimtriuses = cossimvecses(find(triu(true(Ntt),1)),:);
repsimmat = corr(cossimtriuses);
if incldiagrepsim
    repsimvec = repsimmatall(triu(true(Nsessions),1));
else
    repsimvec = repsimmat(triu(true(Nsessions),1));
end

repsimmean = mean(repsimvec);
repsimmedian = median(repsimvec);
repsimsem = std(repsimvec)/sqrt(numel(repsimvec));
repsimq1 = prctile(repsimvec,25);
repsimq3 = prctile(repsimvec,75);

cossimdiagses = mean( cossimvecses(find(eye(Ntt)),:), 1);
if numel(cossimdiagses)~=Nsessions
    error('check cossimdiagses')
end

diagses = cossimdiagses;
diagsesmean = mean(cossimdiagses);
diagsessem = std(cossimdiagses)/sqrt(Nsessions);

if incldiagrepsim
    repsimperses = corr(cossimvecses, cossimavgaggvec);
else
    repsimperses = corr(cossimtriuses, cossimavgaggtriu);
end
repsimses = repsimperses;
repsimsesmean = mean(repsimperses);
repsimsessem = std(repsimperses)/sqrt(Nsessions);



% across disperses
diagtvses = NaN(Nsessions, numel(disperses));
diagtvsesmean = NaN(size(disperses));
diagtvsessem = NaN(size(disperses));

repsimtvses = NaN(Nsessions, numel(disperses));
repsimtvsesmean = NaN(size(disperses));
repsimtvsessem = NaN(size(disperses));

% mean and sem across pairs of sessions
repsimtvmean = NaN(size(disperses));
repsimtvsem = NaN(size(disperses));
repsimtvmedian = NaN(size(disperses));
repsimtvq1 = NaN(size(disperses));
repsimtvq3 = NaN(size(disperses));

for ii = 1:numel(disperses)
    
    tempcossim = cossimtvperses{ii} ;
    
    cossimvecses = reshape(tempcossim,Ntt^2,Nsessions);
    repsimmatall = corr(cossimvecses);
    cossimtriuses = cossimvecses(find(triu(true(Ntt),1)),:);
    repsimmat = corr(cossimtriuses);
    if incldiagrepsim
        repsimvec = repsimmatall(triu(true(Nsessions),1));
    else
        repsimvec = repsimmat(triu(true(Nsessions),1));
    end
    
    repsimtvmean(ii) = mean(repsimvec);
    repsimtvmedian(ii) = median(repsimvec);
    repsimtvsem(ii) = std(repsimvec)/sqrt(numel(repsimvec));
    repsimtvq1(ii) = prctile(repsimvec,25);
    repsimtvq3(ii) = prctile(repsimvec,75);
    
    cossimdiagses = mean( cossimvecses(find(eye(Ntt)),:), 1);
    if numel(cossimdiagses)~=Nsessions
        error('check cossimdiagses')
    end
    
    diagtvses(:,ii) = cossimdiagses;
    diagtvsesmean(ii) = mean(cossimdiagses);
    diagtvsessem(ii) = std(cossimdiagses)/sqrt(Nsessions);
    
    if incldiagrepsim
        repsimperses = corr(cossimvecses, cossimavgaggvec);
    else
        repsimperses = corr(cossimtriuses, cossimavgaggtriu);
    end
    repsimtvses(:,ii) = repsimperses;
    repsimtvsesmean(ii) = mean(repsimperses);
    repsimtvsessem(ii) = std(repsimperses)/sqrt(Nsessions);
    
end

%% cosine similarity and representational similarity of base-subtracted activity pattern
spkcntnsttbs = spkcntICttavg(:,2:end)-spkcntICttavg(:,1); % Nneurons*Nnstt
cossimbsavgagg = normr(spkcntnsttbs')*normr(spkcntnsttbs')';
cossimbsaggvec = reshape(cossimbsavgagg, (Ntt-1)^2,1);
cossimbsaggtriu = cossimbsaggvec(find(triu(true(Ntt-1),1)));

% original data
tempcossimbs = cossimbsperses ;

cossimbsvecses = reshape(tempcossimbs,(Ntt-1)^2,Nsessions);
repsimbsmatall = corr(cossimbsvecses);
cossimbstriuses = cossimbsvecses(find(triu(true(Ntt-1),1)),:);
repsimbsmat = corr(cossimbstriuses);
if incldiagrepsim
    repsimbsvec = repsimbsmatall(triu(true(Nsessions),1));
else
    repsimbsvec = repsimbsmat(triu(true(Nsessions),1));
end

repsimbsmean = mean(repsimbsvec);
repsimbsmedian = median(repsimbsvec);
repsimbssem = std(repsimbsvec)/sqrt(numel(repsimbsvec));
repsimbsq1 = prctile(repsimbsvec,25);
repsimbsq3 = prctile(repsimbsvec,75);

cossimbsdiagses = mean( cossimbsvecses(find(eye(Ntt-1)),:), 1);
if numel(cossimbsdiagses)~=Nsessions
    error('check cossimdiagses')
end

diagbsses = cossimbsdiagses;
diagbssesmean = mean(cossimbsdiagses);
diagbssessem = std(cossimbsdiagses)/sqrt(Nsessions);

if incldiagrepsim
    repsimbsperses = corr(cossimbsvecses, cossimbsaggvec);
else
    repsimbsperses = corr(cossimbstriuses, cossimbsaggtriu);
end
repsimbsses = repsimbsperses;
repsimbssesmean = mean(repsimbsperses);
repsimbssessem = std(repsimbsperses)/sqrt(Nsessions);


% across disperses
diagbstvses = NaN(Nsessions, numel(disperses));
diagbstvsesmean = NaN(size(disperses));
diagbstvsessem = NaN(size(disperses));

repsimbstvses = NaN(Nsessions, numel(disperses));
repsimbstvsesmean = NaN(size(disperses));
repsimbstvsessem = NaN(size(disperses));

% mean and sem across pairs of sessions
repsimbstvmean = NaN(size(disperses));
repsimbstvsem = NaN(size(disperses));
repsimbstvmedian = NaN(size(disperses));
repsimbstvq1 = NaN(size(disperses));
repsimbstvq3 = NaN(size(disperses));

for ii = 1:numel(disperses)
    
    tempcossimbs = cossimbstvperses{ii} ;
    
    cossimbsvecses = reshape(tempcossimbs,(Ntt-1)^2,Nsessions);
    repsimbsmatall = corr(cossimbsvecses);
    cossimbstriuses = cossimbsvecses(find(triu(true(Ntt-1),1)),:);
    repsimbsmat = corr(cossimbstriuses);
    if incldiagrepsim
        repsimbsvec = repsimbsmatall(triu(true(Nsessions),1));
    else
        repsimbsvec = repsimbsmat(triu(true(Nsessions),1));
    end
    
    repsimbstvmean(ii) = mean(repsimbsvec);
    repsimbstvmedian(ii) = median(repsimbsvec);
    repsimbstvsem(ii) = std(repsimbsvec)/sqrt(numel(repsimbsvec));
    repsimbstvq1(ii) = prctile(repsimbsvec,25);
    repsimbstvq3(ii) = prctile(repsimbsvec,75);
    
    cossimbsdiagses = mean( cossimbsvecses(find(eye(Ntt-1)),:), 1);
    if numel(cossimbsdiagses)~=Nsessions
        error('check cossimdiagses')
    end
    
    diagbstvses(:,ii) = cossimbsdiagses;
    diagbstvsesmean(ii) = mean(cossimbsdiagses);
    diagbstvsessem(ii) = std(cossimbsdiagses)/sqrt(Nsessions);
    
    if incldiagrepsim
        repsimbsperses = corr(cossimbsvecses, cossimbsaggvec);
    else
        repsimbsperses = corr(cossimbstriuses, cossimbsaggtriu);
    end
    repsimbstvses(:,ii) = repsimbsperses;
    repsimbstvsesmean(ii) = mean(repsimbsperses);
    repsimbstvsessem(ii) = std(repsimbsperses)/sqrt(Nsessions);
    
end

%%
if Nshuffle>0
anot = sprintf('Open Scope ICwcfg1 high-repetition trials SHUFFLED: representation similarity');
else
anot = sprintf('Open Scope ICwcfg1 high-repetition trials: representation similarity');
end
xl = [disperses(1) disperses(end)];
yl= [0 1];
figure
annotation('textbox', [0.1 0.91 0.9 0.1], 'String', anot, 'edgecolor', 'none');%,'FontSize', 12)
for isp = 1:3
    switch isp
        case 1
            temprepsim = repsimtvpairtrials;
            temporig = repsimpairtrials;
            sptitle = 'trial pairs, off-diagonal';
        case 2
            temprepsim = repsimtvmeantrial;
            temporig = repsimmeantrial;
            sptitle = 'mean vs each trial, off-diagonal';
        case 3
            temprepsim = repsimalltvmeantrial;
            temporig = repsimallmeantrial;
            sptitle = 'mean vs each trial, all matrix elements';
    end
subplot(2,2,isp)
hold all
plot(disperses, temprepsim, '-', 'Color', 0.5*[1 1 1])
errorbar(disperses, mean(temprepsim,1), std(temprepsim,0,1)/sqrt(Nsessions), 'bx', 'LineWidth', 1)
[mv,mi]=max(mean(temprepsim,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot(xl, mean(temporig)*[1 1], 'g-', 'linewidth', 1)
text(xl(1), mean(temporig), sprintf('original mean=%.4f', mean(temporig)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
ylim(yl)
xlim(xl)
xlabel('log(mean) vs log(var) slope')
ylabel('representation similarity')
title(sptitle)
end

%%
fs=10;
xl = [disperses(1) disperses(end)];
if Nshuffle>0
anot = sprintf('Open Scope ICwcfg1 high-repetition trials SHUFFLED: %d %s neurons across n=%d sessions', size(spkcntICttavg,1), neuopt, Nsessions);
else
anot = sprintf('Open Scope ICwcfg1 high-repetition trials: %d %s neurons across n=%d sessions', size(spkcntICttavg,1), neuopt, Nsessions);
end
figure('Position', [0 0 1000 750])
annotation('textbox', [0.1 0.9 0.8 0.1], 'String', anot, 'edgecolor', 'none','FontSize', fs)
subplot(2,3,1)
hold all
plot(disperses, diagtvses, '-', 'Color', 0.5*[1 1 1])
errorbar(disperses, diagtvsesmean, diagtvsessem, 'bx', 'LineWidth', 1)
yl= ylim;
[mv,mi]=max(mean(diagtvsesmean,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot([xl(1) disperses(mi)], mv*[1 1], 'r--')
text(xl(1), mv, num2str(mv), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot(xl, diagsesmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), diagsesmean, sprintf('original mean=%.4f', diagsesmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
ylim(yl)
xlim(xl)
set(gca, 'FontSize', fs)
title('diagonal across sessions')
xlabel('log(mean) vs log(var) slope')
ylabel('cosine similarity')

subplot(2,3,2)
hold all
plot(disperses, repsimtvses, '-', 'Color', 0.5*[1 1 1])
errorbar(disperses, repsimtvsesmean, repsimtvsessem, 'bx', 'LineWidth', 1)
yl= ylim;
[mv,mi]=max(mean(repsimtvsesmean,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot([xl(1) disperses(mi)], mv*[1 1], 'r--')
text(xl(1), mv, num2str(mv), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot(xl, repsimsesmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), repsimsesmean, sprintf('original mean=%.4f', repsimsesmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
ylim(yl)
xlim(xl)
set(gca, 'FontSize', fs)
title('across sessions (against agg mean vec)')
xlabel('log(mean) vs log(var) slope')
ylabel('representation similarity')

subplot(2,3,3)
hold all
errorbar(disperses, repsimtvmean, repsimtvsem, 'bx-', 'LineWidth', 1)
errorbar(disperses, repsimtvmedian, (repsimtvmedian-repsimtvq1), (repsimtvq3-repsimtvmedian), 'ko')
yl= ylim;
[mv,mi]=max(mean(repsimtvmean,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot([xl(1) disperses(mi)], mv*[1 1], 'r--')
text(xl(1), mv, num2str(mv), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
plot(xl, repsimmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), repsimmean, sprintf('original mean=%.4f', repsimmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
ylim(yl)
xlim(xl)
set(gca, 'FontSize', fs)
legend({'mean+/-SEM', 'median+/-IQI'}, 'location', 'best')
title('across pairs of sessions')
xlabel('log(mean) vs log(var) slope')
ylabel('representation similarity')

subplot(2,3,4)
hold all
plot(disperses, diagbstvses, '-', 'Color', 0.5*[1 1 1])
errorbar(disperses, diagbstvsesmean, diagbstvsessem, 'bx', 'LineWidth', 1)
plot(xl, diagbssesmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), diagbssesmean, sprintf('original mean=%.4f', diagbssesmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
xlim(xl)
set(gca, 'FontSize', fs)
title('diagonal across sessions')
xlabel('log(mean) vs log(var) slope')
ylabel('base-subtracted cosine similarity')

subplot(2,3,5)
hold all
plot(disperses, repsimbstvses, '-', 'Color', 0.5*[1 1 1])
errorbar(disperses, repsimbstvsesmean, repsimbstvsessem, 'bx', 'LineWidth', 1)
plot(xl, repsimbssesmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), repsimbssesmean, sprintf('original mean=%.4f', repsimbssesmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
xlim(xl)
set(gca, 'FontSize', fs)
title('across sessions (against agg mean vec)')
xlabel('log(mean) vs log(var) slope')
ylabel('base-subtracted representation similarity')

subplot(2,3,6)
hold all
errorbar(disperses, repsimbstvmean, repsimbstvsem, 'bx-', 'LineWidth', 1)
errorbar(disperses, repsimbstvmedian, (repsimbstvmedian-repsimbstvq1), (repsimbstvq3-repsimbstvmedian), 'ko')
plot(xl, repsimbsmean*[1 1], 'g-', 'linewidth', 1)
text(xl(1), repsimbsmean, sprintf('original mean=%.4f', repsimbsmean), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'g', 'FontWeight', 'bold')
xlim(xl)
set(gca, 'FontSize', fs)
legend({'mean+/-SEM', 'median+/-IQI'}, 'location', 'best')
title('across pairs of sessions')
xlabel('log(mean) vs log(var) slope')
ylabel('base-subtracted representation similarity')


%% trial-by-trial inference cosine similarity 
pltsesavg = true;
PTREasICvsLC = mean(cat(3,PTRE1asIC1vsLC1, PTRE2asIC2vsLC2),3);
PICasXRE = mean(cat(3,PIC1asXRE1vsXRE2, PIC2asXRE2vsXRE1),3);

xl = [disperses(1) disperses(end)];
yl= [0 1];
if Nshuffle>0
anot = sprintf('Open Scope ICwcfg1 high-repetition trials SHUFFLED: trial-by-trial cosine similarity', size(spkcntICttavg,1), neuopt, Nsessions);
else
anot = sprintf('Open Scope ICwcfg1 high-repetition trials: trial-by-trial cosine similarity', size(spkcntICttavg,1), neuopt, Nsessions);
end

figure('Position', [0 0 1000 750])
annotation('textbox', [0.1 0.9 0.8 0.1], 'String', anot, 'edgecolor', 'none','FontSize', 12)
for isp = 1:6
switch isp
    case 1
        tempprct = PTRE1asIC1vsLC1;
        ylab = 'P(cossim(TRE1,IC1)>cossim(TRE1,LC1))';
    case 2
        tempprct = PTRE2asIC2vsLC2;
        ylab = 'P(cossim(TRE2,IC2)>cossim(TRE2,LC2))';
    case 3
        tempprct = PTREasICvsLC;
        ylab = 'P(cossim(TRE,IC)>cossim(TRE,LC))';
    case 4
        tempprct = PIC1asXRE1vsXRE2;
        ylab = 'P(cossim(IC1,XRE1)>cossim(IC1,XRE2))';
    case 5
        tempprct = PIC2asXRE2vsXRE1;
        ylab = 'P(cossim(IC2,XRE2)>cossim(IC2,XRE1))';
    case 6
        tempprct = PICasXRE;
        ylab = 'P(cossim(IC,XREiso)>cossim(IC,XREorth))';
end
subplot(2,3,isp)
hold all
plot(xl, [0.5 0.5], 'k--')
if pltsesavg
plot(disperses, tempprct, '-', 'Color', [0.5 0.5 0.5])
errorbar(disperses, mean(tempprct,1), std(tempprct,0,1)/sqrt(Nsessions), 'bx-', 'linewidth', 1)
[mv,mi]=max(mean(tempprct,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
else
plot(disperses, tempprct, '.-');
end
ylim(yl)
xlim(xl)
set(gca, 'XGrid', 'on', 'FontSize', fs)
xlabel('log(mean) vs log(var) slope')
ylabel(ylab)
end


%% vary dispersion and fano factor (keep total variance per image same): inference cosine similarity across sessions

cossimTRE1asIC1orLC1 = cell(Nsessions, numel(disperses));
cossimTRE2asIC2orLC2 = cell(Nsessions, numel(disperses));
cossimIC1asXRE = cell(Nsessions, numel(disperses));
cossimIC2asXRE = cell(Nsessions, numel(disperses));

Ntt = numel(spkcntICttagg{1});
Nsessions = numel(spkcntICttagg);

% 45sec per iteration
for ii = 1:numel(disperses)
    
    dc = tic;
    for ises = 1:Nsessions
        try
            tempspkcnt = cat(3,spkcntICttagg{ises}{:});
            Nrep = size(tempspkcnt,1);
        catch
            % trial repetitions were not the same across natural scene trial types
            csz = cellfun(@size, spkcntICttagg{ises}, 'UniformOutput', false);
            csz = cat(1,csz{:});
            Nrep = min(csz(:,1));
            if all(csz(:,2)==csz(1,2))
                Nneu = csz(1,2);
            else
                error('check number of neurons in session %d', ises)
            end
            tempspkcnt = NaN(Nrep, Nneu, Ntt);
            for n = 1:Ntt
                tempspkcnt(:,:,n) = spkcntICttagg{ises}{n}(1:Nrep,:);
            end
        end
        
        
        % fit log(mean) vs log(var)
        spkcntICtt = permute(tempspkcnt, [1 3 2]); % Nrep * Nnstt * Nneu
        spkcntres = spkcntICtt - mean(spkcntICtt,1); % Nrep * Nnstt * Nneu
        spkcntmu = mean(spkcntICtt,1); % 1XNimg X Nneurons
        spkcntvar = var(spkcntICtt,0,1); % 1XNimg X Nneurons
        temp = spkcntvar; temp(spkcntvar==0)=NaN;
        totvar = nanmean(temp,2);
        
        tempx = log10(spkcntmu);
        tempx(spkcntmu==0) = NaN;
        meanx = squeeze(nanmean(tempx,3)); % average across neurons: 1XNimg
        
        % when original log(mean) vs log(var) relationship is y=ax+b,
        % and new relationship is y=cx+d
        % d = (a-c)*mean(x)+b —> this keeps mean(y) constant
        % original variance: v0
        % new variance: v1
        % v1 = 10^(c/a * (log10(v0)-b) +d);
        % rescaled residual for every image and every neuron.
        % noise correlation does not change
        Avec = lmlvslope(ises,:);
        Bvec = lmlvyintercept(ises,:);
        Cvec = disperses(ii)*ones(1,Ntt);
        Dvec = (Avec-Cvec).*meanx + Bvec;
        newspkcntvar = 10.^( (Cvec./Avec).*(log10(spkcntvar)-Bvec) + Dvec);
        newspkcntres = spkcntres .* sqrt(newspkcntvar./spkcntvar);
        newspkcntICtt = mean(spkcntICtt,1)+newspkcntres;
        
%max(abs(spkcntmu-mean(newspkcntICtt,1)),[],'all')
IC1avgvec = squeeze(mean(newspkcntICtt(:,hireptt==106,:),1));
LC1avgvec = squeeze(mean(newspkcntICtt(:,hireptt==107,:),1));
LC2avgvec = squeeze(mean(newspkcntICtt(:,hireptt==110,:),1));
IC2avgvec = squeeze(mean(newspkcntICtt(:,hireptt==111,:),1));
XRE1avgvec = squeeze(mean(newspkcntICtt(:,hireptt==1201,:),1));
XRE2avgvec = squeeze(mean(newspkcntICtt(:,hireptt==1299,:),1));

cossimTRE1IC1 = normr( squeeze(newspkcntICtt(:,hireptt==1105,:)) )*normr(IC1avgvec')';
cossimTRE1LC1 = normr( squeeze(newspkcntICtt(:,hireptt==1105,:)) )*normr(LC1avgvec')';
cossimTRE1asIC1orLC1{ises, ii} = cat(2,cossimTRE1IC1, cossimTRE1LC1);

cossimTRE2IC2 = normr( squeeze(newspkcntICtt(:,hireptt==1109,:)) )*normr(IC2avgvec')';
cossimTRE2LC2 = normr( squeeze(newspkcntICtt(:,hireptt==1109,:)) )*normr(LC2avgvec')';
cossimTRE2asIC2orLC2{ises, ii} = cat(2,cossimTRE2IC2, cossimTRE2LC2);

cossimIC1XRE1 = normr( squeeze(newspkcntICtt(:,hireptt==106,:)) )*normr(XRE1avgvec')';
cossimIC1XRE2 = normr( squeeze(newspkcntICtt(:,hireptt==106,:)) )*normr(XRE2avgvec')';
cossimIC1asXRE{ises, ii} = cat(2,cossimIC1XRE1, cossimIC1XRE2);

cossimIC2XRE1 = normr( squeeze(newspkcntICtt(:,hireptt==111,:)) )*normr(XRE1avgvec')';
cossimIC2XRE2 = normr( squeeze(newspkcntICtt(:,hireptt==111,:)) )*normr(XRE2avgvec')';
cossimIC2asXRE{ises, ii} = cat(2,cossimIC2XRE1, cossimIC2XRE2);
        
    end
    
    toc(dc)
end

%% inference cosine similarity across sessions
prctTRE1asIC1vsLC1 = zeros(Nsessions, numel(disperses));
prctTRE2asIC2vsLC2 = zeros(Nsessions, numel(disperses));
prctIC1asXRE1vsXRE2 = zeros(Nsessions, numel(disperses));
prctIC2asXRE2vsXRE1 = zeros(Nsessions, numel(disperses));
for ii = 1:numel(disperses)
    for ises = 1:Nsessions
prctTRE1asIC1vsLC1(ises, ii) = mean( cossimTRE1asIC1orLC1{ises, ii}(:,1)>cossimTRE1asIC1orLC1{ises, ii}(:,2) );
prctTRE2asIC2vsLC2(ises, ii) = mean( cossimTRE2asIC2orLC2{ises, ii}(:,1)>cossimTRE2asIC2orLC2{ises, ii}(:,2) );

prctIC1asXRE1vsXRE2(ises, ii) = mean( cossimIC1asXRE{ises, ii}(:,1)>cossimIC1asXRE{ises, ii}(:,2) );
prctIC2asXRE2vsXRE1(ises, ii) = mean( cossimIC2asXRE{ises, ii}(:,1)<cossimIC2asXRE{ises, ii}(:,2) );
    end
end

prctTREasICvsLC = mean(cat(3,prctTRE1asIC1vsLC1, prctTRE2asIC2vsLC2),3);
prctICasXRE = mean(cat(3,prctIC1asXRE1vsXRE2, prctIC2asXRE2vsXRE1),3);

xl = [disperses(1) disperses(end)];
yl= [0 1];
anot = sprintf('Open Scope ICwcfg1 high-repetition trials: mean vs each trial');
figure('Position', [0 0 1000 750])
annotation('textbox', [0.1 0.9 0.8 0.1], 'String', anot, 'edgecolor', 'none','FontSize', 12)
for isp = 1:6
switch isp
    case 1
        tempprct = prctTRE1asIC1vsLC1;
        ylab = 'P(cossim(TRE1,IC1)>cossim(TRE1,LC1))';
    case 2
        tempprct = prctTRE2asIC2vsLC2;
        ylab = 'P(cossim(TRE2,IC2)>cossim(TRE2,LC2))';
    case 3
        tempprct = prctTREasICvsLC;
        ylab = 'P(cossim(TRE,IC)>cossim(TRE,LC))';
    case 4
        tempprct = prctIC1asXRE1vsXRE2;
        ylab = 'P(cossim(IC1,XRE1)>cossim(IC1,XRE2))';
    case 5
        tempprct = prctIC2asXRE2vsXRE1;
        ylab = 'P(cossim(IC2,XRE2)>cossim(IC2,XRE1))';
    case 6
        tempprct = prctICasXRE;
        ylab = 'P(cossim(IC,XREiso)>cossim(IC,XREorth))';
end
subplot(2,3,isp)
hold all
plot(xl, [0.5 0.5], 'k--')
if pltsesavg
plot(disperses, tempprct, '-', 'Color', [0.5 0.5 0.5])
errorbar(disperses, mean(tempprct,1), std(tempprct,0,1)/sqrt(Nsessions), 'bx-', 'linewidth', 1)
[mv,mi]=max(mean(tempprct,1));
plot(disperses(mi)*[1 1], yl, 'r--')
text(disperses(mi), yl(1), num2str(disperses(mi)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'r')
else
plot(disperses, tempprct, '.-');
end
ylim(yl)
xlim(xl)
set(gca, 'XGrid', 'on', 'FontSize', fs)
xlabel('log(mean) vs log(var) slope')
ylabel(ylab)
end


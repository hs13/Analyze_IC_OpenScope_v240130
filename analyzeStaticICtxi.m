function ICsig = analyzeStaticICtxi(R, trialorder)

Nneurons = size(R, 1);
if size(R,2)~=length(trialorder)
    error('check R dimensions')
end

Palpha = 0.05;
correctmultcomp = true;
verbosity = false;

%% Kruskalwallis statistics
% perhaps do ranksum tests instead of
clk = tic;
tic
% PkwBKI and PmcBKI can be skipped if time is scarce (will save ~15 seconds)
BKtt = [0 106 107 110 111]; % Blank and Kanizsa trial type
% BKItt = [0 106 107 110 111 1301:1308]; % Blank, Kanizsa and Inducer trial type
BItt = [0 1301:1308]; % Blank and Inducer trial type
BICREltt = [0 106 506 111 511]; % Blank, IC and REl trial type
BICREl1tt = [0 106 506];
BICREl2tt = [0 111 511];
for itt = 1:5
    switch itt
        case 1
            temptt = BKtt;
        case 2
            temptt = BItt;
        case 3
            continue
            %             temptt = BKItt;
            temptt = BICREltt;
        case 4
            temptt = BICREl1tt;
        case 5
            temptt = BICREl2tt;
    end
    temptrialsoi = ismember(trialorder, temptt);
    temptt = unique(trialorder(temptrialsoi));
    
        tempPkw = zeros(Nneurons,1);
        tempPmc = NaN(Nneurons, nchoosek(numel(temptt),2));
        for ci = 1:Nneurons
            [p,tbl,stats] = kruskalwallis(R(ci,temptrialsoi), ...
                trialorder(temptrialsoi), 'off');
            tempPkw(ci) = p;
            % if p<0.05
            c = multcompare(stats, 'Display','off');
            tempPmc(ci,:) = c(:,end);
            % end
        end
        tempttpair = [temptt(c(:,1))' temptt(c(:,2))'];

    switch itt
        case 1
            PkwBK = tempPkw;
            PmcBK = tempPmc;
            BKttpair = tempttpair;
        case 2
            PkwBI = tempPkw;
            PmcBI = tempPmc;
            BIttpair = tempttpair;
        case 3
            %             PkwBKI = tempPkw;
            %             PmcBKI = tempPmc;
            %             BKIttpair = tempttpair;
            PkwBICREl = tempPkw;
            PmcBICREl = tempPmc;
            BICRElttpair = tempttpair;
        case 4
            PkwBICREl1 = tempPkw;
            PmcBICREl1 = tempPmc;
            BICREl1ttpair = tempttpair;
        case 5
            PkwBICREl2 = tempPkw;
            PmcBICREl2 = tempPmc;
            BICREl2ttpair = tempttpair;
    end
end
if verbosity
toc
fprintf('did %d kruskal wallis tests\n', itt)
end

blanktrials = trialorder==0;
if verbosity
tic
end
BKttpairrows = find(BKttpair(:,1)==0);
SP_BK = NaN(Nneurons, numel(BKttpairrows) );
Pmww_BK = NaN(Nneurons, numel(BKttpairrows) );
for typi = 1:numel(BKttpairrows)
    Ktrials = trialorder==BKttpair(BKttpairrows(typi),2);
    if nnz(blanktrials)==0 || nnz(Ktrials)==0
        continue
    end
    labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(Ktrials))];
    
    for ci = 1:Nneurons
        scores = [R(ci, blanktrials) R(ci, Ktrials)];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        SP_BK(ci,typi) = AUC;
        
        Pmww_BK(ci, typi) = ranksum(R(ci, blanktrials), R(ci, Ktrials));
    end
end
if verbosity
toc
disp('calculated SP_BK')
end

% sigmcBK = false(size(PmcBK));
% sigmcBK(PmcBK<0.05) = true;
sigmcBK = zeros( Nneurons, size(SP_BK,2) );
sigmcBK(SP_BK>0.5) = 1;
sigmcBK(SP_BK<0.5) = -1;
sigmcBK(PmcBK(:, BKttpair(:,1)==0)>=Palpha) = 0;
sigmcBK(PkwBK>=Palpha, :)=0;

%{
tic
SP_BI = NaN(Nneurons, 8);
for typi = 1:8
    Itt = 1300+typi;
    Itrials = trialorder==Itt;
    
    for ci = 1:Nneurons
        scores = [R(ci, blanktrials) R(ci, Itrials)];
        labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(Itrials))];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        SP_BI(ci,typi) = AUC;
    end
end
SP_BKI = [SP_BK SP_BI];
toc
if verbosity
disp('calculated SP_BKI')
end
sigmcBKI = zeros( Nneurons, size(SP_BKI,2) );
sigmcBKI(SP_BKI>0.5) = 1;
sigmcBKI(SP_BKI<0.5) = -1;
sigmcBKI(PmcBKI(:, 1:12)>=Palpha) = 0;
sigmcBKI(PkwBKI>=Palpha, :)=0;
save(strcat(onlineanalpath, 'SP_BKI.mat'), 'SP_BKI', 'sigmcBKI')
%}

if verbosity
tic
end
SP_Ind = NaN(Nneurons, 4);
Pmww_Ind = NaN(Nneurons, 4);
for typi = 1:4
    intt = 1300+typi;
    outtt = 1304+typi;
    intrials = trialorder==intt;
    outtrials = trialorder==outtt;
    if nnz(outtrials)==0 || nnz(intrials)==0
        continue
    end
    labels = [zeros(1,nnz(outtrials)) ones(1,nnz(intrials))];
    
    for ci = 1:Nneurons
        scores = [R(ci, outtrials) R(ci, intrials)];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        SP_Ind(ci,typi) = AUC;
        
        Pmww_Ind(ci,typi) = ranksum(R(ci, outtrials), R(ci, intrials));
    end
end
if verbosity
toc
disp('calculated SP_Ind')
end
% save(strcat(onlineanalpath, 'SP_Ind.mat'), 'SP_Ind')

%
if verbosity
tic
end
SP_BICREl = NaN(Nneurons, numel(BICREltt)-1);
Pmww_BICREl = NaN(Nneurons, numel(BICREltt)-1);
for typi = 2:numel(BICREltt)
    temptrials = trialorder==BICREltt(typi);
    if nnz(blanktrials)==0 || nnz(temptrials)==0
        continue
    end
    labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(temptrials))];
    
    for ci = 1:Nneurons
        scores = [R(ci, blanktrials) R(ci, temptrials)];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        SP_BICREl(ci,typi-1) = AUC;
        
        Pmww_BICREl(ci,typi) = ranksum(R(ci, blanktrials), R(ci, temptrials));
    end
end
if verbosity
toc
disp('calculated SP_BICREl')
end

% sigmcBICREl = zeros( Nneurons, size(SP_BICREl,2) );
% sigmcBICREl(SP_BICREl>0.5) = 1;
% sigmcBICREl(SP_BICREl<0.5) = -1;
% sigmcBICREl(PmcBICREl(:, 1:4)>=Palpha) = 0;
% sigmcBICREl(PkwBICREl>=Palpha, :)=0;

SP_BICREl1 = SP_BICREl(:, ismember(BICREltt(2:end), BICREl1tt(2:end)));
sigmcBICREl1 = zeros(size(SP_BICREl1));
sigmcBICREl1(SP_BICREl1>0.5) = 1;
sigmcBICREl1(SP_BICREl1<0.5) = -1;
sigmcBICREl1(PmcBICREl1(:, BICREl1ttpair(:,1)==0)>=Palpha) = 0;
sigmcBICREl1(PkwBICREl1>=Palpha, :)=0;

SP_BICREl2 = SP_BICREl(:, ismember(BICREltt(2:end), BICREl2tt(2:end)));
sigmcBICREl2 = zeros(size(SP_BICREl2));
sigmcBICREl2(SP_BICREl2>0.5) = 1;
sigmcBICREl2(SP_BICREl2<0.5) = -1;
sigmcBICREl2(PmcBICREl2(:, BICREl2ttpair(:,1)==0)>=Palpha) = 0;
sigmcBICREl2(PkwBICREl2>=Palpha, :)=0;

if verbosity
tic
end
SP_ICvsRC = NaN(Nneurons, 1);
Pmww_ICvsRC = NaN(Nneurons, 1);
for ci = 1:Nneurons
    if SP_BK(ci,1)>=SP_BK(ci,4)
        tempICtt=106;
    else
        tempICtt=111;
    end
    if SP_BK(ci,2)>=SP_BK(ci,3)
        tempRCtt=107;
    else
        tempRCtt=110;
    end
    tempICtrials = trialorder==tempICtt;
    tempRCtrials = trialorder==tempRCtt;
    
    if nnz(tempICtrials)==0 || nnz(tempRCtrials)==0
        continue
    end
    
    scores = [R(ci, tempRCtrials) R(ci, tempICtrials)];
    labels = [zeros(1,nnz(tempRCtrials)) ones(1,nnz(tempICtrials))];
    [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
    SP_ICvsRC(ci) = AUC;
    
    Pmww_ICvsRC(ci) = ranksum(R(ci, tempRCtrials), R(ci, tempICtrials));
end
if verbosity
toc
disp('calculated SP_ICvsRC')
end

%% IC vs RC encoder neuron
ICencoder = (PkwBK<Palpha) & (sigmcBK(:,1)==1 | sigmcBK(:,4)==1) ...
    & (sigmcBK(:,2)==0) & (sigmcBK(:,3)==0);
RCencoder = (PkwBK<Palpha) & (sigmcBK(:,2)==1 | sigmcBK(:,3)==1) ...
    & (sigmcBK(:,1)==0) & (sigmcBK(:,4)==0);
inducerencoder = (PkwBK<Palpha) & xor(sigmcBK(:,1)==1, sigmcBK(:,4)==1) ...
    & xor(sigmcBK(:,2)==1, sigmcBK(:,3)==1);
inducerresponsive = inducerencoder | ( PkwBK<Palpha & all(sigmcBK(:,1:4)==1, 2) );

% if verbosity
fprintf('IC-encoder %d, RC-encoder %d, inducer-encoder %d inducer resposive %d\n', ...
    nnz(ICencoder), nnz(RCencoder), nnz(inducerencoder), nnz(inducerresponsive) )
% end

% save(strcat(onlineanalpath, 'SP_BK.mat'), 'SP_BK', 'sigmcBK', ...
%     'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive')

%% indin and indout cells
indenc1 = PkwBK<Palpha & ismember(sigmcBK, [1 1 0 0], 'rows');
indenc2 = PkwBK<Palpha & ismember(sigmcBK, [0 0 1 1], 'rows');
indenc3 = PkwBK<Palpha & ismember(sigmcBK, [1 0 1 0], 'rows');
indenc4 = PkwBK<Palpha & ismember(sigmcBK, [0 1 0 1], 'rows');
if nnz(ICencoder)
    if ~all(ismember( (indenc1 | indenc2 | indenc3 | indenc4), ICencoder))
        error('check inducer encoder definition')
    end
end

if verbosity
fprintf('indenc1 %d, indenc2 %d, indenc3 %d indenc4 %d\n', ...
    nnz(indenc1), nnz(indenc2), nnz(indenc3), nnz(indenc4) )
end

for typi = 1:4
    intt = 1300+typi;
    outtt = 1304+typi;
%     tempmcind = BIttpair(:,1)==intt & BIttpair(:,2)==outtt;
%     if ~(nnz(tempmcind)==1)
%         error('check BI comparisons')
%     end
    %     tempindin = SP_Ind(:,typi)>0.5 & PkwBI<Palpha & PmcBI(:,tempmcind)<Palpha;
    %     tempindout = SP_Ind(:,typi)<0.5 & PkwBI<Palpha & PmcBI(:,tempmcind)<Palpha;
    % HS 220518: IMPORTANT CHANGE
    if correctmultcomp
        tempindin = SP_Ind(:,typi)>0.5 & PkwBI<Palpha & Pmww_Ind(:,typi)<Palpha & ~ICencoder & ~RCencoder;
        tempindout = SP_Ind(:,typi)<0.5 & PkwBI<Palpha & Pmww_Ind(:,typi)<Palpha & ~ICencoder & ~RCencoder;
    else
        % below is an even more inclusive definition. use if Nvalind is less than 20
        tempindin = SP_Ind(:,typi)>0.5 & Pmww_Ind(:,typi)<Palpha & ~ICencoder & ~RCencoder;
        tempindout = SP_Ind(:,typi)<0.5 & Pmww_Ind(:,typi)<Palpha & ~ICencoder & ~RCencoder;
    end
    
if verbosity
    fprintf('indin%d N=%d SP values\n', typi, nnz(tempindin))
    %     disp(SP_Ind(tempindin,typi))
    %     disp('valid indin SP values')
    %     disp(SP_Ind(tempvalindin,typi))
    
    fprintf('indout%d N=%d SP values\n', typi, nnz(tempindout))
    %     disp(SP_Ind(tempindout,typi))
    %     disp('valid indout SP values')
    %     disp(SP_Ind(tempvalindout,typi))
end

    switch typi
        case 1
            indin1 = tempindin;
            indout1 = tempindout;
        case 2
            indin2 = tempindin;
            indout2 = tempindout;
        case 3
            indin3 = tempindin;
            indout3 = tempindout;
        case 4
            indin4 = tempindin;
            indout4 = tempindout;
    end
end

%% IC/RC-encoder, RElfaith neurons
ICencoder1 = (PkwBK<Palpha) & (sigmcBK(:,1)==1) & (sigmcBK(:,2)==0) & (sigmcBK(:,3)==0);
ICencoder2 = (PkwBK<Palpha) & (sigmcBK(:,4)==1) & (sigmcBK(:,2)==0) & (sigmcBK(:,3)==0);
RCencoder1 = (PkwBK<Palpha) & (sigmcBK(:,2)==1) & (sigmcBK(:,1)==0) & (sigmcBK(:,4)==0);
RCencoder2 = (PkwBK<Palpha) & (sigmcBK(:,3)==1) & (sigmcBK(:,1)==0) & (sigmcBK(:,4)==0);

indenc13 = (PkwBK<Palpha) & (sigmcBK(:,1)==1) & ~ICencoder1;
indenc24 = (PkwBK<Palpha) & (sigmcBK(:,4)==1) & ~ICencoder2;
indenc14 = (PkwBK<Palpha) & (sigmcBK(:,2)==1) & ~RCencoder1;
indenc23 = (PkwBK<Palpha) & (sigmcBK(:,3)==1) & ~RCencoder2;

if verbosity
fprintf('ICencoder1 %d, RCencoder1 %d, RCencoder2 %d ICencoder2 %d\n', ...
    nnz(ICencoder1), nnz(RCencoder1), nnz(RCencoder2), nnz(ICencoder2) )

fprintf('indenc13 %d, indenc14 %d, indenc23 %d indenc24 %d\n', ...
    nnz(indenc13), nnz(indenc14), nnz(indenc23), nnz(indenc24) )
end


% RElICenc1 = (PkwBICREl1<Palpha) & (sigmcBICREl1(:,1)==1) & (sigmcBICREl1(:,2)==1);
RElfaith1 = (PkwBICREl1<Palpha) & (sigmcBICREl1(:,1)==0) & (sigmcBICREl1(:,2)==1);
% RElICenc2 = (PkwBICREl2<Palpha) & (sigmcBICREl2(:,1)==1) & (sigmcBICREl2(:,2)==1);
RElfaith2 = (PkwBICREl2<Palpha) & (sigmcBICREl2(:,1)==0) & (sigmcBICREl2(:,2)==1);

% RElresp can include IC encoders
RElresp1 = (PkwBICREl1<Palpha) & (sigmcBICREl1(:,2)==1);
RElresp2 = (PkwBICREl2<Palpha) & (sigmcBICREl2(:,2)==1);

if verbosity
fprintf('RElfaith1 %d, RElfaith2 %d\n',nnz(RElfaith1), nnz(RElfaith2) )
end

%%
% select the same number of neurons and compare top IC/RC-responders vs top IC-RC-tuned neurons
% 	- ICresp/RCresp: sort IC/RC-responders based on AUROC discriminating IC/RC vs blank trials (SP)
% 	- ICtuned/RCtuned: sort IC/RC-tuned neurons based on AUROC discriminating most-responsive IC trials vs most-responsive RC trials
% 		- e.g., if a neuron has SP_IC1>SP_IC2 and SP_RC2>SP_RC1, the "tuning" would be determined based on IC1 vs RC2 trials
% 	- MATCH NUMBER and allow overlap between the two groups


[SP_prefICRC, ICRCtuneind] = max(SP_BK,[],2);

if correctmultcomp
    ICresp1 = SP_BK(:,1)>0.5 & Pmww_BK(:,1)<Palpha & PkwBK<Palpha;
    RCresp1 = SP_BK(:,2)>0.5 & Pmww_BK(:,2)<Palpha & PkwBK<Palpha;
    RCresp2 = SP_BK(:,3)>0.5 & Pmww_BK(:,3)<Palpha & PkwBK<Palpha;
    ICresp2 = SP_BK(:,4)>0.5 & Pmww_BK(:,4)<Palpha & PkwBK<Palpha;
    
    ICtuned1 = ICRCtuneind==1 & SP_prefICRC>0.5 & SP_ICvsRC>0.5 & PkwBK<Palpha;
    RCtuned1 = ICRCtuneind==2 & SP_prefICRC>0.5 & SP_ICvsRC<0.5 & PkwBK<Palpha;
    RCtuned2 = ICRCtuneind==3 & SP_prefICRC>0.5 & SP_ICvsRC<0.5 & PkwBK<Palpha;
    ICtuned2 = ICRCtuneind==4 & SP_prefICRC>0.5 & SP_ICvsRC>0.5 & PkwBK<Palpha;
else
    ICresp1 = SP_BK(:,1)>0.5 & Pmww_BK(:,1)<Palpha;
    RCresp1 = SP_BK(:,2)>0.5 & Pmww_BK(:,2)<Palpha;
    RCresp2 = SP_BK(:,3)>0.5 & Pmww_BK(:,3)<Palpha;
    ICresp2 = SP_BK(:,4)>0.5 & Pmww_BK(:,4)<Palpha;
    
    ICtuned1 = ICRCtuneind==1 & SP_prefICRC>0.5 & SP_ICvsRC>0.5;
    RCtuned1 = ICRCtuneind==2 & SP_prefICRC>0.5 & SP_ICvsRC<0.5;
    RCtuned2 = ICRCtuneind==3 & SP_prefICRC>0.5 & SP_ICvsRC<0.5;
    ICtuned2 = ICRCtuneind==4 & SP_prefICRC>0.5 & SP_ICvsRC>0.5;
end
% ICtuned1 = ICRCtuneind==1 & SP_prefICRC>0.5 & SP_ICvsRC>0.5 & Pmww_ICvsRC<Palpha;
% RCtuned1 = ICRCtuneind==2 & SP_prefICRC>0.5 & SP_ICvsRC<0.5 & Pmww_ICvsRC<Palpha;
% RCtuned2 = ICRCtuneind==3 & SP_prefICRC>0.5 & SP_ICvsRC<0.5 & Pmww_ICvsRC<Palpha;
% ICtuned2 = ICRCtuneind==4 & SP_prefICRC>0.5 & SP_ICvsRC>0.5 & Pmww_ICvsRC<Palpha;

if verbosity
fprintf('ICresp1 %d, RCresp1 %d, RCresp2 %d ICresp2 %d\n', ...
    nnz(ICresp1), nnz(RCresp1), nnz(RCresp2), nnz(ICresp2) )
fprintf('ICtuned1 %d, RCtuned1 %d, RCtuned2 %d ICtuned2 %d\n', ...
    nnz(ICtuned1), nnz(RCtuned1), nnz(RCtuned2), nnz(ICtuned2) )
fprintf('tuned & resp overlap (out of tuned): %.0f%% %.0f%% %.0f%% %.0f%%\n', ...
    100*mean(ICtuned1(ICresp1)), 100*mean(RCtuned1(RCresp1)), ...
    100*mean(RCtuned2(RCresp2)), 100*mean(ICtuned2(ICresp2)) )

fprintf('IC/RC-encoder & tuned overlap (out of encoders): %.0f%% %.0f%% %.0f%% %.0f%%\n', ...
    100*mean(ICtuned1(ICencoder1)), 100*mean(RCtuned1(RCencoder1)), ...
    100*mean(RCtuned2(RCencoder2)), 100*mean(ICtuned2(ICencoder2)) )
end

%%
%{
save(strcat(pathpp, 'SP_BKInd.mat'), 'Palpha', 'SP_Ind', 'Pmww_Ind', ...
    'SP_BK', 'sigmcBK', 'Pmww_BK', 'SP_ICvsRC', 'Pmww_ICvsRC', ...
    'SP_BICREl', 'sigmcBICREl1', 'sigmcBICREl2')
save(strcat(pathpp, 'ICsig.mat'), ... % 'whichdef', ...
    'BKtt', 'BKttpair', 'PkwBK', 'PmcBK', ...
    'BItt', 'BIttpair', 'PkwBI', 'PmcBI', ...
    'BICREl1tt', 'BICREl1ttpair', 'PkwBICREl1', 'PmcBICREl1', ...
    'BICREl2tt', 'BICREl2ttpair', 'PkwBICREl2', 'PmcBICREl2', ...
    'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive')
%     'BKItt', 'BKIttpair', 'PkwBKI', 'PmcBKI', ...
%     'BICREltt', 'BICRElttpair', 'PkwBICREl', 'PmcBICREl', ...
save(strcat(pathpp, 'neurongroups.mat'), ...
    'indenc1', 'indenc2', 'indenc3', 'indenc4', ...
    'indenc13', 'indenc14', 'indenc23', 'indenc24', ...
    'ICencoder1', 'RCencoder1', 'RCencoder2', 'ICencoder2', ...
    'indin1', 'indin2', 'indin3', 'indin4', ...
    'indout1', 'indout2', 'indout3', 'indout4', ...
    'RElfaith1', 'RElfaith2', ...
    'ICresp1', 'RCresp1', 'RCresp2', 'ICresp2', ...
    'ICtuned1', 'RCtuned1', 'RCtuned2', 'ICtuned2')
%}

ICsigfields = {'Palpha', 'SP_Ind', 'Pmww_Ind', ...
    'SP_BK', 'sigmcBK', 'Pmww_BK', 'SP_ICvsRC', 'Pmww_ICvsRC', ...
    'SP_BICREl', 'Pmww_BICREl', 'sigmcBICREl1', 'sigmcBICREl2', ...
    'BKtt', 'BKttpair', 'PkwBK', 'PmcBK', ...
    'BItt', 'BIttpair', 'PkwBI', 'PmcBI', ...
    'BICREl1tt', 'BICREl1ttpair', 'PkwBICREl1', 'PmcBICREl1', ...
    'BICREl2tt', 'BICREl2ttpair', 'PkwBICREl2', 'PmcBICREl2', ...
    'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive', ...
    'indenc1', 'indenc2', 'indenc3', 'indenc4', ...
    'indenc13', 'indenc14', 'indenc23', 'indenc24', ...
    'ICencoder1', 'RCencoder1', 'RCencoder2', 'ICencoder2', ...
    'indin1', 'indin2', 'indin3', 'indin4', ...
    'indout1', 'indout2', 'indout3', 'indout4', ...
    'RElfaith1', 'RElfaith2', ...
    'ICresp1', 'RCresp1', 'RCresp2', 'ICresp2', ...
    'ICtuned1', 'RCtuned1', 'RCtuned2', 'ICtuned2'};
ICsig = struct();
for f = 1:numel(ICsigfields)
    ICsig.(ICsigfields{f}) = eval(ICsigfields{f});
end

toc(clk)
end % end of function
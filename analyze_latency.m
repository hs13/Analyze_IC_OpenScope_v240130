% 1. replication of 2p results
% 	1. population average activity for all neurons and ctr-CRF neurons not different between IC and RC across black and white contrasts
% 	2. proportion of IC-encoders and RC-encoders are not significantly different (out of all/visually-responsive neurons and out of ctr-CRF neurons)
% 	3. preferred orientation match between ICs and gratings for IC-encoders and RC-encoders
% 	4. ctrCRF distribution of IC-encoders and RC-encoders are not biased to the illusory region
% 	5. decoder results
% 	6. OSI and iRF
% 2. latency analysis
% 	1. compare latency for ctr-CRF neurons for IC vs REl/REt across visual areas
% 3. synchrony analysis
% 	1. spike synchrony
% 	2. spike spectrum
% 	3. spike-field coherence
% 	4. lfp anaylsis
% 4. functional connectivity analysis
% 	1. CCG based driver/ follower analysis Jia et al 2022 https://www.sciencedirect.com/science/article/pii/S0896627322000848?via%3Dihub#undfig1

% from RFCI block, get RF
% from sizeCI block, get size tuning and orientation tuning

% for each unit, visual block and visual trial type, find the mode in
% 1/5/10ms bins (also do bootstrapped/jack-knifed CI)

addpath(genpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v230821'))

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
% visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
%     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
% probes = {'C', 'D'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

datadir = 'S:\OpenScopeData\00248_v230821\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
    clearvars -except probes visblocks datadir nwbsessions Nsessions ises
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];

    for iprobe = 1:numel(probes)
        tic
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
        %         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'

        vistrialtypes = struct();
        vistrialrep = struct();
        spklatency = struct();
        spklatencyprob = struct();
        spklatencyadj = struct();
        spklatencyadjprob = struct();
        spklatencybins = [1 5 10 25]; % bin size for smoothing histogram
        spklatencyT0s = [0 5 10]; % minimum spike latency
        Nneurons = numel(neuoind);
        for b = 1:numel(visblocks)
            vistrialtypes.(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
            Ntt = numel(vistrialtypes.(visblocks{b}));
            vistrialrep.(visblocks{b}) = zeros(size(vistrialtypes.(visblocks{b})));
            for tt = 1:numel(spklatencyT0s)
            for ibin = 1:numel(spklatencybins)
                whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
                spklatency.(whichbin).(visblocks{b}) = NaN(Nneurons,Ntt);
                spklatencyprob.(whichbin).(visblocks{b}) = zeros(Nneurons,Ntt);
                spklatencyadj.(whichbin).(visblocks{b}) = NaN(Nneurons,Ntt);
                spklatencyadjprob.(whichbin).(visblocks{b}) = zeros(Nneurons,Ntt);
            end
            end
            for ii = 1:Ntt
                trialsoi = vis.(visblocks{b}).trialorder==vistrialtypes.(visblocks{b})(ii);
                vistrialrep.(visblocks{b})(ii) = nnz(trialsoi);
                for ci = 1:Nneurons

                    for tt = 1:numel(spklatencyT0s)
                        temptloi = psthtli>spklatencyT0s(tt) & psthtli<=250;
                        temptli = psthtli(temptloi);

                        % temppsth is Ntimepoints X Ntrials
                        temppsth = squeeze(psth.(visblocks{b})(temptloi,trialsoi,ci));
                        if nnz(temppsth)==0
                            continue
                        end

                        temppsthpre = squeeze(psth.(visblocks{b})(psthtli>-250 & psthtli<=0,trialsoi,ci));
                        [rpre,cpre]=find(temppsthpre & cumsum(temppsthpre,1)==1);
                        hcpre = histcounts(rpre, 0.5:1:nnz(temptloi)+0.5);

                        % find the first spike on each trial
                        [r,c]=find(temppsth & cumsum(temppsth,1)==1);
                        % spklatency.bin1ms.(visblocks{b})(ci,ii) = mode(r);

                        % % alternative way to find the first spike on each trial
                        % %[~,mi] = max(temppsth,[],1);
                        % %isequal(mi(c)', r)

                        hc = histcounts(r, 0.5:1:nnz(temptloi)+0.5);

% NOTE, SMOOTHING IS THE SAME AS HISTOGRAM WITH CORRESPONDING BIN SIZE ONLY
% WHEN BIN SIZE IS AN ODD NUMBER (even numbers are reduced by 1)
% cf. compare_latency_smoothVShistogram

                        %bw = 5; figure; hold all; histogram(r, 'binwidth', bw); plot(bw*smooth(hc,bw), 'm-')
                        for ibin = 1:numel(spklatencybins)
                            whichbin = sprintf('plus%dms_bin%dms', spklatencyT0s(tt), spklatencybins(ibin));
                            hcsm = smooth( hc, spklatencybins(ibin) );
                            [mv,mi] = max(hcsm);
                            if spklatencybins(ibin)==1
                                if mi ~= mode(r)
                                    error('check spike latency calculation for 1ms bin')
                                end
                            end
                            spklatency.(whichbin).(visblocks{b})(ci,ii) = temptli(mi);
                            spklatencyprob.(whichbin).(visblocks{b})(ci,ii) = mv/nnz(trialsoi);

                            hcpresm = smooth( hcpre, spklatencybins(ibin) );
                            hcadjsm = hcsm-hcpresm;
                            %hcnormsm = (hcsm+0.1)./(hcpresm+0.1);
                            [mv,mi] = max(hcadjsm);
                            spklatencyadj.(whichbin).(visblocks{b})(ci,ii) = temptli(mi);
                            spklatencyadjprob.(whichbin).(visblocks{b})(ci,ii) = mv/nnz(trialsoi);
                        end
                        % adjust spike latency by dividing by the prestimulus
                        % distribution (treat prestimulus as 'null')
                    end
                end
            end
        end
        save(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{iprobe}), ...
            'vistrialtypes', 'vistrialrep', 'spklatencybins', ...
            'spklatency', 'spklatencyprob', 'spklatencyadj', 'spklatencyadjprob')
        toc
    end
end


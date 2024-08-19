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
Nneuronsvec = csz(:,2);
Nneurons = sum(csz(:,2));

spkcntICttaccavg = cellfun(@mean, spkcntICttacc, 'uniformoutput', false);
spkcntICttaccvar = cellfun(@var, spkcntICttacc, 'uniformoutput', false);

spkcntICttavg = NaN(Nneurons, Ntt);
spkcntICttvar = NaN(Nneurons, Ntt);
for ii = 1:Ntt
    spkcntICttavg(:,ii) = cat(2,spkcntICttaccavg{:,ii});
    spkcntICttvar(:,ii) = cat(2,spkcntICttaccvar{:,ii});
end

clearvars spkcntICttacc spkcntICttaccavg spkcntICttaccvar

%% effective dimensionality : should we divide by N or no?
% dividing by sqrt(N) seems to eliminate dependence on N
Nneu = max(Nneuronsvec);
Nvec = 10:10:Nneu;
Dmat = NaN(Nsessions, length(Nvec));
for ises = 1:Nsessions % or 5
    Nneu = Nneuronsvec(ises);
    %[2.^(3:floor(log2(Nneu))) Nneu]
    xvec = 10:10:Nneu;
    yvec = NaN(size(xvec));
    for ix = 1:numel(xvec)
        neuoi = randperm(Nneu,xvec(ix));
        tempcovtt = cov(spkcntICttagg{ises}{1}(:,neuoi));

        lambdas = sort( eig(tempcovtt), 'descend');
        yvec(ix) = (sum(lambdas))^2/sum(lambdas.^2);
        %yvec(ix) = (mean(lambdas))^2/mean(lambdas.^2);
    end
    Dmat(ises,1:length(xvec)) = yvec;
end

figure; plot(Nvec, Dmat)

% figure; plot(Nvec, Dmat./Nvec)
% figure; plot(Nvec, Dmat./sqrt(Nvec))

%%

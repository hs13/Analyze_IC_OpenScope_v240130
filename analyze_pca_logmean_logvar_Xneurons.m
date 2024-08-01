% define alpha as slope in log(PC rank) vs log(PC explaned variance)
% define DS (dispersion slope) as slope in log(spike count mean) vs
% log(spike count variance) ACROSS NEURONS, WITHIN TRIAL TYPE
% what is the relationship between alpha and DS?
% alpha seems to be around -1 when DS is 1.5
% in real data (400 repeats) DS seems to be around 1 

% next steps
% need to confirm that DS changes with RR0 and shuftrials
% do numerical simulation with neuropixels natural scenes data (pooled across sessions)

pathpp = 'S:\OpenScopeData\00248_v240130\postprocessed\sub-619296\';
load([pathpp, 'postprocessed.mat'])
load([pathpp, 'qc_units.mat'])

neuV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
neuLM = contains(neuallloc, 'VISl');
neuAL = contains(neuallloc, 'VISal');
neuRS = unit_wfdur>0.4;
neufilt = (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
neuinarea = neuV1 & neuRS & neufilt;

whichblock = 'ICwcfg1_presentations';
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
blocktrialorder = ICtrialtypes(vis.(whichblock).trialorder+1);

figure; 
for itt = 1:numel(ICtrialtypes)
    subplot(4,6,itt)
tempspkcnt= Rall.(whichblock)(blocktrialorder==ICtrialtypes(itt), neuinarea)*0.4;
tempx = log10(mean(tempspkcnt,1));
tempy = log10(var(tempspkcnt,0,1));
plot(tempx, tempy, '.')
valneu = mean(tempspkcnt,1)>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
plot(xl, a*xl+b, 'g-')
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
title(ICtrialtypes(itt))
xlabel('log10(mean)')
ylabel('log10(var)')
end


indsoi = 1:floor(nnz(neuinarea)/2);
figure; 
for itt = 1:numel(ICtrialtypes)
    subplot(4,6,itt)
tempspkcnt= Rall.(whichblock)(blocktrialorder==ICtrialtypes(itt), neuinarea)*0.4;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(tempspkcnt);
tempx = log10(indsoi);
tempy = log10(EXPLAINED(indsoi));
plot(tempx, tempy, '.')
LM = fitlm(tempx, tempy ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
plot(xl, a*xl+b, 'g-')
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
title(ICtrialtypes(itt))
xlabel('PC rank')
ylabel('Explained Variance')
end

% Rows of R111 correspond to observations and columns to variables.
spkcnttt = Rall.(whichblock)(blocktrialorder==111, neuinarea)*0.4;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(spkcnttt);
% Each column of COEFF contains coefficients for one principal component.
% PC1 coefficients COEFF(:,1)

figure; plot(mean(spkcnttt,1), MU, '.')
figure; plot(COEFF(:,1), MU, '.')

spkcnttt = Rall.(whichblock)(blocktrialorder==111, neuinarea)*0.4;
res = randn(size(spkcnttt));

%% vary fano factor and see if it changes PCA within trial type -- nope
ffs = [0.1 1 1.5 10];

% lognormal firing rate, mean firing rate = 1.83
nneu = 2000;
nrep = 10001;
spkcntmu = 3.^randn(1,nneu);
res = randn(nrep, nneu);

indsoi = 1:floor(nneu/2);
figure; 
for ii = 1:numel(ffs)
Rpoiss = spkcntmu + sqrt( spkcntmu*ffs(ii) ).*res;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(Rpoiss);

subplot(2,numel(ffs), ii)
tempx = log10(mean(Rpoiss,1));
tempy = log10(var(Rpoiss,0,1));
plot(tempx, tempy, '.')
valneu = mean(Rpoiss,1)>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
plot(xl, a*xl+b, 'g-')
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis([xl yl])
title(ffs(ii))
xlabel('log10(mean)')
ylabel('log10(var)')

subplot(2,numel(ffs), ii+numel(ffs))
plot(log10(indsoi), log10(EXPLAINED(indsoi)))
LM = fitlm(log10(indsoi), log10(EXPLAINED(indsoi))); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis([xl yl])
xlabel('PC rank')
ylabel('Explained Variance')
end

%% when log mean and log variance has slope of 1.5, PCA alpha is -1
% spkcntavg = mean(Rall.(whichblock)(blocktrialorder==0, neuinarea)*0.4, 1);
% 10.^mean(log10(spkcntavg(spkcntavg>0)))

disperses = [0.5 1 1.5 2];

% lognormal firing rate, mean firing rate = 1.83
nneu = 2000;
nrep = 10001;
spkcntmu = 3.^randn(1,nneu);
res = randn(nrep, nneu);

indsoi = 1:floor(nneu/2);
figure; 
for ii = 1:numel(disperses)
    tic
Rpoiss = spkcntmu + sqrt( 10.^( log10(spkcntmu)*disperses(ii) ) ).*res;
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(Rpoiss);

subplot(2,numel(disperses), ii)
tempx = log10(mean(Rpoiss,1));
tempy = log10(var(Rpoiss,0,1));
plot(tempx, tempy, '.')
valneu = mean(Rpoiss,1)>0;
LM = fitlm(tempx(valneu), tempy(valneu) ); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
plot(xl, a*xl+b, 'g-')
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis([xl yl])
title(disperses(ii))
xlabel('log10(mean)')
ylabel('log10(var)')

subplot(2,numel(disperses), ii+numel(disperses))
plot(log10(indsoi), log10(EXPLAINED(indsoi)))
LM = fitlm(log10(indsoi), log10(EXPLAINED(indsoi))); % y = ax+b
b=LM.Coefficients.Estimate(1);
a=LM.Coefficients.Estimate(2);
xl=xlim;yl=ylim; hold on;
text(xl(1), yl(1), sprintf('y = %.2fx + %.2f', a,b), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
axis([xl yl])
xlabel('PC rank')
ylabel('Explained Variance')
toc
end
% load('S:\OpenScopeData\00248_v240130\postprocessed\sub-620333\postprocessed.mat')
whichblock = 'ICwcfg1_presentations';
trialsoi = vis.(whichblock).trialorder==0;
neuinV1 = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
spkcnt0 = Rall.(whichblock)(trialsoi,neuinV1)*400/1000;

c2p = randperm(size(spkcnt0,2));
figure
for ii = 1:24
ci = c2p(ii);
[lambdamv, numv] = cmpFromMeanVariance(mean(spkcnt0(:,ci)), var(spkcnt0(:,ci)));

xbe = -0.5+(0:1:max(spkcnt0(:,ci))+1);
xbc = (xbe(1:end-1)+xbe(2:end))/2;
h = histcounts(spkcnt0(:,ci), 'BinEdges', xbe);

[lambda, nu, negLogL] = fitCMPDistribution(xbc, h);
[poisslam, poissexpected] = fitPoissonDistribution(xbc, h);

ymv = cmpTerm(xbc, lambdamv, numv);
ybc = cmpTerm(xbc, lambda, nu);
ypoiss = poisspdf(xbc,poisslam);

[lambdamv2, numv2] = cmpFromMeanVariance(mean(spkcnt0(:,ci)), var(spkcnt0(:,ci))^2);
ymv2 = cmpTerm(xbc, lambdamv2, numv2);

subplot(4,6,ii)
hold all
h = histogram(spkcnt0(:,ci), 'BinEdges', xbe, 'normalization', 'probability');
plot(xbc, ymv/sum(ymv), 'yo-', 'linewidth', 1)
plot(xbc, ybc/sum(ybc), 'g--', 'linewidth', 1)
plot(xbc, ypoiss/sum(ypoiss), 'b:', 'linewidth', 1)
plot(xbc, ymv2/sum(ymv2), 'r-', 'linewidth', 1)
legend({'histogram', 'CMP from mean & var', 'CMP fit', 'Poisson fit', 'CMP from mean & var^2'})

title(sprintf('Neuron%d mean %.4f variance %.4f variance^2 %.4f', ci, mean(spkcnt0(:,ci)), var(spkcnt0(:,ci)), var(spkcnt0(:,ci))^2 ))
end

%%
function [lambda, expected] = fitPoissonDistribution(x, y)
% fitPoissonDistribution fits a Poisson distribution to count data.
%
%   [lambda, expected] = fitPoissonDistribution(x, y)
%
% Inputs:
%   x - vector of count values (nonnegative integers)
%   y - corresponding observed frequencies for each count in x
%
% Outputs:
%   lambda   - estimated Poisson parameter (mean)
%   expected - expected frequencies computed from the Poisson model
%
% The maximum likelihood estimator for the Poisson mean (λ) is:
%   lambda = sum(x.*y) / sum(y)
%
% Once λ is computed, the expected frequencies for each x value are given by:
%   expected = sum(y) * poisspdf(x, lambda)
%
% Example:
%   x = 0:10;
%   y = [5, 15, 20, 10, 5, 3, 2, 1, 1, 0, 0];
%   [lambda, expected] = fitPoissonDistribution(x, y);
%   fprintf('Estimated lambda: %f\n', lambda);
%
% Note: This function uses MATLAB's poisspdf function available in the 
% Statistics Toolbox.

    % Check that x and y are vectors and have the same length
    if ~isvector(x) || ~isvector(y)
        error('Both x and y must be vectors.');
    end
    if length(x) ~= length(y)
        error('Vectors x and y must be of the same length.');
    end

    % Calculate the total number of observations
    totalCount = sum(y);
    
    % Estimate lambda using the MLE formula:
    %   lambda = (sum of (x * frequency)) / (total frequency)
    lambda = sum(x .* y) / totalCount;
    
    % Compute expected frequencies using the Poisson probability mass function
    expected = totalCount * poisspdf(x, lambda);
end

%%
function [lambda, nu, negLogL] = fitCMPDistribution(x, y)
% fitCMPDistribution fits a Conway-Maxwell-Poisson (CMP) distribution to count data.
%
%   [lambda, nu, negLogL] = fitCMPDistribution(x, y)
%
% Inputs:
%   x - vector of count values (nonnegative integers)
%   y - corresponding observed frequencies for each count in x
%
% Outputs:
%   lambda  - estimated rate parameter (λ) of the CMP distribution.
%   nu      - estimated dispersion parameter (ν) of the CMP distribution.
%   negLogL - negative log-likelihood value at the optimum.
%
% The CMP probability mass function is defined as:
%   P(X=k) = (lambda^k / (k!)^nu) / Z(lambda, nu)
% where Z(lambda, nu) = sum_{j=0}^∞ lambda^j/(j!)^nu.
%
% This function uses maximum likelihood estimation (MLE). The log-likelihood is:
%
%   logL = sum( y .* [ x*log(lambda) - nu*log(x!) ] ) - (sum(y)) * log(Z(lambda,nu))
%
% We use gammaln(x+1) in place of log(x!) for numerical stability.
%
% Example:
%   x = 0:10;
%   y = [5, 15, 20, 10, 5, 3, 2, 1, 1, 0, 0];
%   [lambda, nu, negLogL] = fitCMPDistribution(x, y);
%   fprintf('Estimated lambda: %f\nEstimated nu: %f\n', lambda, nu);

    % Validate inputs
    if ~isvector(x) || ~isvector(y)
        error('Both x and y must be vectors.');
    end
    if length(x) ~= length(y)
        error('Vectors x and y must be of the same length.');
    end

    % Use the Poisson case as an initial guess:
    % For Poisson, ν = 1 and λ equals the weighted mean.
    totalCount = sum(y);
    lambda0 = sum(x .* y) / totalCount;
    nu0 = 1;

    % Optimize in log-space to ensure parameters stay positive.
    initParams = [log(lambda0); log(nu0)];

    % Set optimization options; adjust tolerances if needed.
    options = optimset('Display', 'iter', 'TolX', 1e-6, 'TolFun', 1e-6);

    % Minimize the negative log-likelihood.
    [optParams, negLogL] = fminsearch(@(params) negLogLikelihood(params, x, y), initParams, options);

    % Convert back from log-scale.
    lambda = exp(optParams(1));
    nu = exp(optParams(2));
end

%% Local function: Negative Log-Likelihood for CMP distribution
function nll = negLogLikelihood(params, x, y)
    % Unpack parameters (they are in log scale).
    log_lambda = params(1);
    log_nu = params(2);
    lambda = exp(log_lambda);
    nu = exp(log_nu);

    totalCount = sum(y);
    
    % Compute the normalization constant Z(lambda, nu)
    Z = cmpZ(lambda, nu);
    logZ = log(Z);
    
    % Use gammaln(x+1) for log(x!) to ensure numerical stability.
    logFactorial = gammaln(x+1);
    
    % Compute log-likelihood:
    % logL = sum( y .* ( x*log(lambda) - nu*gammaln(x+1) ) ) - totalCount*log(Z)
    logL = sum( y .* ( x * log(lambda) - nu .* logFactorial ) ) - totalCount * logZ;
    
    % Return negative log-likelihood (for minimization).
    nll = -logL;
end

%% Local function: Compute normalization constant Z(lambda, nu)
function Z = cmpZ(lambda, nu)
    % Compute Z(lambda, nu) = sum_{j=0}^∞ lambda^j/(j!)^nu by summing until the term is negligible.
    maxIter = 1000;  % maximum number of terms to sum
    tol = 1e-10;     % tolerance for terminating the series sum
    Z = 0;
    j = 0;
    term = cmpTerm(j, lambda, nu);
    while j < maxIter && term > tol
        Z = Z + term;
        j = j + 1;
        term = cmpTerm(j, lambda, nu);
    end
    if j == maxIter
        warning('cmpZ: Maximum iterations reached while computing the normalization constant.');
    end
end

%% Local function: Compute a single term of the CMP normalization series
function t = cmpTerm(j, lambda, nu)
    % Compute term = lambda^j / (j!)^nu using gammaln for numerical stability.
    % log(t) = j*log(lambda) - nu * gammaln(j+1)
    t = exp( j * log(lambda) - nu * gammaln(j+1) );
end

%%
function [lambda, nu] = cmpFromMeanVariance(desiredMean, desiredVar)
    % cmpFromMeanVariance estimates the parameters of a Conway-Maxwell-Poisson
    % distribution given a desired mean and variance.
    %
    %   [lambda, nu] = cmpFromMeanVariance(desiredMean, desiredVar)
    %
    % Inputs:
    %   desiredMean - the target mean (E[X])
    %   desiredVar  - the target variance (Var[X])
    %
    % Outputs:
    %   lambda      - estimated rate parameter (λ > 0)
    %   nu          - estimated dispersion parameter (ν > 0)
    %
    % The CMP probability mass function is given by:
    %   P(X = k) = (lambda^k / (k!)^nu) / Z(lambda, nu)
    % where Z(lambda, nu) = sum_{j=0}^∞ (lambda^j/(j!)^nu)
    %
    % This function defines the error functions based on the computed moments,
    % and uses fsolve to find lambda and nu that satisfy:
    %   computedMean - desiredMean = 0   and   computedVar - desiredVar = 0
    
    % Set initial guesses. For a Poisson (ν=1), mean equals lambda.
    lambda0 = desiredMean;
    % Adjust initial nu guess based on dispersion:
    if desiredVar > desiredMean
        nu0 = 0.8; % overdispersion typically corresponds to ν < 1
    elseif desiredVar < desiredMean
        nu0 = 1.2; % underdispersion typically corresponds to ν > 1
    else
        nu0 = 1;   % equidispersion (Poisson)
    end
    
    % Combine initial guesses into a vector.
    initialGuess = [lambda0; nu0];
    
    % Options for fsolve (you can adjust tolerances as needed)
    options = optimoptions('fsolve', 'Display', 'off', ...
        'TolFun', 1e-8, 'TolX', 1e-8);
    
    % Define the system of equations as a function handle.
    fun = @(x) cmpMomentEquations(x, desiredMean, desiredVar);
    
    % Solve the system for [lambda; nu].
    sol = fsolve(fun, initialGuess, options);
    lambda = sol(1);
    nu = sol(2);
end

%% Subfunction: Define the moment equations error
function F = cmpMomentEquations(x, desiredMean, desiredVar)
    lambda = x(1);
    nu = x(2);
    [computedMean, computedVar] = cmpMoments(lambda, nu);
    F = [computedMean - desiredMean; computedVar - desiredVar];
end

%% Subfunction: Compute CMP moments by summing the series
function [mu, variance] = cmpMoments(lambda, nu)
    % Set maximum number of terms and tolerance for series truncation
    maxIter = 1000;
    tol = 1e-10;
    
    % Initialize sums for normalizing constant, first and second moments.
    Z = 0;
    sum_k = 0;
    sum_k2 = 0;
    
    k = 0;
    term = cmpTerm(k, lambda, nu);
    
    % Sum until the term is very small or max iterations reached.
    while (term > tol || k < 10) && k < maxIter
        Z = Z + term;
        sum_k = sum_k + k * term;
        sum_k2 = sum_k2 + k^2 * term;
        k = k + 1;
        term = cmpTerm(k, lambda, nu);
    end
    
    % Compute mean and variance.
    mu = sum_k / Z;
    variance = sum_k2 / Z - mu^2;
end

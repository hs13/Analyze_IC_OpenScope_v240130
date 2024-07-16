% make a random lognormal^2 connectivity matrix
% 1. start with a square matrix where elements are random and normally distributed
% (centered at 0), then make exponential distribution
% 2. multiply a row vector and a column vector that are both log-normally
% distributed
% 3. randomly select ~20% of the columns to be inhibitory (just put a negative
% sign on the weights)


% 4. find the eigenvector of this matrix -- should be log-normally
% distributed

% WHY DOES DP CORRELATE WITH SPONTANEOUS FIRING RATE?
% select ~20% of neurons as 'SP' neurons (sensory neurons, include both E and I)
% on every 'trial', select ~50% of these neurons to be stimulated
% (note, eventually we want to base these numbers on my detection dataset)
% the overall network activity should correlate with activity for high
% firing rate neurons
% this explains why DP is most prominent in high firing rate neurons
% another prediction of this is that DP neurons (and/or high firing rate
% neurons) should show trial-by-trial correlation with overall firing rate
% in the network (to avoid circularity, remove the neuron being probed from
% the network...).

% POPULATION COUPLING, SYNAPTIC INPUT STRENGTH AND SPONTANEOUS FIRING RATE
% note -- this is effectively what the choroister paper shows! population
% coupling measures how well a neuron's activity correlates with the rest of the
% population. Okun et al. showed that population coupling correlates well
% with the input amount. We need to show that population coupling
% correaltes with spontaneous firing rate

% CHECK CONNECTIVITY ASSUMPTIONS WITH OUR CCG DATA
% use CCG peak height matrix (jitter mean
% subtracted, normalized by geometric mean of rate) to check:
% i) approximate connectivity weight as CCG peak
% height, and see if they follow a log normal distribution
% ii) check whether firing rate correlates with the eigenvector of the
% peak height matrix

% AVENUE FOR FUTURE EXPLORATION: PATTERN COMPLETION IN THIS LOGNORMAL NETWORK
% Hopfield network
% specific inputs lead to specific attractor states
% which inputs? perhaps non-principal eigenvectors...?
% for these attractor states, can we demonstrate that pattern completion
% can be driven by activating the neuron with the most output connectivity?


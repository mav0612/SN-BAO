function [chain, loglikes] = mcmc3(initial, loglike, stddevs, Nsteps, frac)
% MCMC sampler using Metropolis-Hastings algorithm
% Inputs:
%   - initial: vector of initial parameter values
%   - loglike: function handle to compute log-likelihood: loglike(params)
%   - stddevs: vector of standard deviations for proposal steps
%   - Nsteps: total number of MCMC steps
%   - frac: fraction of steps to store (0 < frac <= 1)
%
% Outputs:
%   - chain: matrix of stored samples (Nsaved x Nparams)
%   - loglikes: vector of corresponding log-likelihood values

    % Initialization
    current = initial(:)';
    ndim = length(current);
    logL_current = loglike(current);
    chain = [];
    loglikes = [];
    
    % Determine number of saved steps and saving interval
    Nsave = round(frac * Nsteps);
    stride = round(1 / frac);
    chain = zeros(Nsave, ndim);
    loglikes = zeros(Nsave, 1);
    
    % Proposal covariance matrix
    proposal_cov = diag(stddevs.^2);
    
    % Cholesky decomposition for efficient sampling
    chol_cov = chol(proposal_cov, 'lower');

    save_idx = 0;
    for step = 1:Nsteps
        % Generate proposal
        proposal = current + (chol_cov * randn(ndim,1))';
        logL_proposal = loglike(proposal);
        
        % Metropolis acceptance rule
        if log(rand) < logL_proposal - logL_current
            current = proposal;
            logL_current = logL_proposal;
        end
        
        % Save every 'stride' steps
        if mod(step, stride) == 0
            save_idx = save_idx + 1;
            chain(save_idx, :) = current;
            loglikes(save_idx) = logL_current;
        end
    end
end
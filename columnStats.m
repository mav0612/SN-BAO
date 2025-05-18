function [means, stds] = columnStats(data)
%COLUMNSTATS  Compute mean and standard deviation of each column
%
%   [means, stds] = columnStats(data)
%
%   Inputs:
%     data  - an M×N numeric matrix, where each of the M rows is an
%             observation and each of the N columns is a variable.
%
%   Outputs:
%     means - 1×N row vector of column means
%     stds  - 1×N row vector of column standard deviations
%
%   By default, std uses the unbiased estimator (normalizing by M–1).
%   If you prefer the population standard deviation (normalizing by M),
%   replace the call to std with std(data,1,1).

    % Validate input
    if isempty(data)
        % No data → return empty row vectors
        means = zeros(1, size(data,2));
        stds  = zeros(1, size(data,2));
        return;
    end

    % Compute column means
    means = mean(data, 1);

    % Compute column standard deviations (unbiased, i.e. normalization by M-1)
    stds  = std(data, 0, 1);
end
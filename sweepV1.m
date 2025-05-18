function results = sweepV1(V1min, V1max, nSteps, initSeed)
% sweepD2combo012MIN  Scan V2 with fminsearch on D2combo012MIN
%
%   results = sweepD2combo012MIN(V2min, V2max, nSteps, initSeed)
%
% Inputs:
%   V2min    – minimum value of global V2
%   V2max    – maximum value of global V2
%   nSteps   – number of V2 increments (must be ≥1)
%   initSeed – 1×4 vector, initial guess for D2combo012MIN
%
% Output:
%   results  – nSteps×6 array. Each row is
%              [V2, x1, x2, x3, x4, fval]

    % Preallocate
    results = zeros(nSteps, 5);
    
    % Create sweep vector
    V1vals = linspace(V1min, V1max, nSteps);
    
    % Initialize seed
    xSeed = initSeed(:);

    % Declare global
    global V1

    for i = 1:nSteps
        % Set current V2
        V1 = 1- V1vals(i);

        % Define objective
        obj = @(x) D2combo1635(x);

        % Run fminsearch
        [xOpt, fval] = fminsearch(obj, xSeed);

        % Store: [V2, xOpt(1:4), fval]
        results(i, :) = [V1, xOpt(:)', fval];

        % Update seed
        xSeed = xOpt;
    end
end
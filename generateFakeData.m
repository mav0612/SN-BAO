function [fake_BAO, fake_SN] = generateFakeData(inputs)
% generateFakeData  Simulate BAO & SN data for a given cosmology
%
%   [fake_BAO, fake_SN] = generateFakeData(inputs)
%
%   inputs = [OM, H0, V0, V1]
%   Globals (must be defined in workspace):
%     c_km_per_s   – speed of light [km/s]
%     BAOdata      – BAO table: [z, Obs, type]
%     BAOinv       – inverse covariance matrix for BAOdata(:,2)
%     data         – SN table: columns include z_CMB (col 3), z_Hel (col 5), mu_obs (col 6)
%     covinv       – inverse covariance matrix for SN mu values
%
%   fake_BAO is an N_bao×1 vector of simulated BAO observables,
%   fake_SN  is an N_sn×1 vector of simulated SN distance moduli.

    global c_km_per_s BAOdata BAOinv data covinv

    % unpack parameters
    OM = inputs(1);
    H0 = inputs(2);
    V0 = inputs(3);
    V1 = inputs(4);

    %--- compute theory BAO ---
    rd = 1;                      % sound horizon scale (set =1)
    Pre = c_km_per_s / H0 / rd;
    tspan = [1 0.25];
    ics = [0; sqrt(2*(1 - OM - V0)); 0];
    sol = ode45(@(a,y) myODE(a,y,OM,V0,V1), tspan, ics);

    a_bao = 1 ./ (1 + BAOdata(:,1));
    etamodel = deval(sol, a_bao, 3);
    phimodel = deval(sol, a_bao, 1);
    dphimodel= deval(sol, a_bao, 2);
    Hmodel = sqrt( OM./a_bao.^3 + 0.5*(dphimodel.^2) + V0 + V1*phimodel );

    Dm = Pre * etamodel;
    Dh = Pre ./ Hmodel;
    Dv = Pre * ( BAOdata(:,1) .* (etamodel.^2) ./ Hmodel ).^(1/3);

    % assemble the theory vector in the same order/type as BAOdata(:,3)
    type = BAOdata(:,3);
    TheoBAO = zeros(size(type));
    TheoBAO(type==1) = Dh(type==1);
    TheoBAO(type==2) = Dm(type==2);
    TheoBAO(type==3) = Dv(type==3);

  %--- compute theory SN distance moduli (fixed orientations) ---
    zcmb = data(:,3)';      % make 1×N
    zhel = data(:,5)';      % make 1×N
    a_sn = 1 ./ (1 + zcmb); % 1×N
    eta_sn = deval(sol, a_sn, 3);                         % 1×N
    mu_th  = (5*log10( eta_sn .* (1 + zhel) ) + 25)';      % → N×1 column
    
    %--- invert inverse-covariances to get full covariances ---
    BAOcov = inv(BAOinv);
    SNcov  = inv(covinv);

    %--- draw random errors with the correct covariance ---
    Lbao = chol(BAOcov, 'lower');   % BAO: N_bao×N_bao
    Lsn  = chol(SNcov,  'lower');   % SN:  N_sn ×N_sn

    rand_BAO = Lbao * randn(size(BAOdata,1),1);
    rand_SN  = Lsn  * randn(size(data,1),1);

    %--- fake observed data ---
    fake_BAO = TheoBAO + rand_BAO;
    fake_SN  = mu_th   + rand_SN;
end


function dy = myODE(a,y,OM,V0,V1)
    % same background ODE as in D2comboB
    E = sqrt( OM/a^3 + 0.5*y(2)^2 + V0 + V1*y(1) );
    dy = zeros(3,1);
    dy(1) = -y(2) / (a * E);
    dy(2) = V1/(a * E) - 3*y(2)/a;
    dy(3) = -1/(a^2 * E);
end
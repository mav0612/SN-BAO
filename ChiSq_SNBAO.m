function cstot = ChiSq_SNBAO(inputs)
    % D2comboB computes the combined log-likelihood (chi-squared) from BAO and
    % supernova data for a linear potential dark energy model.

    % Global constants and data
    global c_km_per_s    % Speed of light in km/s
    global BAOdata       % BAO data array: [z, observed_distance, type]
    global BAOinv        % Inverse covariance matrix for BAO measurements
    global data          % Supernova data array: columns including redshift and mu
    global covinv        % Inverse covariance matrix for supernovae

    % Unpack input parameters
    OM = inputs(1);      % Matter density parameter \Omega_m
    H0 = inputs(2);      % Hubble constant today (km/s/Mpc)
    V0 = inputs(3);      % Linear potential parameter V0
    V1 = inputs(4);      % Linear potential slope V1

    rd = 1;  % Sound horizon scale (set to unity for normalization)

    % Parameter validity check: ensure physical parameter range
    if OM < 0 || OM + V0 > 1
        % Return large negative log-likelihood for invalid parameters
        cstot = -4000;
    else
        % Pre-factor for distance calculations: c/(H0 * rd)
        Pre = c_km_per_s / H0 / rd;

        % Set up ODE integration for the scalar field and "eta" function
        tspan = [1 0.25];  % Integrate from scale factor a=1 down to a=0.25
        % Initial conditions [phi(a=1); dphi/da(a=1); eta(a=1)]
        ics = [0; sqrt(2*(1 - OM - V0)); 0];
        % Solve ODE system defined in myODE
        sol = ode45(@(t,y) myODE(t,y,OM,V0,V1), tspan, ics);

        % ----- BAO likelihood calculation -----
        % Convert BAO redshifts to scale factor a
        aD2 = 1 ./ (1 + BAOdata(:,1));
        % Evaluate solution at each a: eta, phi, dphi/da
        etamodel  = deval(sol, aD2, 3);
        phimodel  = deval(sol, aD2, 1);
        dphimodel = deval(sol, aD2, 2);
        % Compute H(a) via Friedmann eqn including scalar field terms
        Hmodel = sqrt( OM ./ (aD2'.^3) + 0.5*(dphimodel.^2) + V0 + V1*phimodel );

        % BAO distance measures:
        Dm_theory = Pre * etamodel;                                  % comoving distance
        Dh_theory = Pre ./ Hmodel;                                   % Hubble distance
        Dv_theory = Pre * ( BAOdata(:,1)' .* (etamodel.^2) ./ Hmodel ).^(1/3); % volume distance

        % Select the appropriate observable for each BAO data point
        type = BAOdata(:,3);   % 1=Dh, 2=Dm, 3=Dv
        Theory = zeros(size(type));
        Theory(type==1) = Dh_theory(type==1);
        Theory(type==2) = Dm_theory(type==2);
        Theory(type==3) = Dv_theory(type==3);

        % Compute residuals and BAO log-likelihood: -1/2 * chi^2
        Obs = BAOdata(:,2);
        Res = Theory - Obs;
        cstot = -0.5 * (Res' * (BAOinv * Res));

        % ----- Supernovae likelihood calculation -----
        % Convert SN redshifts to scale factor a
        aData = 1 ./ (1 + data(:,3)');
        % Evaluate eta at SN a-values
        etamodel = deval(sol, aData, 3);
        % Compute model distance modulus mu = 5 log10(d_L / 10pc)
        % Note: d_L = eta * (1+z)
        mumodel = 5*log10(etamodel .* (1 + data(:,5)')) + 25;
        mudata  = data(:,6)';  % Observed distance moduli

        % Residuals for supernovae
        del = mumodel - mudata;
        Nsn = length(del);
        Onevec = diag(ones(Nsn));  % Vector of 1s used for analytic marginalization

        % Compute covariance terms
        CDC = del * (covinv * del');  % chi^2 term
        CD  = del * (covinv * Onevec);% cross term with nuisance
        C   = Onevec' * (covinv * Onevec); % normalization

        % Marginalized log-likelihood including nuisance parameter
        cs = CDC - CD^2/C + log(C / (2*pi));
        % Add SN contribution to total log-likelihood
        cstot = cstot - 0.5 * cs;
    end
end

% ODE system defining phi(a), dphi/da, and eta(a)
function dy = myODE(a, y, OM, V0, V1)
    dy = zeros(3,1);
    % Denominator: H(a) = sqrt(OM/a^3 + scalar field energy)
    denom = sqrt(OM/a^3 + 0.5*y(2)^2 + V0 + V1*y(1));
    % d(phi)/da
    dy(1) = -y(2) / (a * denom);
    % d(dphi/da)/da
    dy(2) = V1/(a * denom) - 3*y(2)/a;
    % d(eta)/da = -1/(a^2 H(a))
    dy(3) = -1 / (a^2 * denom);
end

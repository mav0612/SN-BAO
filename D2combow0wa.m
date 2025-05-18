function cstot = D2combow0wa(inputs)
    global c_km_per_s
    global BAOdata
    global BAOinv
    global data;
    global covinv;
    OM = inputs(1);
    H0 = inputs(2);
    w0 = inputs(3);
    wa = inputs(4);
    rd = 1;
    %lh = hcalc(H0);
if OM < 0 
    cstot = -200;
else
    Pre =  (c_km_per_s/H0/rd);
    tspan = [1 0.25];
    ics = [0];
    sol = ode45(@(t,y) myODE(t,y,OM,w0,wa),tspan,ics);

    aD2 = 1 ./ (1 + BAOdata(:,1));
    etamodel = deval(sol,aD2,1);
    Hmodel = (OM ./ (aD2' .^ 3) + (1 - OM) ./ (aD2' .^ (3*(1 + w0 + wa))) .* exp(3*wa*(aD2'-1))) .^(1/2);
  
    Dm_theory = Pre * etamodel;
    Dh_theory = Pre ./ Hmodel; 
    Dv_theory = Pre * ( BAOdata(:,1)' .* (etamodel .^2 ) ./ Hmodel) .^ (1/3);
    
    type = BAOdata(:,3);
    Theory = zeros(size(type));       
    mask1 = (type==1);
    mask2 = (type==2);
    mask3 = (type==3);

    Theory(mask1) = Dh_theory(mask1);
    Theory(mask2) = Dm_theory(mask2);
    Theory(mask3) = Dv_theory(mask3);

    Obs = BAOdata(:,2);
    Res = Theory - Obs;
    cstot = -1/2*( Res' * (BAOinv * Res ));
    
    %Supernovae

    aData = 1 ./ (1 + data(:,3)');
    etamodel = deval(sol,aData,1);
    mumodel = 5*log10(etamodel .* (1 + data(:,5)')) + 25;
    mudata = data(:,6)';
    Onevec = diag(ones(1829));
    del = mumodel - mudata;
    CDC = del * (covinv * del');
    CD = del*( covinv * Onevec);
    C = Onevec'*( covinv * Onevec);
    cs = CDC - CD^2/C + log(C/2/3.141592);
    cstot = cstot -1/2*(cs - 1654);
   
end
end

function dy = myODE(a,y,OM,w0,wa)
  dy = zeros(1,1);
  dy(1) = -1/a^2/sqrt(OM/a^3 + (1 - OM)/a^(3*(1 + w0 + wa))*exp(3*wa*(a-1)));
end
function cstot = D2comboF(inputs)
    global c_km_per_s
    global fBAO;
    global BAOinv;
    global fSN;
    global covinv;
    OM = inputs(1);
    H0 = inputs(2);
    V0 = inputs(3);
    V1 = inputs(4);
    rd = 1;
    %lh = hcalc(H0);
if OM < 0 ||  OM + V0 > 1 
    cstot = -4000;
else
    Pre =  (c_km_per_s/H0/rd);
    tspan = [1 0.25];
    ics = [0;sqrt(2*(1 - OM - V0));0];
    sol = ode45(@(t,y) myODE(t,y,OM,V0,V1),tspan,ics);

    aD2 = 1 ./ (1 + fBAO(:,1));
    etamodel = deval(sol,aD2,3);
    phimodel = deval(sol,aD2,1);
    dphimodel = deval(sol,aD2,2);
    Hmodel = (OM ./ (aD2' .^ 3) + 1/2*(dphimodel .^ 2) + V0 + V1* phimodel) .^(1/2);
    Dm_theory = Pre * etamodel;
    Dh_theory = Pre ./ Hmodel; 
    Dv_theory = Pre * ( fBAO(:,1)' .* (etamodel .^2 ) ./ Hmodel) .^ (1/3);
    
    type = fBAO(:,3);
    Theory = zeros(size(type));       
    mask1 = (type==1);
    mask2 = (type==2);
    mask3 = (type==3);

    Theory(mask1) = Dh_theory(mask1);
    Theory(mask2) = Dm_theory(mask2);
    Theory(mask3) = Dv_theory(mask3);

    Obs = fBAO(:,2);
    Res = Theory - Obs;
    cstot = -1/2*( Res' * (BAOinv * Res ));
    
    %Supernovae

    aData = 1 ./ (1 + fSN(:,3)');
    etamodel = deval(sol,aData,3);
    mumodel = 5*log10(etamodel .* (1 + fSN(:,5)')) + 25;
    mudata = fSN(:,6)';
    Onevec = diag(ones(1829));
    del = mumodel - mudata;
    CDC = del * (covinv * del');
    CD = del*( covinv * Onevec);
    C = Onevec'*( covinv * Onevec);
    cs = CDC - CD^2/C + log(C/2/3.141592);
    cstot = cstot -1/2*(cs);
    
end
end

function dy = myODE(a,y,OM,V0,V1)
  dy = zeros(3,1);
  dy(1) = -y(2)/a/sqrt(OM/a^3 + 1/2*y(2)^2 + V0 + V1*y(1));
  dy(2) = V1/a/sqrt(OM/a^3 + 1/2*y(2)^2 + V0 + V1*y(1)) - 3*y(2)/a;
  dy(3) = -1/a^2/sqrt(OM/a^3 + 1/2*y(2)^2 + V0 + V1*y(1));
end
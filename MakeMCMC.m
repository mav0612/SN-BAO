loglike=@(m)ChiSq_SNBAO(m);
minit=[0.304786029303   9967.291121505626 0.6 1.5];
m = mcmc3(minit,loglike,[0.03 50 0.06 0.1],200000,0.05);
m(1:1000,:)=[]; %crop drift





fakedat = zeros(16,8);
for i = 1:16
    [fakeBAO,fakeSN] = generateFakeData([0.31,10054,0.69,0]);
    fBAO = [BAOdata(:,1),fakeBAO,BAOdata(:,3)];
    fSN = [data(:,1:5),fakeSN,data(:,7)];
    testmcmc3;
    [a,b] = columnStats(m);
    fakedat(i,:) = [a,b];
    columnStats(fakedat(1:i,:))
end
save fake.mat

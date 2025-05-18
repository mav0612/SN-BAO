global data
global covinv
global c_km_per_s
global BAOdata
global BAOinv

%Load DES 5yr supernova data and covariance matrix
%Source: https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT

data = readmatrix('SN2.txt');
covvec = readmatrix('SNcov.txt', 'NumHeaderLines', 1);
cd = diag(data(:,7));
covmat = reshape(covvec,1829,1829);
covmat = covmat + cd * cd;
covinv = inv(covmat);

%Load DESI DR2 data and covariance matrix
%Source: https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2

c_km_per_s = 299792.458;
BAOdata = readmatrix('BAOdat.txt');
BAOcov = readmatrix('BAOcov.txt');
BAOinv = inv(BAOcov);

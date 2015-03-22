%
% cris_test2 -- test two versions of iasi2cris
% 

% set paths to asl libs
addpath /asl/packages/ccast/source
addpath /asl/packages/iasi_decon

% get IASI user-grid parameters
iasi = iasi_params;
ifrq1 = iasi.freq;

% load some AIRS/IASI matchups
afile = 'data/airs_iasi_sample2.mat';
d1 = load(afile);
irad1 = d1.ri1000(:, 1:1000);
clear d1

% translation options
opt1 = struct;
opt1.hapod = 0;     
opt1.nguard = 2;
opt1.resmode = 'hires2';

% do the new translation
tic
[crad1, cfrq1] = iasi2cris(irad1, ifrq1, opt1);
toc

% do the old translation
tic
[crad2, cfrq2] = iasi2cris_old(irad1, ifrq1, opt1);
toc

isequal(crad1, crad2)



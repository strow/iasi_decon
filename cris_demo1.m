%
% cris_demo1 - IASI to CrIS translation demo
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
irad1 = d1.ri1000(:, 1:800);
clear d1

% translation options
opt1 = struct;
opt1.hapod = 0;     
opt1.nguard = 0;
opt1.resmode = 'hires2';

% do the translation
tic
[crad1, cfrq1] = iasi2cris(irad1, ifrq1, opt1);

% timing report
[m,n] = size(irad1);
fprintf(1, 'translated %d obs in %.2f seconds\n', n, toc)

% translation with guard chans
opt1.nguard = 4;
[crad2, cfrq2] = iasi2cris(irad1, ifrq1, opt1);

% convert to brightness temps
ibt1 = real(rad2bt(ifrq1, irad1));
cbt1 = real(rad2bt(cfrq1, crad1));
cbt2 = real(rad2bt(cfrq2, crad2));

% plot IASI and IASI to CrIS 
iobs = 800;
figure(1); clf
plot(ifrq1, ibt1(:,iobs), cfrq1, cbt1(:,iobs))
title(sprintf('IASI and IASI CrIS, obs %d', iobs))
legend('IASI', 'IASI CrIS', 'location', 'north')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on

% IASI to CrIS with and without guard chans
figure(2); clf
plot(cfrq2, cbt2(:,iobs), 'r', cfrq1, cbt1(:,iobs), 'g')
title(sprintf('IASI CrIS with guard chans, obs %d', iobs))
legend('guard chans', 'no guard cnans', 'location', 'south')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on



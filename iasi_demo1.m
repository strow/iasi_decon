%
% iasi_demo1 - IASI to AIRS translation demo
% 

% set paths to asl libs
addpath /asl/matlib/h4tools
addpath /asl/packages/iasi_decon

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% get IASI user-grid parameters
iasi = iasi_params;
ifrq1 = iasi.freq;

% get AIRS channel frequencies
afile = 'data/airs-iasi-sample.mat';
d1 = load(afile);
afrq1 = d1.fairs;

% load some AIRS/IASI matchups
afile = 'data/airs_iasi_sample2.mat';
d1 = load(afile);
arad1 = d1.ra1000(:, 1:400);
irad1 = d1.ri1000(:, 1:400);
clear d1

% do the translation
tic
[arad2, afrq2] = iasi2airs(irad1, ifrq1, sfile, afrq1);

% timing report
[m,n] = size(irad1);
fprintf(1, 'translated %d obs in %.2f seconds\n', n, toc)

% convert to brightness temps
ibt1 = real(rad2bt(ifrq1, irad1));
abt1 = real(rad2bt(afrq1, arad1));
abt2 = real(rad2bt(afrq2, arad2));

% select a sample obs
iobs = 21;

% sample obs overview
figure(1); clf
plot(ifrq1, ibt1(:,iobs), afrq1, abt1(:,iobs), afrq2, abt2(:,iobs))
ax = axis; ax(3) = 140; ax(4) = 300; axis(ax);
title(sprintf('IASI, true AIRS, and IASI AIRS, obs %d', iobs))
legend('IASI', 'true AIRS', 'IASI AIRS', 'location', 'south')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on

% sample obs residual
figure(2); clf
plot(afrq1, abt2(:, iobs) - abt1(:, iobs))
ax = axis; ax(3) = -8; ax(4) = 8; axis(ax);
title(sprintf('IASI AIRS minus true AIRS, obs %d', iobs))
xlabel('wavenumber')
ylabel('brightness temp diff')
grid on; zoom on

% mean all obs residual
figure(3); clf
plot(afrq1, mean(abt2 - abt1, 2))
ax = axis; ax(3) = -8; ax(4) = 8; axis(ax);
title('IASI AIRS minus true AIRS mean')
xlabel('wavenumber')
ylabel('brightness temp diff')
grid on; zoom on

% std all obs residual
figure(4); clf
plot(afrq1, std(abt2 - abt1, 0, 2))
ax = axis; ax(3) = -2; ax(4) = 10; axis(ax);
title('IASI AIRS minus true AIRS std')
xlabel('wavenumber')
ylabel('brightness temp diff')
grid on; zoom on


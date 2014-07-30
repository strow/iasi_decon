%
% demo1 - IASI to AIRS translation demo
% 

% set paths to standard libs
addpath /asl/matlib/h4tools

% specify an SRF tabulation file
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% specify and load the IASI data
afile = 'airs-iasi-sample.mat';
d1 = load(afile);
ifrq1 = d1.fiasi;
irad1 = d1.radiasi;

% sort AIRS data by channel frequency
[afrq1, ix] = sort(d1.fairs);
arad1 = d1.radairs(ix, :);

% do the translation
[arad2, afrq2] = iasi2airs(irad1, ifrq1, sfile, afrq1);

% convert to brightness temps
ibt1 = real(rad2bt(ifrq1, irad1));
abt1 = real(rad2bt(afrq1, arad1));
abt2 = real(rad2bt(afrq2, arad2));

% select an obs
iobs = 3;

% overview plot
figure(1); clf
plot(ifrq1, ibt1(:,iobs), afrq1, abt1(:,iobs), afrq2, abt2(:,iobs))
title(sprintf('IASI, true AIRS, and IASI AIRS, obs %d', iobs))
legend('IASI', 'true AIRS', 'IASI AIRS', 'location', 'south')
xlabel('wavenumber')
ylabel('brightness temp')
grid on; zoom on

% residual plot
[j1, j2] = seq_match(afrq1, afrq2);
afrq1 = afrq1(j1);
afrq2 = afrq2(j2);
abt1 = abt1(j1, iobs);
abt2 = abt2(j2, iobs);
figure(2); clf
plot(afrq1, abt2 - abt1)
title(sprintf('IASI AIRS minus true AIRS, obs %d', iobs))
xlabel('wavenumber')
ylabel('brightness temp diff')
grid on; zoom on

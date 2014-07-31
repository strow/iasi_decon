% 
% iasi_test1 -- compare IASI to AIRS with true AIRS
% 
% reference truth: start with kcarta radiances, convolve these to 
% AIRS channel radiances, and call the result “true AIRS”.
% 
% deconvolution: start with kcarta radiances, convolve to IASI
% channel radiances (“true IASI”), deconvolve to an intermediate
% grid, e.g. 0.05 1/cm spacing, and reconvolve to the AIRS grid
% (“IASI AIRS”).  Then compare IASI AIRS vs true AIRS over the
% set of kcarta test radiances
%
% comparison tests: (1) IASI deconvolution to an intermediate 
% grid, followed by an AIRS convolution, (2) IASI interpolation 
% to an intermediate grid, followed by an AIRS convolution, and 
% (3) simple interpolation from IASI to AIRS
%

%-----------------
% test parameters
%-----------------

addpath /asl/matlib/h4tools

dvb = 0.1;       % deconvolution frequency step
dvk = 0.0025;    % kcarta frequency spacing
fig = 'png';     % plot type

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% AIRS 1C channel frequencies
cfreq = load('data/freq2645.txt');  

% specify an AIRS SRF tabulation
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% build the AIRS convolution matrix
[sconv, sfreq, tfreq] = mksconv1(sfile, cfreq, dvk);

%-----------------------------------
% true IASI and true AIRS radiances
%-----------------------------------

% loop on kcarta files
rad1 = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % kcarta to IASI channel radiances
  [rtmp, ftmp] = kc2iasi(rkc, vkc);
  rad1 = [rad1, rtmp];

  % kcarta to AIRS channel radiances
  [ix, jx] = seq_match(sfreq, vkc);
  rtmp = zeros(length(sfreq), 1);
  rtmp(ix) = rkc(jx);
  rtmp = sconv * rtmp;
  rad2 = [rad2, rtmp];

  fprintf(1, '.');
end
fprintf(1, '\n')
frq1 = ftmp(:);     % from kc2iasi
frq2 = tfreq(:);    % from mksconv2
clear d1 vkc rkc

%-------------------------------
% convolution and interpolation
%-------------------------------

% deconvolve the AIRS radiances
[rad3, frq3] = iasi_decon(rad1, frq1, dvb);

% IASI to AIRS via deconvolution and AIRS convolution
[rad4, frq4] = iasi2airs(rad1, frq1, sfile, cfreq, dvb);

% IASI direct interpolation to AIRS
rad5 = interp1(frq1, rad1, frq2, 'spline', 'extrap');

% IASI to AIRS via interpolaton and AIRS convolution
[rad6, frq6] = iasi2airsX(rad1, frq1, sfile, cfreq, dvb);

%-----------------
% stats and plots
%-----------------

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % true IASI
bt2 = real(rad2bt(frq2, rad2));   % true AIRS
bt3 = real(rad2bt(frq3, rad3));   % deconvolved IASI
bt4 = real(rad2bt(frq4, rad4));   % IASI decon AIRS conv
bt5 = real(rad2bt(frq2, rad5));   % IASI interp to AIRS
bt6 = real(rad2bt(frq6, rad6));   % IASI interp AIRS conv

% IASI and AIRS overview
figure(1); clf; j = 1; 
plot(frq1, bt1(:,j), frq2, bt2(:,j), frq3, bt3(:,j), frq4, bt4(:,j))
axis([600, 2800, 140, 320])
legend('true IASI', 'true AIRS', 'IASI decon', 'IASI AIRS', ...
       'location', 'south')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('IASI and AIRS profile %d', j));
grid on; zoom on
saveas(gcf, 'test1_fig_1', fig)

% IASI decon AIRS conv minus true AIRS
figure(2); clf
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2800, -1.0, 1.0])
xlabel('wavenumber'); ylabel('dBT')
title('IASI AIRS minus true AIRS mean');
grid on; zoom on
saveas(gcf, 'test1_fig_2', fig)

% IASI interp to AIRS minus true AIRS
figure(3); clf
plot(frq2, mean(bt5 - bt2, 2))
xlabel('wavenumber'); ylabel('dBT')
title('interpolated AIRS minus true AIRS mean');
grid on; zoom on
saveas(gcf, 'test1_fig_3', fig)

% IASI interp AIRS conv minus true AIRS 
figure(4); clf
plot(frq2, mean(bt6 - bt2, 2))
xlabel('wavenumber'); ylabel('dBT')
title('IASI interp AIRS conv  minus true AIRS mean');
grid on; zoom on
saveas(gcf, 'test1_fig_4', fig)


% 
% airs_test1 -- compare IASI AIRS with true AIRS
% 
% reference truth: start with kcarta radiances, convolve these 
% to AIRS channel radiances, and call the result "true AIRS".
% 
% deconvolution: start with kcarta radiances, convolve to IASI
% channel radiances ("true IASI"), deconvolve to an intermediate
% grid, typically at a 0.1 cm-1 spacing, and reconvolve to AIRS
% ("IASI AIRS").  Compare IASI AIRS vs true AIRS.
%

%-----------------
% test parameters
%-----------------

addpath /asl/matlib/h4tools

dvb = 0.1;       % deconvolution frequency step
dvk = 0.0025;    % kcarta frequency spacing

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% AIRS 1C channel frequencies
cfreq = load('data/freq2645.txt');  

% specify an AIRS SRF tabulation
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';

% build the AIRS convolution matrix
[sconv, sfreq, tfreq] = mksconv1(sfile, cfreq, dvk);

%-----------------------------
% get true IASI and true AIRS
%-----------------------------

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
frq2 = tfreq(:);    % from mksconv1
clear d1 vkc rkc

%------------------------
% transform IASI to AIRS
%------------------------

% deconvolve the IASI radiances
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
rtmp = bandpass(frq1, rad1, v1+5, v2-5, 5);
[rad3, frq3] = iasi_decon(rtmp, frq1, dvb);

% call the IASI to AIRS user app
[rad4, frq4] = iasi2airs(rad1, frq1, sfile, cfreq, dvb);

%-----------------
% stats and plots
%-----------------

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % true IASI
bt2 = real(rad2bt(frq2, rad2));   % true AIRS
bt3 = real(rad2bt(frq3, rad3));   % deconvolved IASI
bt4 = real(rad2bt(frq4, rad4));   % IASI decon AIRS conv

% IASI and AIRS spectra
figure(1); clf; j = 1; 
plot(frq1, bt1(:,j), frq2, bt2(:,j), frq3, bt3(:,j), frq4, bt4(:,j))
axis([600, 2800, 140, 320])
legend('true IASI', 'true AIRS', 'IASI decon', 'IASI AIRS', ...
       'location', 'south')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('IASI and AIRS profile %d', j));
grid on; zoom on

% IASI AIRS minus true AIRS mean
figure(2); clf
subplot(2,1,1)
plot(frq2, mean(bt4 - bt2, 2))
axis([600, 2800, -0.3, 0.3])
xlabel('wavenumber'); ylabel('dBT')
title('IASI AIRS minus true AIRS mean');
grid on; zoom on

% IASI AIRS minus true AIRS std
subplot(2,1,2)
plot(frq2, std(bt4 - bt2, 0, 2))
axis([600, 2800, -0.2, 0.2])
xlabel('wavenumber'); ylabel('dBT')
title('IASI AIRS minus true AIRS std');
grid on; zoom on


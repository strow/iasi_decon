% 
% cris_test1 -- compare IASI CrIS with true CrIS
% 
% reference truth: start with kcarta radiances, convolve to the 
% CrIS user grid, and call the result "true CrIS".
% 
% deconvolution: start with kcarta radiances, convolve to IASI
% radiances ("true IASI"), and translate to CrIS with iasi2cris, 
% to get "IASI CrIS".  Compare IASI CrIS to true CrIS.
%

%-----------------
% test parameters
%-----------------

addpath /asl/packages/ccast/source
addpath /asl/packages/airs_decon/test
addpath /asl/packages/airs_decon/source

% test params
band = 'LW';
hapod = 0;

% opts for inst_params and iasi2cris
opt1 = struct;
opt1.hapod = hapod;
opt1.resmode = 'hires2';

% kcarta test data
kcdir = '/home/motteler/cris/sergio/JUNK2012/';
flist =  dir(fullfile(kcdir, 'convolved_kcart*.mat'));

% get CrIS inst and user params
wlaser = 773.1301;  % nominal value
[inst, user] = inst_params(band, wlaser, opt1);

%-----------------------------
% get true IASI and true CrIS
%-----------------------------

% loop on kcarta files
rad1 = []; rad2 = [];
for i = 1 : length(flist)
  d1 = load(fullfile(kcdir, flist(i).name));
  vkc = d1.w(:); rkc = d1.r(:);

  % kcarta to IASI channel radiances
  [rtmp, ftmp] = kc2iasi(rkc, vkc);
  rad1 = [rad1, rtmp];
  frq1 = ftmp(:);

  % kcarta to CrIS channel radiances 
  [rtmp, ftmp] = kc2cris(user, rkc, vkc);
  rad2 = [rad2, rtmp];
  frq2 = ftmp(:);

  fprintf(1, '.');
end
fprintf(1, '\n')
clear d1 vkc rkc

%------------------------
% transform IASI to CrIS
%------------------------

[rad4, frq4] = iasi2cris(rad1, frq1, opt1);

% option to apodize true CrIS
if hapod
  rad2 = hamm_app(rad2);
end

%-----------------
% stats and plots
%-----------------

% take radiances to brightness temps
bt1 = real(rad2bt(frq1, rad1));   % true IASI
bt2 = real(rad2bt(frq2, rad2));   % true CrIS
bt4 = real(rad2bt(frq4, rad4));   % IASI CrIS

% plot parameters
[i2, i4] = seq_match(frq2, frq4); 
pv1 = min(frq2(i2)) - 10; 
pv2 = max(frq2(i2)) + 10;

% IASI and CrIS spectra
figure(1); clf; j = 1; 
plot(frq1, bt1(:,j), frq2, bt2(:,j), frq4, bt4(:,j))
ax(1)=pv1; ax(2)=pv2; ax(3)=200; ax(4)=300; axis(ax)
legend('true IASI', 'true CrIS', 'IASI CrIS', ...
       'location', 'southeast')
xlabel('wavenumber'); ylabel('brighness temp')
title(sprintf('IASI and CrIS %s profile %d', band, j));
grid on; zoom on
% saveas(gcf, sprintf('iasi_cris_spec_%s', band), 'png')

% IASI CrIS minus true CrIS mean
figure(2); clf
subplot(2,1,1)
[i2, i4] = seq_match(frq2, frq4);
plot(frq2(i2), mean(bt4(i4,:) - bt2(i2,:), 2))
ax = axis; ax(1)=pv1; ax(2)=pv2; axis(ax);
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('IASI CrIS minus true CrIS %s mean', band));
grid on; zoom on

% IASI CrIS minus true CrIS std
subplot(2,1,2)
plot(frq2(i2), std(bt4(i4,:) - bt2(i2,:), 0, 2))
ax = axis; ax(1)=pv1; ax(2)=pv2; axis(ax);
xlabel('wavenumber'); ylabel('dBT')
title(sprintf('IASI CrIS minus true CrIS %s std', band));
grid on; zoom on
% saveas(gcf, sprintf('iasi_cris_diff_%s', band), 'png')


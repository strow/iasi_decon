%
% nedn_test1 -- noise estimate for the IASI to CrIS translation
%

addpath /asl/packages/ccast/source

% get a sample ccast CrIS NEdN estimate
d = load('/asl/data/cris/ccast/sdr60_hr/2016/018/SDR_d20160118_t0801033.mat');
nedn_ccast = [squeeze(d.nLW(:,5,1)); squeeze(d.nMW(:,5,1)); squeeze(d.nSW(:,5,1))];
freq_ccast = [d.vLW; d.vMW; d.vSW];
clear d

% load Larrabee's IASI info
load iasi_noise
nedn_carm = 1e3 * sqrt(iasi_noise_carmine);
nedn_cnes = 1e3 * sqrt(iasi_noise_cnes);

% test param
nrad = 200;          % simulated radiance obs per set
nset = 80;           % number of simulated radiance sets

% IASI params
v1 = 645;            % iasi band low
v2 = 2760;           % iasi band high
dv = 0.25;           % iasi dv
nchan = round((v2 - v1) / dv) + 1;
freq_iasi = v1 + (0:nchan-1)' * dv;

% get IASI expected ICT radiance
r_280K = bt2rad(freq_iasi, 280) * ones(1, nrad);

% loop on radiance sets
for i = 1 : nset

  % add noise scaled by the IASI NEdN spec
  r_iasi = r_280K + randn(nchan, nrad) .* (nedn_carm * ones(1, nrad));

  % IASI to CrIS translation
  opt1 = struct;
  opt1.resmode = 'hires2';
  [r_cris, freq_cris] = iasi2cris(r_iasi, freq_iasi, opt1);
  r_cris = real(r_cris);

  % initialize tables on the first iteration
  if i == 1
    tab_iasi = zeros(length(freq_iasi), nset);
    tab_cris = zeros(length(freq_cris), nset);
  end

  % measure and save the simulated noise (as a check)
  tab_iasi(:, i) = std(r_iasi, 0, 2);

  % measure and save the translated simulated noise
  tab_cris(:, i) = std(r_cris, 0, 2);

  fprintf(1, '.')
end
fprintf(1, '\n')

% take means over the std sets
nedn_iasi = mean(tab_iasi, 2);
nedn_cris = mean(tab_cris, 2);

% plot the results
plot(freq_iasi, nedn_carm, freq_iasi, nedn_iasi, ...
     freq_cris, nedn_cris, freq_ccast, nedn_ccast);

axis([600, 2600, 0, 0.4])
title('IASI to CrIS NEdN estimates')
legend('IASI Carmine', 'IASI simulated', 'IASI to CrIS', 'CrIS ccast')
xlabel('wavenumber')
ylabel('NEdN')
grid on; zoom on

% saveas(gcf, 'iasi2cris_nedn', 'png')
% save iasi2cris_nedn ...
%      freq_iasi nedn_carm nedn_iasi freq_cris nedn_cris freq_ccast nedn_ccast


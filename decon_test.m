
addpath /home/motteler/cris/ccast/motmsc/utils
addpath /home/motteler/cris/ccast/source

load airs-iasi-sample.mat
ifrq = fiasi;
irad = radiasi(:, 3);

[brad, bfrq] = iasi_decon(irad, ifrq);

ibt = real(rad2bt(ifrq, irad));
bbt = real(rad2bt(bfrq, brad));

figure(1); clf
plot(ifrq, ibt, bfrq, bbt)
xlabel('wavenumber')
ylabel('degrees K')
legend('iasi', 'decon')
grid on; zoom on


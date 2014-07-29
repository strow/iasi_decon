
addpath /home/motteler/cris/ccast/motmsc/utils
addpath /home/motteler/cris/ccast/source

load airs-iasi-sample.mat
[brad2, bfrq2] = iasi_decon2(radiasi, fiasi);

for i = 1 : 9

  ifrq = fiasi;
  irad = radiasi(:, i);

  [brad, bfrq] = iasi_decon(irad, ifrq);

  isequal(brad, brad2(:, i))
  isequal(bfrq, bfrq2)

end

return

ibt = real(rad2bt(ifrq, irad));
bbt = real(rad2bt(bfrq, brad));

figure(1); clf
plot(ifrq, ibt, bfrq, bbt)
xlabel('wavenumber')
ylabel('degrees K')
legend('iasi', 'decon')
grid on; zoom on


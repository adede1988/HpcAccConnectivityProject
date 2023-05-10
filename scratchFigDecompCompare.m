

figure
subplot 221
imagesc(squeeze(mean(pow,3))')
caxis([-15,15])
xticks([100:100:900])
xticklabels(mulTim([100:100:900]))
ylabel('trials')
xlabel('time (ms)')
title('multi-taper HFB all trials')

subplot 223
plot(squeeze(mean(pow, [2,3])))
xticks([100:100:900])
xticklabels(mulTim([100:100:900]))
xlim([1, size(pow,1)])
ylabel('mean z-score')
xlabel('time (ms)')
title('multi-taper HFB mean')


subplot 222
imagesc(squeeze(mean(powOG,3))')
caxis([-15,15])
xticks([500:500:4501])
xticklabels(chanDat.enctim([500:500:4501]))
ylabel('trials')
xlabel('time (ms)')
title('wavelet decomposition HFB trials')

subplot 224
plot(squeeze(mean(powOG, [2,3])))
xticks([500:500:4501])
xticklabels(chanDat.enctim([500:500:4501]))
xlim([1, size(powOG,1)])
ylabel('mean z-score')
xlabel('time (ms)')
title('wavelet decomposition HFB mean')


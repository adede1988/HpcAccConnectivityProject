



%lines 1-124 of publicationFigs.m have already been run: 


postFigPath = 'G:\My Drive\Johnson\CNS2024/';
test = load([headFilesHFB(12).folder '/' headFilesHFB(12).name]).statInfo;

reactChan = load([datPre 'CHANDAT\finished\chanDat_IR84_034.mat']).chanDat; 
testSub = test.hitChi(cellfun(@(x) strcmp(x, 'IR84'), test.hitSub));
nonReact = load([datPre 'CHANDAT\finished\chanDat_IR84_035.mat']).chanDat;

set(0,'defaultfigurewindowstyle', 'docked')


for ii = 1:71
figure
subplot 211
plot(reactChan.enctim, reactChan.enc(:,ii), 'linewidth', 1)
xlim([-450, 1250])
title(ii)
subplot 212 
HFB = bandpass(reactChan.enc(:,ii), [70,150], 1000).*4;
plot(reactChan.enctim, HFB)
hold on 
HFB = abs(hilbert(HFB)); 
plot(reactChan.enctim, HFB, 'linewidth', 3)
xlim([-450, 1250])
end


figure
ii = 56;
subplot 211
plot(reactChan.enctim, reactChan.enc(:,ii), 'linewidth', 1)
xlim([-450, 1250])
title(ii)
subplot 212 
plot(reactChan.enctim, reactChan.enc(:,ii), 'linewidth', 1)
xlim([150, 600])
ylim([-20, 60])
hold on 
HFB = bandpass(reactChan.enc(:,ii), [70,150], 1000);
HFB = HFB.*abs(hilbert(HFB).*1.5); 
plot(reactChan.enctim, HFB)
plot(reactChan.enctim, abs(hilbert(HFB)), 'linewidth', 3)

figure
ii = 54;
subplot 211
plot(reactChan.enctim, reactChan.enc(:,ii), 'linewidth', 1)
xlim([-450, 1250])
title(ii)
subplot 212 
plot(reactChan.enctim, reactChan.enc(:,ii), 'linewidth', 1)
xlim([150, 600])
ylim([-20, 60])
hold on 
HFB = bandpass(reactChan.enc(:,ii), [70,150], 1000);
HFB = HFB.*abs(hilbert(HFB).*1.5); 
plot(reactChan.enctim, HFB)
plot(reactChan.enctim, abs(hilbert(HFB)), 'linewidth', 3)
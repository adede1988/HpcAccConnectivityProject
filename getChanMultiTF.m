function [trialTF, outTim, frexOut] = ...
    getChanMultiTF(trialDat, frex, srate, times, timeStepms)

%trialDat:            timepoints X trials data points
frex = linspace(70,150,10); 

% highfrex = linspace(70, 150, 81); 
% highnumfrex = length(highfrex); 
% highstds = logspace(log10(10),log10(20),highnumfrex)./(2*pi*highfrex);
% [trialTF] = getChanTrialTF(trialDat, highfrex, highnumfrex, highstds, srate); 
% pow = abs(trialTF).^2; 
% trialTF = pow; 
% outTim = times(1:timeStepms:end); 
% trialTF = trialTF(1:5:end, :, :);



% pow = arrayfun(@(x) myChanZscore(pow(:,:,x), ...
%     [find(times>=-450,1), find(times>=-50,1)] ), 1:size(pow,3), 'UniformOutput',false ); %z-score
%   
% pow = cell2mat(pow); %organize
% pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
% pow = mean(pow,3);
% pow = pow(1:5:end,:);
% 
n_data               = size(trialDat,1); 
n_conv_pow2          = 2^(nextpow2(n_data));



[spectrum, ~, frexOut, outTim] = ft_specest_mtmconvol(trialDat', ...
        times/srate, 'freqoi', frex, ...
        'timeoi', times(1:timeStepms:end)/srate,...
        'timwin', 5./frex, ...
        'taper', 'dpss', ...
        'tapsmofrq', 0.4 *frex, ...
        'pad', n_conv_pow2 / srate);
    outTim = outTim * 1000; 
    powft = abs(spectrum).^2;
    powft = squeeze(mean(powft, 1, 'omitnan')); 
    trialTF = permute(powft, [3,1,2]);
%     powft = arrayfun(@(x) myChanZscore(powft(:,:,x), ...
%     [find(outTim>=-450,1), find(outTim>=-50,1)] ), 1:size(powft,3), 'UniformOutput',false ); %z-score
%     powft = cat(3, powft{:});
%     powft = mean(powft,3);

% 
% 
% figure('visible', false)
% lat1 = gausLat(squeeze(pow), ...
%                                outTim,...
%                                ones(size(pow,2))*1800);
% [~, order] = sort(lat1); 
% subplot 221
% imagesc(pow(:, order)')
% caxis([-10, 15])
% 
% lat2 = gausLat(squeeze(powft), ...
%                                outTim,...
%                                ones(size(powft,2))*1800);
% [~, order] = sort(lat2); 
% subplot 222
% imagesc(powft(:, order)')
% caxis([-10, 15])
% 
% w = 100; 
% 
% LatIdx = arrayfun(@(x) find(x==outTim), lat1); 
% LatIdx(LatIdx>length(outTim)-w) = length(outTim) - w; 
% LatIdx(LatIdx<=w) =  w+1; 
% vals1 = pow; %just look at one channel 
% peaks1 = arrayfun(@(x) vals1(LatIdx(x)-w:...
%                              LatIdx(x)+w, x), [1:length(LatIdx)], ...
%                              'uniformoutput', false);
% peaks1 = cell2mat(peaks1);
% 
% LatIdx = arrayfun(@(x) find(x==outTim), lat2); 
% LatIdx(LatIdx>length(outTim)-w) = length(outTim) - w; 
% LatIdx(LatIdx<=w) =  w+1; 
% vals2 = powft; %just look at one channel 
% peaks2 = arrayfun(@(x) vals2(LatIdx(x)-w:...
%                              LatIdx(x)+w, x), [1:length(LatIdx)], ...
%                              'uniformoutput', false);
% peaks2 = cell2mat(peaks2);
% 
% subplot 223
% plot([-500:5:500], mean(peaks1,2), 'color', 'blue', 'linewidth', 2)
% hold on 
% plot([-500:5:500], mean(peaks2,2), 'color', 'red', 'linewidth', 2)
% 
% subplot 224
% scatter(lat1, lat2)
% export_fig(['G:\My Drive\Johnson\MTL_PFC_networkFigs\HFB_method_testing\testFig' num2str(ii) '.jpg'], '-r300')

%% old version 


% 
% 
% %initialize indexing
% timewinms = 300; %time window in milliseconds
% % timeStepms = 25; %steps in terms of milliseconds independent of sampling rate 
% timewin = round(timewinms / (1000/srate)); %timewin length in data point steps
% hz = linspace(0, srate, timewin);
% numfrex = sum(hz>min(frex) & hz<max(frex));
% L = size(trialDat,1); 
% ii = 1; %index variable for output
% winStart = L - round(timewin/2); %starting data index
% [curTim, timi] = min(times); %starting time and time index
% 
% %pad the data to prevent edge artifacts
% padDat = mirrorPad(trialDat); 
% 
% %initialize output variables
% trialTF = zeros([size(trialDat), numfrex]);
% outTim = zeros(size(times)); 
% frexOut = hz(hz>min(frex) & hz<max(frex));
% 
% 
% %prepare the tapers
% curDat = padDat(winStart:winStart+timewin-1, :);
% tapers    = dpss(size(curDat,1),3);
% tapers = tapers(:,1:end-1); %discard the last taper
% 
% 
% 
% 
% % Perform multitaper spectral estimation with matlab builtin
% [power, f] = pmtm(curDat, [],[],1000);
% trialTF2 = zeros([size(trialTF, [1,2]), ...
%                  sum(f>min(frex) & f<max(frex)) ]) ;
% 
% 
% check = true; 
% 
% while check
%     curDat = padDat(winStart:winStart+timewin-1, :);
% %     tapSpect = zeros(size(tapers,2),timewin, size(curDat,2)); 
% %     for tapi = 1:size(tapers,2)
% %         test = abs(fft(curDat.*tapers(:,tapi))).^2;
% %         tapSpect(tapi,:,:) = test(1:timewin,:); 
% %     end
% %     
% %     %store to the outputs
% %     spect = squeeze(mean(tapSpect(:,hz>min(frex) & hz<max(frex),:), 1));
% %     trialTF(ii,:,:) = spect';
% %     outTim(ii) = curTim; 
% 
%     %alternate calculation
%     [power, f] = pmtm(curDat, 3,[],1000);
%     trialTF2(ii, :, :) = power(f>min(frex)&f<max(frex),:)';
%     frexout2 = f(f>min(frex)&f<max(frex));
%     %update index variables
%     ii = ii+1; 
%     timStep = find(times(timi:end)>=curTim+timeStepms,1)-1;
%     if  isempty(timStep) %out beyond the end of time
%         check = false; 
%     else
%         timi = timi+timStep; 
%         curTim = times(timi);
%         winStart = winStart+timStep; 
%     end
%  
% end
% 
% %trim down outTim and trialTF
% outTim(ii:end) = []; 
% trialTF(ii:end, :, :) = []; 
% trialTF2(ii:end,:,:) = [];
% trialTF = trialTF2; %adjustment here to use pmtm for core calculation! 
% frexOut = frexout2;

%% end old version
end

%% testing scratch code below here

% pow = arrayfun(@(x) myChanZscore(trialTF(:,:,x), ...
%         [1, length(outTim)] ), ...
%         1:size(trialTF,3), 'UniformOutput',false ); %z-score
% highnumfrex = length(frexOut); 
% pow = cell2mat(pow); %organize
% pow = reshape(pow, size(pow,1), size(pow,2)/highnumfrex, []); %organize
% powMean = mean(pow,3); 
% 
% 
% 
% 
% pow2 = arrayfun(@(x) myChanZscore(trialTF2(:,:,x), ...
%         [1, length(outTim)] ), ...
%         1:size(trialTF2,3), 'UniformOutput',false ); %z-score
% highnumfrex = length(frexout2); 
% pow2 = cell2mat(pow2); %organize
% pow2 = reshape(pow2, size(pow2,1), size(pow2,2)/highnumfrex, []); %organize
% powMean2 = mean(pow2, 3); 
% 

% 
% RUN THE MIKE CODE BLOCK BELOW THIS BEFORE THIS SECTION!!!  
%
% pow3 = arrayfun(@(x) myChanZscore(tf(:,:,x), ...
%         [1, length(outTim)] ), ...
%         1:size(tf,3), 'UniformOutput',false ); %z-score
% highnumfrex = length(frexOut3); 
% pow3 = cell2mat(pow3); %organize
% pow3 = reshape(pow3, size(pow3,1), size(pow3,2)/highnumfrex, []); %organize
% powMean3 = mean(pow3, 3); 
% 
% lat3 = gausLat(squeeze(powMean3), ...
%                                outTim,...
%                                ones(size(powMean3,2))*1800);
% [~, order] = sort(lat3); 
% subplot 313
% imagesc(powMean3(:, order)')
% caxis([-10, 15])

%% Mike's code using present variable names to extract TF
% srate = 1000;
% 
% timewin    = 300; % in ms
% timewinidx = round(timewin/(1000/srate));
% tapers     = dpss(timewinidx,3); % this line will crash without matlab signal processing toolbox
% 
% 
% times2save   = min(times):timeStepms:max(times);
% 
% % convert time points to indices
% times2saveidx = dsearchn(times',times2save'); 
% 
% 
% % define frequencies for FFT
% hz = linspace(0,srate/2,timewinidx/2+1);
% 
% 
% % initialize output matrix
% tf = zeros(size(trialDat,2), floor(timewinidx/2)+1,length(times2save));
% 
% % loop through time bins
% for ti=1:length(times2saveidx)
%     
%     % initialize power vector (over tapers)
%     taperpow = zeros(floor(timewinidx/2)+1,size(trialDat,2));
%     
%     % loop through tapers
%     for tapi = 1:size(tapers,2)-1
%         
%         % get data from this time window and taper
%         tempEEG  = squeeze(padDat(L+times2saveidx(ti)-floor(timewinidx/2)+1:...
%                                   L+times2saveidx(ti)+ceil(timewinidx/2),:));
%         data     = bsxfun(@times,tempEEG,tapers(:,tapi));
%         
%         % compute FFT and extract power
%         powTemp      = fft(data)/timewinidx;
%         powTemp      = powTemp(1:floor(timewinidx/2)+1,:);
%         taperpow = taperpow + abs(powTemp).^2;
%     end
%     
%     % divide by N tapers for average
%     tf(:,:,ti) = taperpow'/tapi;
% end
% 
% 
% frexOut3 = hz(hz>min(frex) & hz<max(frex)); 
% tf = tf(:, hz>min(frex) & hz<max(frex), :);
% tf = permute(tf, [3,1,2]); 


%% end scratch testing code

% 
% 
% 
% winStarts = [1:timeStepms:size(zeroPadTD,1)-timewin];
% %time resolution  = 5ms (.05 seconds)
% 
% %time X trials X freq
% trialTF = zeros([length(winStarts), size(trialDat,2), numfrex]); 
% outTim = times([1:timeStepms:length(times)]);
% 
%  
% 
% 
% 
% for tt = 1:length(winStarts)
%     curDat = zeroPadTD(winStarts(tt):winStarts(tt)+timewin-1, :);
%     curDat = detrend(curDat); 
%     tapers    = dpss(size(curDat,1),5);
%     tapers = tapers(:,1:end-1); %discard the last taper
%     tapSpect = zeros(size(tapers,2),timewin, size(curDat,2)); 
%     for tapi = 1:size(tapers,2)
%         test = abs(fft(curDat.*tapers(:,tapi))).^2;
%         tapSpect(tapi,:,:) = test(1:timewin,:); 
%     end
%     
%     spect = squeeze(mean(tapSpect(:,hz>min(frex) & hz<max(frex),:), 1));
%     trialTF(tt,:,:) = spect'; 
% 
%     
% 
% end
%     
% 
% 





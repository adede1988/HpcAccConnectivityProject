function [ispc] = getChanISPC(trialDat, trialDat2, frex, numfrex, stds, ...
                                srate, di, trials, HFBlats, tim)

%trialDat:            timepoints X trials data points
%trialDat2:            timepoints X trials data points



%ispc stats are: ispc raleigh's Z, ppc raw, ispc raleigh's Z at HFB peak,
%ppc at HFB peak

ispc = zeros([length(di), numfrex, 4]); 

padDat = mirrorPad(trialDat(:,trials)); 
padDat2 = mirrorPad(trialDat2(:,trials)); 


time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate

n_wavelet            = length(time);
n_data               = prod(size(padDat)); 
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% hz = linspace(0, srate, size(n_convolution,1));
    
  
fftDat = fft(reshape(padDat,1,numel(padDat)),n_conv_pow2);
fftDat2 = fft(reshape(padDat2,1,numel(padDat2)),n_conv_pow2);
   


   
    for fi = 1:numfrex
       

        f  = frex(fi); % frequency of wavelet in Hz
        s  = stds(fi); 
        % and together they make a wavelet
        wavelet = sqrt(1/(s*sqrt(pi))) * ...
            exp(2*1i*pi*f.*time) .* exp(-time.^2./(2*(s^2))); 
        
        wavelet = fft(wavelet, n_conv_pow2); 

        %convolve and ifft, then extract key part
        eegconv = ifft(wavelet.*fftDat);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        temppower = reshape(eegconv,size(padDat));
        outTF = temppower(size(trialDat,1)+1:size(trialDat,1)*2, :); 
        %convolve and ifft, then extract key part
        eegconv = ifft(wavelet.*fftDat2);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        temppower = reshape(eegconv,size(padDat));
        outTF2 = temppower(size(trialDat,1)+1:size(trialDat,1)*2, :); 


    
        %downsample for speed and compactness 
        fDat = angle(outTF(di,:)); 
        fDat2 = angle(outTF2(di,:)); 


      
        %get ispc
        ispc(:,fi,1) = size(fDat,2)*arrayfun(@(x) ...
            abs(sum(exp(1i * (fDat(x,:) - fDat2(x,:))))./size(fDat,2)) ,...
            1:size(fDat,1) ).^2;
       

        %implementing pairwise phase consistency from Vink et al., 2010 
        difs = fDat - fDat2; 
        N = size(difs,2); 
        %equation 14
        ispc(:,fi,2) = mean(...
                    cell2mat(arrayfun(@(j) ...
                        cell2mat(arrayfun(@(k) ...
                            cos(difs(:,j)).*cos(difs(:,k)) + sin(difs(:,j)).*sin(difs(:,k)), ...
                        j+1:N, 'uniformoutput', false)), ...
                    1:N-1, 'uniformoutput', false)),...
                2);

       
        HFBidx = arrayfun(@(x) find(tim >= x, 1), HFBlats); 
        HFBidx(HFBidx<21) = 21; 
        HFBidx(HFBidx>length(tim)-20) = length(tim)- 20;
        %get the +- 20 points around the HFB peaks
        fDat = arrayfun(@(x, y) fDat(x-20:x+20, y), HFBidx', ...
            1:length(HFBidx), 'UniformOutput',false );
        fDat = cat(2, fDat{:});
        fDat2 = arrayfun(@(x, y) fDat2(x-20:x+20, y), HFBidx', ...
            1:length(HFBidx), 'UniformOutput',false );
        fDat2 = cat(2, fDat2{:});

        
         ispc(1:41,fi,3) = size(fDat,2)*arrayfun(@(x) ...
            abs(sum(exp(1i * (fDat(x,:) - fDat2(x,:))))./size(fDat,2)) ,...
            1:size(fDat,1) ).^2;
       

        %implementing pairwise phase consistency from Vink et al., 2010 
        difs = fDat - fDat2; 
        N = size(difs,2); 
        %equation 14
        ispc(1:41,fi,4) = mean(...
                    cell2mat(arrayfun(@(j) ...
                        cell2mat(arrayfun(@(k) ...
                            cos(difs(:,j)).*cos(difs(:,k)) + sin(difs(:,j)).*sin(difs(:,k)), ...
                        j+1:N, 'uniformoutput', false)), ...
                    1:N-1, 'uniformoutput', false)),...
                2);
       


    end










end





















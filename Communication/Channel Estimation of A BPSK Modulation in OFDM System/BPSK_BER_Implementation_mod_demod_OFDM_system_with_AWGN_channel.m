clearvars;
close all;
clc;

% Generating random bitstream
n = 1; m = 1280; % bitstream dimensions

subcarrierNumber = 10; % the number of subcarrier will be generated by multiplexing the bitstream.

SNR = 0:2:28; % SNR Range from 0 to 20 [in dB]

% Generation of random binary bitstream n by m = (i.e. 1 x 1280 )

bitstream = randi([0 1], n,m);

% BPSK Modulation

Serial_Modulated_Data=2*bitstream-1;

% Serial to parallel data conversion of the bitstream

Paralel_Modulated_Data(:,1) = Serial_Modulated_Data(1:m/subcarrierNumber);

for j = 2:subcarrierNumber
    startIndex = (m/subcarrierNumber)*(j-1) + 1;
    endIndex = (m/subcarrierNumber)*j;
    Paralel_Modulated_Data(:,j) = Serial_Modulated_Data(startIndex:endIndex);
end

% Pilot Insertion randomly

pilot1 = ones((m/subcarrierNumber),1);

% Pilot Signals: Randomly selected location

% pilot2 = randi([2,11]);
% pilot3 = randi([pilot1+1,12]);
% pilot4 = randi([pilot2+1,13]);
% pilot5 = randi([pilot3+1,14]);

% Pilot Signals: Determined location
pilot2 = 4;
pilot3 = 7;
pilot4 = 10;
pilot5 = 13;

pilotNumber = 5;

pilotInsertedData = [pilot1,Paralel_Modulated_Data];
twoPilotInsertedData = [pilotInsertedData(:,1:pilot2-1) pilot1 pilotInsertedData(:,pilot2:11)];
threePilotInsertedData = [twoPilotInsertedData(:,1:pilot3-1) pilot1 twoPilotInsertedData(:,pilot3:12)];
fourPilotInsertedData = [threePilotInsertedData(:,1:pilot4-1) pilot1 threePilotInsertedData(:,pilot4:13)];
fivePilotInsertedData = [fourPilotInsertedData(:,1:pilot5-1) pilot1 fourPilotInsertedData(:,pilot5:14)];

% Inverse FFT
dataInTimeDomain = ifft(fivePilotInsertedData);

% Cyclic Prefix Insertion
cyclic_prefix_rate = 1/8;
subcarrierSize = length(bitstream(:))/subcarrierNumber;
cyclic_prefix_Number = subcarrierSize * cyclic_prefix_rate;

out_with_cyclicPrefix = [dataInTimeDomain(((m/subcarrierNumber)-(m/subcarrierNumber * cyclic_prefix_rate)+1):(m/subcarrierNumber),:);dataInTimeDomain];

%Generating a channel its impulse response
h = [0.10 0.10 0.90 0.10 0.10]; % 5-tap filter

% Convolution with channel h
h_conv_x = ones(length(out_with_cyclicPrefix(:,1))+length(h)-1, subcarrierNumber+pilotNumber);

for i = 1:(subcarrierNumber+pilotNumber)
    h_conv_x(:,i) = conv(out_with_cyclicPrefix(:,i),h);
end

% Generated noisy output
iterNum = 1000;

% Avoiding variable dimension settings for the matrices below in a loop.
without_cyclic_prefix = ones(subcarrierSize,subcarrierNumber+length(h));
estimated_ErrorNumber = ones(iterNum,length(SNR));
perfect_ErrorNumber = ones(iterNum,length(SNR));

% Monte Carlo Simulation
for iter=1:iterNum
    for i = 1:length(SNR)
        y_noise = awgn(h_conv_x,SNR(i),'measured');
    
        %Removed cyclic prefix and additional transition bits
        for k=1:(m/subcarrierNumber)
            for j=1:(subcarrierNumber+pilotNumber)
                without_cyclic_prefix(k,j)= y_noise(k+cyclic_prefix_Number,j);
            end
        end

        %FFT of Removed Cyclic Prefix
        fft_wcycpre = fft(without_cyclic_prefix);
    
        %Creating new pilot vector
        pilot_new = zeros(m/subcarrierNumber,1);
    
        for t=1: (m/subcarrierNumber)
            pilot_new(t,:) = ( fft_wcycpre(t,1) + fft_wcycpre(t,pilot2) + fft_wcycpre(t,pilot3) + fft_wcycpre(t,pilot4) + fft_wcycpre(t,pilot5) ) / pilotNumber;
        end
        
        %Removing pilots from original data
         fft_wcycpre(:,1) = [];
         fft_wcycpre(:,(pilot2-1)) = [];
         fft_wcycpre(:, (pilot3-2)) = [];
         fft_wcycpre(:, (pilot4-3)) = [];
         fft_wcycpre(:, (pilot5-4)) = [];
   
        % %FFT of the channel to estimate perfect data
        h_fft = fft(h,m/subcarrierNumber);
        h_fft_t = transpose(h_fft);
        estimated_perfect = zeros(m/subcarrierNumber,subcarrierNumber);
        estimated_data = zeros(m/subcarrierNumber,subcarrierNumber);

        %Estimating data by dividing elements of column to the elements of the pilot column elemenetwisely
        for x = 1:(subcarrierNumber)
            for y = 1:(m/subcarrierNumber)
                estimated_data(y,x) = fft_wcycpre(y,x) ./ pilot_new(y,1);
                estimated_perfect(y,x) = fft_wcycpre(y,x) ./ h_fft_t(y,1);
            end
        end

        estimated_perfect_serial = zeros(1,m);
        estimated_serial = zeros(1,m);

        %Parallel to Serial Data Conversion
        estimated_serial(1,1:m/subcarrierNumber) = transpose(estimated_data(:,1));

        estimated_perfect_serial(1,1:m/subcarrierNumber) = transpose(estimated_perfect(:,1));

        for j = 2: (subcarrierNumber)
            estimated_serial(1,(m/subcarrierNumber*(j-1)+1):m/subcarrierNumber*j) = transpose(estimated_data(:,j));
            estimated_perfect_serial(1,(m/subcarrierNumber*(j-1)+1):m/subcarrierNumber*j) = transpose(estimated_perfect(:,j));
        end

        % BPSK Demodulation
        demodulated_data = zeros(1,m);
        demodulated_perfect_data = zeros(1,m);

        for u = 1:m
            if real(estimated_serial(u)) < 0
                demodulated_data(u) = 0;
            else
                demodulated_data(u) = 1;
            end
        end

        for u = 1:m
            if real(estimated_perfect_serial(u)) < 0
                demodulated_perfect_data(u) = 0;
            else
                demodulated_perfect_data(u) = 1;
            end
        end

        % Calculate BER vs SNR
        estimated_ErrorNumber(iter,i)=biterr(bitstream,demodulated_data);
        perfect_ErrorNumber(iter,i)=biterr(bitstream,demodulated_perfect_data);
    end
end

BER_estimated=mean(estimated_ErrorNumber)/length(bitstream);
BER_perfect=mean(perfect_ErrorNumber)/length(bitstream);

semilogy(SNR,BER_estimated,'b');
hold on;
semilogy(SNR,BER_perfect,'r');
title('5-Pilot & 5-Tap with Cyclic Prefix Insertion');
xlabel('SNR(dB)'); ylabel('BER'); legend('Estimated Data','Perfect Data');
xlim([0, max(SNR)]);
ylim([1.0e-6, 1]);
grid on

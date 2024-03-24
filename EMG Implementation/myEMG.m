clearvars;
close all;
clc;


%% EMG Data Okuma
lines = textread('emgData1.txt','%s','delimiter','\r');

emgData = [];
index = 1; % emgData dizisine ekleme yapmak için indeks tutar

for ith_line = 1:length(lines)
    line = lines{ith_line};
    data = str2num(line);
    
    % Eğer veri varsa ve en az 2 elemana sahipse, emgData'ya ekle
    if ~isempty(data) && length(data) > 1
        emgData(index,:) = data(1:2);
        index = index + 1; % Bir sonraki indeksi güncelle
    end
end

clear line;
clear lines;
clear ith_line;
clear data;

%% Verileri ve Zaman bilgilerini ayıklama
data1 = emgData(:,2);
time1 = emgData(:,1);

%data1 = data1(1:8000);
%time1 = time1(1:8000);
% data1 = data1(21700:29200);
% time1 = time1(21700:29200);
% time1 = linspace(0,30,length(time1));
% time1 = transpose(time1);
time1 = linspace(0,30.161,60323);
time1 = transpose(time1);

data1 = data1(1:60001);
time1 = time1(1:60001);

% Örnekleme frekansı fs = 2kHz, 30sn lik kayıt sinyali için 60000 örnek 
% 0 periyotu içinde 1 örnek;



%% Step 1: Raw veri Plot
figure;
subplot(1,1,1);
plot(time1(:), data1(:)); 
xlabel('Zaman (s)','fontsize', 14); ylabel('Sinyal Genliği (mV)', 'fontsize', 14); 
title('İşlenmemiş EMG Verisi- Sağlıklı', 'fontsize', 14);set(gca,'FontSize',14);
ylim([0 1]);

%% Step 2: Perform Fast Fourier Transform and Shift Zero-component to center

fq = 2000; % örnekleme frekansı

sEMG1 = data1;
n = length(sEMG1);

% Hızlı Fourier Dönüşümü (FFT) gerçekleştir
fft_sEMG1 = fft(sEMG1);

% FFT sonucunu merkeze taşı
fft_sEMG1 = fftshift(fft_sEMG1);

% Sıfır-merkezli frekans aralığını oluştur
fq_axis1 = (-n/2:n/2-1)*(fq/n);

% FFT sonucunun mutlak değerini al
abs_fft_sEMG1 = abs(fft_sEMG1);

% Grafiği çiz
figure;
plot(fq_axis1, abs_fft_sEMG1);
axis([fq_axis1(1) fq_axis1(end) 0 40]); 
xlabel('Frekans (Hz)','fontsize', 14); 
ylabel('|X(jw)| büyüklüğü','fontsize', 14); 
title('Frekans Alanı - Sağlıklı','fontsize', 14);
set(gca,'FontSize',14);


%% Step 3: Generate Bandpass filter 
% can change to your desired cutoff frequency
highpass = 5;   
lowpass = 10;    

s = 12.5;

% get the index of -10 to -5 and 5 to 10Hz. 
cutoff1 = ceil((s-highpass)/(fq/length(sEMG1))); cutoff2 = ceil((s-lowpass)/(fq/length(sEMG1)));
cutoff3 = ceil((highpass+s)/(fq/length(sEMG1))); cutoff4 = ceil((lowpass+s)/(fq/length(sEMG1)));
% cutoff1 = -lowpass; cutoff2 = -highpass;
% cutoff3 = lowpass; cutoff4 = highpass;

H = zeros(length(sEMG1),1);
H(cutoff2:cutoff1) = 1; % take only the -10 to -5Hz
H(cutoff3:cutoff4) = 1; % take only the 5 to 10Hz
% H(6000:18000) = 1; % take only the -10 to -5Hz
% H(42000:54001) = 1; % take only the 5 to 10Hz

fq_axis2 = linspace(-s, s, 751);
h = ifftshift(H);
h = ifft(H);
plot(time1,real(h));
figure; plot(fq_axis2, H(1:751)); set(gca,'YLim',[0 2]); xlabel('Freqeuncy/Hz','fontsize', 14); ylabel('Amplitude','fontsize', 14);
title('Bandpass filter','fontsize', 14);set(gca,'FontSize',14);

fq_axis_new = linspace(min(fq_axis2), max(fq_axis2), 60001); % Yeni frekans ekseni
H = interp1(1:751, H(1:751), linspace(1, 751, 60001)); % H'yi upsample et
plot(fq_axis_new, H);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Upsampled Frequency Response');



%% Step 4: Perform Bandpass filter, inverse shifting and inverse Fourier transform

cutoff1 = ceil((s-highpass)/(fq/length(sEMG1))); cutoff2 = ceil((s-lowpass)/(fq/length(sEMG1)));
cutoff3 = ceil((highpass+s)/(fq/length(sEMG1))); cutoff4 = ceil((lowpass+s)/(fq/length(sEMG1)));

H = zeros(length(sEMG1),1); 
% H(cutoff2:cutoff1) = 1; % take only the -10 to -5Hz
% H(cutoff3:cutoff4) = 1; % take only the 5 to 10Hz
H(6000:18000) = 1; % take only the -10 to -5Hz
H(42000:54001) = 1; % take only the 5 to 10Hz

yt1 = ifftshift(fft_sEMG1.*H); 
yt1 = ifft(yt1);

figure;
subplot(1,1,1);
plot(time1, real(yt1),'r');xlabel('Zaman (s)','fontsize', 14); ylabel('Genlik','fontsize', 14); 
title('BPF çıkış Sinyali - Sağlıklı','fontsize', 14);set(gca,'FontSize',14);

%% Ağırlık Düzeltme
weightAdj = 22;%7.8420853796751038911975821684926;
yt1 = yt1 * weightAdj;

%% Doğrultma

yt2 = abs(yt1);

figure;
subplot(1,1,1);
plot(time1, real(yt2),'r');xlabel('Zaman (s)','fontsize', 14); ylabel('Genlik','fontsize', 14); 
title('BPF çıkış Sinyali - Sağlıklı','fontsize', 14);set(gca,'FontSize',14);


%% Step 5: Feature extraction - Moving Average Filter. (Smoothing)

filterWindowSize = 1000;

% Hareketli Ortalama için filtre oluştur
filter = ones(filterWindowSize,1) / filterWindowSize;

% Hareketli Ortalama (Moving Average) uygula
MovAvgAppliedData = movmean(abs(real(yt2)), filterWindowSize);




%% Step 6: Plot EMG

figure;
subplot(1,1,1);
plot(time1,MovAvgAppliedData,'linewidth',2); hold on;
plot(time1,real(yt2),'r'); hold off;
legend('Hareketli Ortalama','EMG sinyali', 'location', 'Northeast');
xlabel('Zaman (s)'); ylabel('Sinyal Genliği (mV)'); title('Sağlıklı EMG Verisi'); grid on;
ylim([min(real(yt2)) max(real(yt2))]);
set(gca,'FontSize',14);

%% Step 7: RMS Value of the EMG signal
RMSstep = 32;
n = length(yt2);

% RMS değerlerini hesapla
emgRMS = zeros(1, ceil(n / RMSstep));

for i = 1:length(emgRMS)
    start_index = (i - 1) * RMSstep + 1;
    end_index = min(i * RMSstep, n);
    emgRMS(i) = rms(real(yt2(start_index:end_index)));
end

% Tekrarlama yaparak EMG RMS değerlerini genişlet
emgRms = repelem(emgRMS, RMSstep);

% Son RMS değerini ekleyin
emgRms(end+1:n) = emgRMS(end);

% time1 ve emgRms vektörlerinin uzunluklarını kontrol edin ve ayarlayın
if length(time1) > length(emgRms)
    time1 = time1(1:length(emgRms));
elseif length(emgRms) > length(time1)
    emgRms = emgRms(1:length(time1));
end

% Grafiği çiz
figure;
plot(time1, MovAvgAppliedData, 'linewidth', 2);
hold on;
plot(time1, emgRms, 'r', 'linewidth', 2);
hold off;
legend('Hareketli Ortalama', 'RMS EMG sinyali', 'location', 'Northeast');
xlabel('Zaman (s)');
ylabel('Sinyal Genliği (mV)');
title('Sağlıklı EMG Verisi');
grid on;
ylim([min(yt2/2) max(yt2/2)]);
set(gca, 'FontSize', 14);

%% Contraction Detection or Classification
ClassCoef = 0.05;

% EMG sinyalini sınıflandırma
emgContraction = zeros(size(emgRms)); % Sınıflandırılmış EMG sinyali

% Sınırlayıcı değeri karşılaştır
emgContraction(emgRms > ClassCoef) = 0.2;


figure;
subplot(1,1,1);
plot(time1,real(yt1),'linewidth',2); hold on;
plot(time1,MovAvgAppliedData,'linewidth',2); hold on;
plot(time1,emgRms,'r','linewidth',2); hold on;
plot(time1,emgContraction,'b','linewidth',2); hold off;
legend('EMG Sinyali','Hareketli Ortalama','RMS EMG sinyali','EMG Kasılma Sinyali', 'location', 'Northeast');
xlabel('Zaman (s)'); ylabel('Sinyal Genliği (mV)'); title('Sağlıklı EMG Verisi'); grid on;
ylim([min(real(yt1)) max(real(yt1))]);
set(gca,'FontSize',14);

%%%% Contraction Signal 2 for 1,0 -valued classification

ClassCoef = 0.05;

% EMG sinyalini sınıflandırma
emgContraction2 = double(emgRms > ClassCoef);


%% Modulation of Contraction and MOV signals
for i= 1:length(emgRms)
    emgSignal(i) = MovAvgAppliedData(i) * emgContraction2(i);
end

%% Gain Adjustment
MovAvgDataGain = 1;
emgRmsGain = 1;
emgContractionGain = 1;
emgGain = 2;
MovAvgAppliedData = MovAvgDataGain * MovAvgAppliedData;
emgRms = emgRmsGain * emgRms;
emgContraction = emgContractionGain * emgContraction;
emgSignal = emgGain * emgSignal;

%% Plot EMG
figure;
subplot(1,1,1);
plot(time1,real(yt1),'linewidth',2); hold on;
plot(time1,MovAvgAppliedData,'linewidth',2); hold on;
plot(time1,emgRms,'r','linewidth',2); hold on;
plot(time1,emgContraction,'k','linewidth',2); hold on;
plot(time1,emgSignal,'linewidth',2); hold off;
legend('EMG Sinyali','Hareketli Ortalama','RMS EMG sinyali','EMG Kasılma Sinyali','EMG çıkış Sinyali', 'location', 'Northeast');
xlabel('Zaman (s)'); ylabel('Sinyal Genliği (mV)'); title('Sağlıklı EMG Verisi'); grid on;
%ylim([min(real(yt1)) max(real(yt1))]);
set(gca,'FontSize',14);

figure;
subplot(1,1,1);
plot(time1,real(yt1),'linewidth',2); hold on;
plot(time1,emgContraction,'k','linewidth',4); hold on;
plot(time1,emgSignal,'linewidth',2); hold off;
legend('EMG Sinyali','EMG Kasılma / Gevşeme Sinyali','EMG çıkış Sinyali', 'location', 'Northeast');
xlabel('Zaman (s)'); ylabel('Sinyal Genliği (mV)'); title('Sağlıklı EMG Verisi'); grid on;
%ylim([min(real(yt1)) max(real(yt1))]);
set(gca,'FontSize',14);

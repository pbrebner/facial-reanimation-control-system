%Set up data for system identification toolbox
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 5000;             % Length of signal
t = (0:L-1)*T;        % Time vector

load EMG_Whisk_LongImplant.mat
trial='1';
animal ='1';
eval(['Stim=H',animal,'Stimulus' trial ';']);
eval(['EMG=H',animal,'EMG' trial ';' ]);
eval(['Whisk=H',animal,'Whisk' trial ';' ]);

hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',10,'PassbandRipple',0.2, ...
         'SampleRate',Fs);   %High Pass Filter
%fvtool(hpFilt)

bpFilt1 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter           
%fvtool(bpFilt1)

bpFilt2 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',126,'HalfPowerFrequency2',130, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter           
%fvtool(bpFilt2)

bpFilt3 = designfilt('bandpassiir','FilterOrder',8, ...
               'HalfPowerFrequency1',55,'HalfPowerFrequency2',65, ...
               'DesignMethod','butter','SampleRate',Fs);   % Band Pass filter           
%fvtool(bpFilt3)

Stim=decimate(Stim,10);
EMG=decimate(EMG,10);

EMG = EMG - mean(EMG);
EMG = detrend(EMG);
Whisk = Whisk-mean(Whisk);

%EMG = filtfilt(bpFilt1,EMG);
%EMG = filtfilt(bpFilt2,EMG);
%EMG = filtfilt(hpFilt,EMG);

%%
%Zoom in on Stim Burst
figure(1)
plot(tWhisk(:,2450:2600),Stim(2450:2600,:))
title('Stimulus Burst','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Voltage (V)','Fontsize',15)

%%

%Plot
figure(2)
subplot(3,1,1)
plot(tWhisk,Stim)
title('Stimulus','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Voltage (V)','Fontsize',15)

subplot(3,1,2)
plot(tWhisk,EMG)
title('EMG','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Voltage (V)','Fontsize',15)

subplot(3,1,3)
plot(tWhisk,Whisk)
title('Whisker Displacement','Fontsize',18)
xlabel('Time (s)','Fontsize',15)
ylabel('Angle (degrees)','Fontsize',15)


%%
%PDF of EMG
figure(3)
subplot(2,1,1)
histogram(EMG,'Normalization','pdf')

xlabel('Volts (V)','Fontsize',20)
ylabel('Density','Fontsize',20)
title('EMG Distribution','Fontsize',24)
%legend('Estimated','Theoretical')

%Power Spectrum of EMG
Nfft = 1000;
[Pxx,f] = pwelch(EMG,gausswin(Nfft),Nfft/2,Nfft,Fs);

% Plot frequency spectrum
subplot(2,1,2)
plot(f,Pxx);
title('Power Spectrum of EMG Signal','Fontsize',24);
ylabel('PSD','Fontsize',20); 
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

%PDF of Whisker Output
figure(4)
subplot(2,1,1)
histogram(Whisk,'Normalization','pdf')

xlabel('Volts (V)','Fontsize',20)
ylabel('Density','Fontsize',20)
title('Whisk Displacement Distribution','Fontsize',24)
%legend('Estimated','Theoretical')

%Power Spectrum of Whisker Output
Nfft = 1000;
[Pxx,f] = pwelch(Whisk,gausswin(Nfft),Nfft/2,Nfft,Fs);

% Plot frequency spectrum
subplot(2,1,2)
plot(f,Pxx);
title('Power Spectrum of Whisk Displacement','Fontsize',24);
ylabel('PSD','Fontsize',20); 
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

%PDF of Stimulus
figure(5)
subplot(2,1,1)
histogram(Stim,'Normalization','pdf')

xlabel('Volts (V)','Fontsize',20)
ylabel('Density','Fontsize',20)
title('Stimulus Distribution','Fontsize',24)
%legend('Estimated','Theoretical')

%Power Spectrum of Whisker Output
Nfft = 1000;
[Pxx,f] = pwelch(Stim,gausswin(Nfft),Nfft/2,Nfft,Fs);

% Plot frequency spectrum
subplot(2,1,2)
plot(f,Pxx);
title('Power Spectrum of Stimulus','Fontsize',24);
ylabel('PSD','Fontsize',20); 
xlabel('Frequency (Hz)','Fontsize',20);
grid on;


%%
%FFT
U = fft(EMG);
P2 = abs(U/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(6)
plot(f,P1) 
title('Mean Spectra of Input','Fontsize',32)
xlabel('Frequency (Hz)','Fontsize',30)
ylabel('|P1(f)|','Fontsize',30)


%%
%Relationship between the stimulus and EMG
Stim_loss = Stim; 
Noise = EMG/Stim_loss;
figure(7)
plot(Noise)
title('Noise')

Filtered_EMG = filtfilt(bpFilt3,EMG);
figure(8)
plot(Filtered_EMG)
title('Filtered EMG')

%Autocorrelation of Stimulus
figure(9)
maxlag = L*0.05;
[Rxx,lags] = xcorr(Filtered_EMG,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,1),plot(lags,Rxx)
title('Auto-Correlation Function of Stimulus','Fontsize',30)
xlabel('Time (s)','Fontsize',24)
ylabel('Correlation','Fontsize',22)
grid on

%Auto Correlation of EMG
[Ryy,lags] = xcorr(Stim_loss,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,2),plot(lags,Ryy)
title('Auto-Correlation Function of EMG','Fontsize',30)
xlabel('Time (s)','Fontsize',24)
ylabel('Correlation','Fontsize',22)
grid on

%Cross Correlation of Stimulus/EMG
[Rxy,lags] = xcorr(Stim_loss,Filtered_EMG,maxlag,'coeff');
lags = lags/Fs;
subplot(3,1,3),plot(lags,Rxy)
title('Cross-Correlation Function (Stim/EMG)','Fontsize',30)
xlabel('Time (s)','Fontsize',24)
ylabel('Correlation','Fontsize',22)
grid on


%%
%Take Absolute of EMG
EMG = abs(EMG);

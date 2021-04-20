%Set up data for system identification toolbox
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 50000;             % Length of signal
t = (0:L-1)*T;        % Time vector

load EMG_Whisk_NEW2.mat
trial='5';
eval(['Stim=Stimulus' trial ';']);
eval(['EMG=EMG' trial ';' ]);
eval(['Whisk=Whisk' trial ';' ]);

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

Stim=decimate(Stim,20);
EMG=decimate(EMG,20);
%EMG = filtfilt(bpFilt1,EMG);
%EMG = filtfilt(bpFilt2,EMG);
%EMG = filtfilt(hpFilt,EMG);

EMG = EMG - mean(EMG);
EMG = detrend(EMG);
Whisk = Whisk-mean(Whisk);

%%
%PDF
figure(1)
subplot(2,1,1)
histogram(EMG,'Normalization','pdf')

xlabel('Volts (mv)')
ylabel('Density')
title('EMG Distribution')
%legend('Estimated','Theoretical')

%Power Spectrum
Nfft = 1000;
[Pxx,f] = pwelch(EMG,gausswin(Nfft),Nfft/2,Nfft,Fs);

% Plot frequency spectrum
subplot(2,1,2)
plot(f,Pxx);
title('Power Spectrum of Signal','Fontsize',24);
ylabel('PSD','Fontsize',20); 
xlabel('Frequency (Hz)','Fontsize',20);
grid on;

%FFT
U = fft(EMG);
P2 = abs(U/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure(2)
plot(f,P1) 
title('Mean Spectra of Input','Fontsize',32)
xlabel('Frequency (Hz)','Fontsize',30)
ylabel('|P1(f)|','Fontsize',30)

%Plot
figure(3)
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

%Take Absolute of EMG
EMG = abs(EMG);
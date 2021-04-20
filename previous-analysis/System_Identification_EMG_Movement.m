%System Identification of EMG-Movement Model
%Based on Real Data
clear all
clc

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 5000;             % Length of signal
t = (0:L-1)*T;        % Time vector

%Load experimental data
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

bsFilt1 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',55,'HalfPowerFrequency2',65, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter           
%fvtool(bsFilt1)

EMG = filtfilt(bsFilt1,EMG);
EMG = filtfilt(hpFilt,EMG);

z=cat(2, Stim, EMG);
Z=nldat(z,'domainIncr',1/10000,'chanNames', {'Stimulus' 'EMG'});
Z=decimate(Z,10);
W=nldat(Whisk,'domainIncr',.001, 'chanNames', {'Whisk'});
Z=cat(2,Z,W);
Z1=Z(:,[2 3]);
Z1=Z1-mean(Z1(1:1000,:)); 
Z1(:,1)=abs(Z1(:,1));

figure(25);
plot(Z1);

%%
%System Identification between EMG and Whisker Displacement
Zcur=detrend(Z1);

NHK=nlbl;
set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);
% Set umber of lags in IRF
i=NHK{1,2};
set(i,'nLags',150);
NHK{1,2}=i;

NHK=nlident(NHK,Zcur);
figure(26); 
plot(NHK); 
figure(27);
[R, V, yp] = nlid_resid(NHK,Zcur);

figure(28)
subplot(3,1,1)
plot(R)
title('Residuals of EMG/Whisk Hammerstein Model','Fontsize',20)
xlabel('Time (s)','Fontsize',18)
ylabel('Whisk (rad)','Fontsize',18)

subplot(3,1,2)
p = pdf(R);
plot(p)
xlabel('Volts (V)','Fontsize',18)
ylabel('Density','Fontsize',18)
title('Residual Distribution','Fontsize',20)

S = spect(R);
subplot(3,1,3)
plot(S);
title('Power Spectrum of Residuals','Fontsize',20);
ylabel('PSD','Fontsize',18); 
xlabel('Frequency (Hz)','Fontsize',18);
grid on

disp(NHK.idMethod)

%%
%Frequency response of IRF

F = fresp(NHK{1,2});
figure(29)
plot(F)

%%
%Below this is not currently used

%%
%m = nlbl(Zcur,'idMethod','subspace','nDelayInput',50,'orderLE',2);
%m = nlbl(Zcur,'idMethod','subspace');
%m = nlbl(Zcur);

%figure(2)
%plot(m);


%%
%Validate the model (%VAF)
wPre=nlsim(m,Zcur(:,1));
v=vaf(Zcur(:,2),wPre);
    
plot (Zcur(:,2));
h=line(wPre);set(h,'color','r');
title(['VAF= ' num2str(v)],'Fontsize',20);
xlabel('Time (s)','Fontsize',18)
ylabel('Whisk Displacement','Fontsize',18)
legend('Expected','Predicted','Fontsize',15);


%%
%Analyze the Residuals

figure(4)
R = Zcur(:,2) - wPre;
plot(R)
title('Residuals of EMG/Whisk Hammerstein Model')
    
figure(5)
subplot(2,1,1)
p = pdf(R);
plot(p)
xlabel('Volts (V)','Fontsize',20)
ylabel('Density','Fontsize',20)
title('Residual Distribution','Fontsize',24)

% Plot frequency spectrum
S = spect(R);
subplot(2,1,2)
plot(S);
title('Power Spectrum of Residuals','Fontsize',24);
ylabel('PSD','Fontsize',20); 
xlabel('Frequency (Hz)','Fontsize',20);
grid on
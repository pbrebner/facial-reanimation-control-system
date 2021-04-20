%System Identification of Stimulus-Movement Model
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
Z2=Z(:,[1 3]);
Z2=Z2-mean(Z2(1:1000,:)); 
%Z1(:,1)=abs(Z1(:,1));

figure(50);
plot(Z2);

%%
%System Identification between EMG and Whisker Displacement
%Zcur = Z2;

Stim_decimated = decimate(Stim,10);
Whisk_mean=Whisk-mean(Whisk(1:1000,:)); 
Zcur = iddata(Whisk_mean,Stim_decimated,Fs);

Orders = [3 4 1];

InputNL = sigmoidnet;
InputNL.NumberOfUnits = 20;

OutputNL = sigmoidnet;
OutputNL.NumberOfUnits = 20;

opt = nlhwOptions;
opt.Display = 'on';
%opt.SearchOptions.MaxIterations = 80;

%NLN = nlhw(Zcur,Orders,InputNL,OutputNL,opt);
NLN = nlhw(Zcur,Orders,opt);

figure(51)
compare(Zcur,NLN)

%%
%System Identification between EMG and Whisker Displacement
Zcur = Z2;

LNL=lnlbl;
set(LNL,'idMethod','hk','nLags1',150,'nLags2',150, 'polyOrderMax',4);

LNL=nlident(LNL,Zcur);
figure(52); 
plot(LNL); 
figure(53);
[R, V, yp] = nlid_resid(LNL,Zcur);

figure(54)
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
plot(S(1:250,:));
title('Power Spectrum of Residuals','Fontsize',20);
ylabel('PSD','Fontsize',18); 
xlabel('Frequency (Hz)','Fontsize',18);
grid on

disp(LNL.idMethod)


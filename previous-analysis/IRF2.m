%% EMG Anal
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
               'HalfPowerFrequency1',55,'HalfPowerFrequency2',65, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter           
%fvtool(bpFilt1)

bpFilt2 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',126,'HalfPowerFrequency2',130, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter           
%fvtool(bpFilt2)

EMG = filtfilt(bpFilt1,EMG);
%EMG = filtfilt(bpFilt2,EMG);
EMG = filtfilt(hpFilt,EMG);

z=cat(2, Stim, EMG);
Z=nldat(z,'domainIncr',1/10000,'chanNames', {'Stimulus' 'EMG'});
Z=decimate(Z,10);
W=nldat(Whisk,'domainIncr',.001, 'chanNames', {'Whisk'});
Z=cat(2,Z,W);
Z1=Z(:,[2 3]);
Z1=Z1-mean(Z1(1:1000,:)); 
Z1(:,1)=abs(Z1(:,1));
figure(1);
plot(Z1);

Z2 = cat(2,Z(:,1),Z1);
Z3 = Z2(:,[1 3]);
Z2 = Z2(:,[1 2]);
figure(2)
plot(Z2)

% Gain changes with amplitude so look at 10 s intervals. 

%%
%Relationship between EMG and Whisker Displacement
iNum=0;
gVal=[];
V = [];

for iStart=1
    iNum=iNum+1;
    figure(3); clf
    iLen=5000;
    iEnd=iStart+iLen-1;
    Zcur=Z1(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags',200);
    subplot(2,2,[1 2]);
    plot(i);
    title('IRF (EMG/Whisk)','Fontsize',20)
    xlabel('Time (s)','Fontsize',18)
    ylabel('IRF','Fontsize',18)

    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre);
    

    plot (Zcur(:,2));
    h=line(wPre);set(h,'color','r');
    title(['VAF= ' num2str(v)],'Fontsize',20);
    xlabel('Time (s)','Fontsize',18)
    ylabel('Whisk Displacement','Fontsize',18)
    legend('Expected','Predicted','Fontsize',15);
    
    %Residuals
    figure(4)
    R = Zcur(:,2) - wPre;
    plot(R)
    title('Residuals of EMG/Whisk IRF')
    
    figure(5)
    subplot(2,1,1)
    p = pdf(R);
    plot(p)

    xlabel('Volts (V)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('Residual Distribution','Fontsize',24)
    %legend('Estimated','Theoretical')

    %Power Spectrum of EMG
    Nfft = 1000;
    %[Pxx,f] = pwelch(R,gausswin(Nfft),Nfft/2,Nfft,Fs);
    S = spect(R);

    % Plot frequency spectrum
    subplot(2,1,2)
    plot(S);
    title('Power Spectrum of Residual','Fontsize',24);
    ylabel('PSD','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on

end

%%
%Relationship between Stimulus and Whisker Displacement
iNum=0;
gVal=[];
V = [];

for iStart=1
    iNum=iNum+1;
    figure(6); clf
    iLen=5000;
    iEnd=iStart+iLen-1;
    Zcur=Z3(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags',200);
    subplot(2,2,[1 2]);
    plot(i);
    title('IRF (Stimulus/Whisk)','Fontsize',20)
    xlabel('Time (s)','Fontsize',18)
    ylabel('IRF','Fontsize',18)

    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre);

    plot (Zcur(:,2));
    h=line(wPre);set(h,'color','r');
    title(['VAF= ' num2str(v)],'Fontsize',20);
    xlabel('Time (s)','Fontsize',18)
    ylabel('Whisk Displacement','Fontsize',18)
    legend('Expected','Predicted','Fontsize',15);
    
    %Residuals
    figure(7)
    R = Zcur(:,2) - wPre;
    plot(R)
    title('Residuals of Stimulus/Whisk IRF')

    figure(8)
    subplot(2,1,1)
    p = pdf(R);
    plot(p)

    xlabel('Volts (V)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('Residual Distribution','Fontsize',24)
    %legend('Estimated','Theoretical')

    %Power Spectrum of EMG
    Nfft = 1000;
    %[Pxx,f] = pwelch(R,gausswin(Nfft),Nfft/2,Nfft,Fs);
    S = spect(R);

    % Plot frequency spectrum
    subplot(2,1,2)
    plot(S);
    title('Power Spectrum of Residual','Fontsize',24);
    ylabel('PSD','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on    
    
end

%%
%Relationship between Stimulus and EMG
iNum=0;
gVal=[];
V = [];

for iStart=1
    iNum=iNum+1;
    figure(9); clf
    iLen=5000;
    iEnd=iStart+iLen-1;
    Zcur=Z2(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags',200);
    subplot(2,2,[1 2]);
    plot(i);
    title('IRF (Stimulus/EMG)','Fontsize',20)
    xlabel('Time (s)','Fontsize',18)
    ylabel('IRF','Fontsize',18)

    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre);

    plot (Zcur(:,2));
    h=line(wPre);set(h,'color','r');
    title(['VAF= ' num2str(v)],'Fontsize',20);
    xlabel('Time (s)','Fontsize',18)
    ylabel('Whisk Displacement','Fontsize',18)
    legend('Expected','Predicted','Fontsize',15);
    
    %Residuals
    figure(10)
    R = Zcur(:,2) - wPre;
    plot(R)
    title('Residuals of Stimulus/EMG IRF')

    figure(11)
    subplot(2,1,1)
    p = pdf(R);
    plot(p)

    xlabel('Volts (V)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('Residual Distribution','Fontsize',24)
    %legend('Estimated','Theoretical')

    %Power Spectrum of EMG
    Nfft = 1000;
    %[Pxx,f] = pwelch(R,gausswin(Nfft),Nfft/2,Nfft,Fs);
    S = spect(R);

    % Plot frequency spectrum
    subplot(2,1,2)
    plot(S);
    title('Power Spectrum of Residual','Fontsize',24);
    ylabel('PSD','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on
    
end
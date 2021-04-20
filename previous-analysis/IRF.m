%% EMG Anal
clear all
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 5000;             % Length of signal
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
Z=nldat(z,'domainIncr',1/20000,'chanNames', {'Stimulus' 'EMG'});
Z=decimate(Z,20);
W=nldat(Whisk,'domainIncr',.001, 'chanNames', {'Whisk'});
Z=cat(2,Z,W);
Z1=Z(:,[2 3]);
Z1=Z1-mean(Z1(1:1000,:));
plot(Z1)
Z1(:,1)=abs(Z1(:,1));
figure(1);
plot(Z1);

%Zcur=Z1(iStart:iEnd,:);
figNum=2;
% Gain changes with amplitude so look at 10 s intervals. 
iNum=0;
gVal=[];
V=[];
for iStart=1:10000:50000
    iNum=iNum+1;
    figure(figNum); clf
    iLen=10000;
    iEnd=iStart+iLen-1;
    Zcur=Z1(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags',200);
    subplot(2,2,[1 2]);
    plot(i);
    title('IRF','Fontsize',20)
    xlabel('Time (s)','Fontsize',18)
    ylabel('IRF','Fontsize',18)

    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre);
    V(iNum)=v;
    plot (Zcur(:,2));
    h=line(wPre);set(h,'color','r');
    title(['VAF= ' num2str(v)],'Fontsize',20);
    xlabel('Time (s)','Fontsize',18)
    ylabel('Whisk Displacement','Fontsize',18)
    legend('Expected','Predicted','Fontsize',15);
    figNum=figNum+1;
end
figure(figNum);clf
subplot (1,2,2); 
plot (V,'o-');
title('VAF of EMG-Postion IRF','Fontsize',20);
xlabel('Segment','Fontsize',16);
ylabel('%VAF','Fontsize',16); 




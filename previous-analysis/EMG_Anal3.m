%% EMG Anal 3

%trials 3,4,5 are new trials, with animal behaving normally and whiksing
%well. EMG and whisker displacemnt recorded from right side of animal only.
%Analog notch filter was used for trials 1,2.  Analog notch filter NOT used
%for trials 3-5. Trial 3 and 4 were recorded from electrode pair 1
%(differential).  Trial 5 was recroded from a seperate implanted electrode
%pair (differential).  (The animal was implanted with two 2-channel
%epimysial electrode arrays into the right face, underlying the whisker
%musculature). 

clear Z Z1 EMG Whisk
load EMG_Whisk_SPONT3.mat
trial='5';
animal='13';


Fs = 1000;   % Sampling frequency after factor 10 decimation of EMG
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);   % Notch filter

           
EMG = decimate(eval(['H', animal, 'EMG', trial ]),10);
EMG = filtfilt(d,EMG);
Whisk = eval(['H',animal, 'Whisk',trial ]);

nFFT = 200;
Z=nldat(cat (2,EMG,Whisk),'domainIncr',0.001, 'chanNames', {'EMG' 'Whisk'});
Z=ddt(Z);
Z=Z-mean(Z);
%Z=detrend(Z);
%Z=smo(Z,10);

Z1=Z;
%Z1=Z1-mean(Z1(1:1000,:));
plot(Z1)
Z1(:,1)=abs(Z1(:,1));
figure(1);
plot(Z1);

 figure(2)
 subplot (2,2,1); 
 plot (pdf(Z(:,1))); 
 subplot(2,2,2);
 plot (pdf(Z(:,2))); 
 subplot (2,2,3);
 plot (spect(Z(:,1),'nFFT',nFFT));
  set(gca,'xlim',[0 400]); 
  subplot (2,2,4);
 plot (spect(Z(:,2),'nFFT',nFFT));
  set(gca,'xlim',[0 20]); 
 
 
%Zcur=Z1(iStart:iEnd,:);
figNum=3;
% Gain changes with amplitude so look at 10 s intervals. 
iNum=0;
gVal=[];
V = [];
for iStart=1:10000:100000;
    iNum=iNum+1;
    figure(figNum); clf
    iLen=10000;
    iEnd=iStart+iLen-1;
    Zcur=Z1(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags', 200);
    subplot(2,2,1);
    plot(i);
    subplot (2,2,2);
    gain=(cumsum(i)*.001);
    plot(gain);
    gVal(iNum)=mean(gain(150:200)); 
    title(['Gain = ' num2str(double(gVal(iNum)))])
    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre)
    V(iNum)=v;
    plot (Zcur(:,2));
    h=line(wPre);set(h,'color','r');
    title(['VAF= ' num2str(v)]);
    legend('Exp','Pre');
    figNum=figNum+1;
end
figure(figNum);clf
subplot (1,2,1); 
plot(-gVal,'o-')
xlabel('Segment');
ylabel('Gain'); 
title('Gain of EMG-Postion IRf');
subplot (1,2,2); 
plot (V,'o-');
title('VAF of EMG-Postion IRf');
xlabel('Segment');
ylabel('%VAF'); 


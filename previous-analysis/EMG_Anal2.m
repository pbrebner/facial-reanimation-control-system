%% EMG Anal

load EMG_Whisk_LongImplant.mat
trial='1';
eval(['Stim=H1Stimulus' trial ';']);
eval(['EMG=H1EMG' trial ';' ]);
eval(['Whisk=H1Whisk' trial ';' ]);

z=cat(2, Stim, EMG);
Zo=nldat(z,'domainIncr',1/10000,'chanNames', {'Stimulus' 'EMG'});
Z=decimate(Zo,10);
W=nldat(Whisk,'domainIncr',.001, 'chanNames', {'Whisk'});
Z=cat(2,Z,W);
Z1=Z(:,[2 3]);
Z1=Z1-mean(Z1(1:1000,:)); %Why just the mean from the first 1000 points?
plot(Z1)
Z1(:,1)=abs(Z1(:,1));
figure(1);
plot(Z1);

%Zcur=Z1(iStart:iEnd,:);
figNum=2;
% Gain changes with amplitude so look at 10 s intervals. 
iNum=0;
gVal=[];
V = [];
for iStart=1:5000:5000
    iNum=iNum+1;
    figure(figNum); clf
    iLen=4000;
    iEnd=iStart+iLen-1;
    Zcur=Z1(iStart:iEnd,:);
    Zcur=detrend(Zcur);
    i=irf(Zcur,'nLags',200);
    subplot(2,2,1);
    plot(i);
    subplot (2,2,2);
    gain=(cumsum(i)*.001);
    plot(gain);
    gVal(iNum)=mean(gain(150:200)); 
    title(['Gain = ' num2str(double(gVal(iNum)))])
    subplot (2,2,[3 4]);
    wPre=nlsim(i,Zcur(:,1));
    v=vaf(Zcur(:,2),wPre);
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






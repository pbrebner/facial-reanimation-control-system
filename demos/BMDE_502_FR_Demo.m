function mod10Demo ( inputCutOff, noiseLevel, nSamp, nFFT)
%%  Frequency Response  demos

if nargin < 1,
    inputCutOff=input_d('Input bandwidth',.1,0,1);
end
if nargin <2
    noiseLevel=input_d('Noise level', .1, 0, 100);
end
if nargin <3
    nLags=input_d('Number of lags for IRF', 50, 10,1000);
end
if nargin <4
    nSamp=input_d('Number of samples', 20000, 100,10^6);
end
nFFT=input_d('nFFT',500,100,5000);

nSides=2;
nLags=50;

delete (get(0,'children'));


%% Generate input signal
u=randn(nSamp,1);
if inputCutOff<1,
    [b,a]=butter(2,inputCutOff/2, 'low');
    u=filter(b,a,u);
    
end
figure(1);clf
U=nldat(u,'domainIncr',.005,'comment','Input signal');
subplot (2,1,1);
plot (U);
set (gca,'xlim',[0 10]);
ylabel('Input Signal');
fullTitle=['Input BW= ' num2str(inputCutOff) '; Noise:' num2str(noiseLevel) '; nSamp=' num2str(nSamp) ';NFFT:' num2str(nFFT)];
title(fullTitle);
subplot (2,1,2);
plot (spect(U));
ylabel('Spectrum');  title(''); 
set(gcf,'name',['Input signal BW=' num2str(inputCutOff)]);
figMod(1,'label_size',18,'lineWidth',1.5,'tick_label_size',16);




%% Static Linear
z=nlid_sim('static_linear',U,'noise_level', noiseLevel);
titleStr='static Linear' ;
plotFresp(2)
figMod(2,'label_size',18,'lineWidth',1.5,'tick_label_size',16);

%% No relation
z(:,2)=randn(length(z),1);
titleStr='No relation' ;
plotFresp(3);
figMod(3,'label_size',18,'lineWidth',1.5,'tick_label_size',16);


%% Dynamic LowPass
z=nlid_sim('L1',U,'noise_level', noiseLevel);
nSides=1; nLags=32
titleStr='Dynamic LowPass' ;
plotFresp(4)
title(''); 
figMod(4,'label_size',18,'lineWidth',1.5,'tick_label_size',16);

%% Dynamic LowPass with 60Hz noise;
z=nlid_sim('L1',U,'noise_level', noiseLevel);
t=domain(z);
noise60=sin(2*pi*60*t);
z(:,2)=z(:,2)+noise60;
titleStr='Dynamic LowPass with 60Hz noise' ;
plotFresp(5)
figMod(5,'label_size',18,'lineWidth',1.5,'tick_label_size',16);

%% Dynamic HighPass
z=nlid_sim('H2',U,'noise_level', noiseLevel);
nSides=2; nLags=16
titleStr='Dynamic high Pass' ;
plotFresp(6)

%% Static Linear with Delay
z=nlid_sim('static_linear',U,'noise_level', noiseLevel,'delay_time', .100);
nSides=1; nLags=50;
titleStr='Statc Linear with Delay' ;
plotFresp(7)
%% LowPass with Delay
z=nlid_sim('L1',U,'noise_level', noiseLevel,'delay_time', .100);
nSides=1; nLags=50;
titleStr='Low Pass with Delay' ;
plotFresp(8)

% low pass system with noise - large input 
figNum=8;
figNum=figNum+1; 
z=nlid_sim('L1',100*U,'noise_level', noiseLevel/10);
titleStr='Dynamic LowPass;Large input amplitude' ;
plotFresp(figNum)



%% nonlinear system 
z=nlid_sim('N3L',U,'noise_level', noiseLevel,'delay_time', .100);0;
titleStr='NL System; Low input amplitude' ;
figNum=figNum+1;figure(figNum); clf
plot(U); set(gca,'xlim',[0 25]);
title ('Low amplitude Input'); 
figMod(figNum,'label_size',18,'lineWidth',1.5,'tick_label_size',16);

figNum=figNum+1;
plotFresp(figNum)

figNum=figNum+1;
figure(figNum); clf
plot(8*U); set(gca,'xlim',[0 25]); 
title ('High amplitude Input'); 
figMod(figNum,'label_size',18,'lineWidth',1.5,'tick_label_size',16);



[z,M]=nlid_sim('N3L',10*U,'noise_level', noiseLevel,'delay_time', .100);0;
titleStr='NL System; Input*10' ;
figNum=figNum+1; plotFresp(figNum)


figNum=figNum+1; figure(figNum );
subplot(1,2,1);
P=M{1,1};
set(P,'polyRange',[-10 10]);
plot(P);
subplot (1,2,2); 
I=M{1,2};
plot(I)
streamer('Nonlinear System');
figMod(figNum,'label_size',18,'lineWidth',1.5,'tick_label_size',16);


function plotFresp(figNum)
z1=z;figure(figNum);clf;
fullTitle=['Input BW= ' num2str(inputCutOff) '; Noise:' num2str(noiseLevel) '; nSamp=' num2str(nSamp) ';NFFT:' num2str(nFFT)];
fullTitle=[titleStr '; ' fullTitle];
set (gcf,'name',fullTitle);
set(z1,'comment',fullTitle);
z1=z1-mean(z1); 
F=fresp(z1,'nFFT',nFFT);
set(F,'comment',fullTitle,'domainName','Hz')
plot(F);
figMod(figNum,'label_size',18,'lineWidth',1.5,'tick_label_size',16);

end
end



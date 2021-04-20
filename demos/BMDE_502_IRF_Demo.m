function BMDE_502_IRF_Demo (inputCutOff, noiseLevel, nLags, nSamp)
% mod9Demo (inputCutOff, noiseLevel, nLags, nSamp)
% Demonstrate IRF estimation for various systems
% inputCutOff - normalized cutoff for input singal (0-1)
% noiseLevel - ration of noise STD to output SRD (0 -100);
% nLags - lenght of IRF
% nSamp - number of samples.

%%  IRF demos

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
 
delete(get(0,'children'));
set(0,'DefaultFigureWindowStyle','docked') 

nSides=2;

%% Select bandwidth


%% Generate input signal
u=randn(nSamp,1);
if inputCutOff<1,
    [b,a]=butter(2,inputCutOff/2, 'low');
    u=filter(b,a,u);
end
figure(1); clf
U=nldat(u,'domainIncr',.01,'comment','Input');
subplot (2,1,1);
plot (U(1:1000));
title('');
ylabel('Input');
set(gca,'xtick',[0:2:10]);

subplot (2,1,2);
SU=spect(U,'nFFT',256);
set(U,'comment','Input Spectrum','domainName','sec'); 
plot (SU);
ylabel('Spectrum'); 
set(gca,'ytick',[0:.01:.03]) 
streamer(['Input signal BW=' num2str(inputCutOff)],.9);
figMod(1,'title_size',14,'lineWidth',2);
figure(2);
plotHessian(U);

% Generate noise signal
rNoise = randvar;
noise =nlsim(rNoise,domain(U));
noise=noise;

%% Static Linear
[z,m]=nlid_sim('static_linear',U,'noise_level',0);
nSides=2;
titleStr='Static Linear' ;
plotIRF(2)
%% Dynamic LowPass
z=nlid_sim('L1',U,'noise_level', 0);
nSides=1;
titleStr='Dynamic LowPass' ;
plotIRF(3)
%% Dynamic HighPass
z=nlid_sim('H1',U,'noise_level', 0);
nSides=2;
titleStr='Dynamic high Pass' ;
plotIRF(4)
%% Static Linear with Delay
z=nlid_sim('static_linear',U,'noise_level', 0,'delay_time', .100);
nSides=1;
titleStr='Statc Linear with Delay' ;
plotIRF(5)
%% LowPass with Delay
z=nlid_sim('L1',U,'noise_level', 0,'delay_time', .100);
nSides=1;
titleStr='Low Pass with Delay' ;
plotIRF(6)
%% 

z=nlid_sim('N3L',U,'noise_level',0);
nSides=1;
titleStr='Hammerstein System ' ;
plotIRF(7) 

function plotIRF(figNum)
figure(figNum);clf 
disp(titleStr);
fullTitle=[titleStr '; Input BW= ' num2str(inputCutOff) '; Noise:' num2str(noiseLevel) '; nSamp=' num2str(nSamp)];
set (gcf,'name',fullTitle);
zIn=z(:,1);acIn=cor(zIn,'nSides',2,'nLags',nLags);
% Add noise
stdZout=std(double(z(:,2)));
stdNoise=std(double(noise));
gain =noiseLevel*stdZout/stdNoise; 
zOut=z(:,2)+(noise*gain);
Z=cat(2,zIn,zOut); 
acOut=cor(zOut,'nSides',2,'nLags',nLags);
xCor=cor(Z,'nSides',nSides,'nLags',nLags);
I=irf(Z,'nSides',nSides,'nLags',nLags);
zPre=nlsim(I,zIn);
subplot (2,3,2);
plot (acIn);
title('');
xlabel('Input autocorrelation');

subplot (2,3,3);
plot (acOut);
ylabel('Output autocorrelation');
title(''); 
subplot (2,3,4);
plot (xCor);
ylabel('Cross Correlation');
title(''); 
subplot (2,3,5);
plot (I)
title('');
ylabel('IRF'); 
subplot (2,3,6);
plot (smo(I,3));
title(''); 
ylabel('smoothed IRF');

iVAF=double(vaf(zOut,zPre));
% Plot input, output, predicted output 
zPlot = cat(2,Z, zPre);
subplot (6,3,1);
plot (zPlot(50:550,1)); 
title (['VAF=' num2str(iVAF)]);
ylabel('Input');
xlabel(''); 
 
subplot (6,3,4);
plot (zPlot(50:550,2),'linecolor','r'); 
ylabel('Output');
title('');
xlabel(''); 

subplot (6,3,7);
plot (zPlot(50:550,3),'linecolor','g'); 
ylabel('Predicted');
title('');




streamer(fullTitle,.90);
figMod(figNum,'title_size',12,'linewidth',2);

end
end
function plotHessian(U)
        c=cor(U,'nLags',32);
        T=toeplitz(double(c));
        H=T'*T;
        mesh(H)
        condNum=cond(H);
        title(['Hessian. Condition number=' num2str(condNum)]);
;
    end
        
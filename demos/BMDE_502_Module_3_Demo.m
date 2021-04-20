function BMDE_502_Mod3_Demo (sysType, inputCutOff, noiseLevel, nLags, nSamp)
% mod9Demo (inputCutOff, noiseLevel, nLags, nSamp)
% Demonstrate IRF estimation for various systems
% inputCutOff - normalized cutoff for input singal (0-1)
% noiseLevel - ration of noise STD to output SRD (0 -100);
% nLags - lenght of IRF
% nSamp - number of samples.

%%  IRF demos
if nargin < 2,
    sysTypeList = { 'LP1' 'HP1'  'L1' 'L2'  'LNRect' 'LN2' 'LN3' 'N3L' 'N2L' ...
        'LNL' 'PC' 'POLY' 'Static_Linear' 'Cuber' 'N3HP'};
    i=menu('Option',sysTypeList);
    sysType=sysTypeList{i};
end
if nargin < 2,
    inputCutOff=input_d('Input bandwidth',.1,0,1);
end
if nargin <3
    noiseLevel=input_d('Noise level', .1, 0, 100);
end
if nargin <4
    nLags=input_d('Number of lags for IRF', 50, 10,1000);
end
if nargin <5
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

%%
[z,m]=nlid_sim(sysType,U,'noise_level',0);
plotIRF(3,m)
plotFresp(10,m);

    function plotIRF(figNum,m)
        figure(figNum);clf
        titleStr=sysType;
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
        %% Plot IRF and correlation functions
        subplot (2,3,1);
        plot (acIn);
        title('Input autocorrelation');
        subplot (2,3,2);
        plot (acOut);
        title('Output autocorrelation');
        
        subplot (2,3,3);
        plot (xCor);
        title('Cross Correlation');
        subplot (2,3,4);
        plot (I)
        title('IRF Estimate');
        subplot (2,3,5);
        plot (smo(I,3));
        title('smoothed IRF')
        subplot (2,3,6);
        plot(m)
        
        title('True system ');
        
        iVAF=double(vaf(zOut,zPre));
        % Plot input, output, predicted output
        f=figNum+1;
        figure(f);clf
        zPlot = cat(2,Z, zPre);
        zResid=Z(:,2)-zPre;
        
        subplot (2,2,1);
        plot (zPlot(50:550,1));
        
        title('Input');
        xlabel('');
        
        subplot (2,2,3);
        plot (zPlot(50:550,2),'linecolor','r');
        title('Output');
        xlim=get(gca,'xlim');
        xlabel('');
        
        subplot (2,2,2);
        plot (zPlot(50:550,3),'linecolor','g');
        title(['Predicted VAF=' num2str(iVAF)]);
        subplot (2,2,4);
        plot (zResid,'linecolor','g');
        set(gca,'xlim',xlim);
        title('Residual');
        streamer('Signals' ,.90);
        
        %% Propeerties of Residuals
        
        f=f+1;
        figure(f); clf
        
        % Autocorrelation
        subplot (2,2,1);
        plot (zResid);
        title('Residual');
        rAC=cor(zResid);
        
        subplot (2,2,2);
        plot (rAC)
        title('AytoCorrelation');
        % Power Spectrum
        rSpect=spect(zResid);
        subplot (2,2,3);
        plot (rSpect);
        title('Power Spectrum');
        % PDF
        rPDF=pdf(zResid,'pdfType','Frequency');
        subplot (2,2,4) ;
        plot (rPDF);
        title('PDFl');
        streamer('Properties of Residuals' ,.90);
        
        %% Different options for IRF etimation
        i2=irf(Z,'nLags',nLags,'irfIdMethod','corr');
        i3=irf(Z,'nLags',nLags,'irfIdMethod','pseudo','irfPseudoInvMode','full' );
        i4=irf(Z,'nLags',nLags,'irfIdMethod','pseudo','irfPseudoInvMode','auto' )
        
        i5=irf(Z,'nLags',nLags,'irfIdMethod','pseudo','irfPseudoInvMode','manual' ,'irfFigNum',99);
        f=f+1;
        figure(f);clf
        subplot (2,2,1); plot(i2);title('Correlation Method');
        subplot (2,2,2); plot(i3);title('Full pseduoinverse');
        subplot (2,2,3); plot(i4);title('Pseduoinverse - automatic');
        subplot (2,2,4); plot(i5);title('Pseduoinverse -manual');;
        streamer ('Diffeernt IRF options');
        ;
        
        
        
        
        figMod(figNum,'title_size',12,'linewidth',2);
    end

    function plotHessian(U)
        c=cor(U,'nLags',32);
        T=toeplitz(double(c));
        H=T'*T;
        mesh(H)
        condNum=cond(H);
        title(['Hessian. Condition number=' num2str(condNum)]);
    end;
    
    %%
    function plotFresp(figNum,m)
        nFFT=input_d('nFFT',500,100,5000);
        z1=z;
        % vAdd noise
        stdZout=std(double(z(:,2)));
        stdNoise=std(double(noise));
        gain =noiseLevel*stdZout/stdNoise;
        zIn=z(:,1);
        zOut=z(:,2)+(noise*gain);
        Z=cat(2,zIn,zOut);
        Z=Z-mean(Z);
        F=fresp(Z,'nFFT',nFFT);
        set(F,'comment','Frquency Response'); 
        zPre=nlsim(F,zIn);
        
        figure(figNum);clf;
        fullTitle='Fequency Response';
        fullTitle=[sysType '; ' fullTitle];
        set (gcf,'name',fullTitle);
        set(z1,'comment',fullTitle); 
        set(F,'comment',fullTitle,'domainName','Hz')  
        figMod(figNum,'label_size',18,'lineWidth',1.5,'tick_label_size',16);
        
        f=figNum+1;figure(f); clf
       
        subplot (3,1,1); 
        plot (zOut); 
        title('Ooutput'); 
        subplot (3,1,2);
        plot (zPre); 
        title('FRESP predicted output');
        subplot (3,1,3);
        plot (zOut-zPre);
        title('Residuals'); 
        figMod(f); 
        
    end
end

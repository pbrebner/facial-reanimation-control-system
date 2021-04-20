function tvmDemo
% tvmDemo - demonstrate tvm propoerties
%% TV IRF
%%Generate an TV IRF
figNum=1;
disp ('tvIRF identification');
disp('Generate TV IRF');
gain=ones(200,1);
gain=cat(1,gain,linspace(0,5,600)'); %Concatenate arrays
gain=cat(1,gain, ones(200,1)); %cat(DIM, A, B), cat(1,A,B) is the same as [A;B]
plot(gain)
tvI={};
 for i=1:1000
     tvTemp=irf2(irf,'g',gain(i)); %irf2 - Generate a second order IRF with Gain, Damping and Natural Frequency
     %irf - impulse response function class for NLID toolbox
     tvI{i}=tvTemp;
 end
 TVirf=tvm; %Time varying Model Class
set(TVirf,'elements',tvI,'tvStart',0,'tvIncr',.01);
figNum=figNum+1; figure(figNum);clf;
plot(TVirf); title('Simuulated TV IRF'); 
%% simulate TV responses
disp('Simulate at TV Response');
x=randn(1000,1,200);
X=nldat(x,'domainIncr',.01);%data class for NLID toolbox
Y=nlsim(TVirf,X); %Simulate response of tvm to input data set
Z=cat(2,X,Y); %cat(2,A,B) is equal to [A,B]
set(Z,'chanNames',{ 'Input' 'OutPut'},'comment','Simulated TV trial');
figNum=figNum+1; figure(figNum);clf;

plot(Z);

%% Estimate a TI IRF
disp('TI IRF');
figNum=figNum+1; figure(figNum);clf;
I=irf; %impulse response function class for NLID toolbox
set(I,'nLags',100);
tiIRF=nlident(I,Z);
[R,V,yp]=nlid_resid(tiIRF,Z); %Compute and Display Prediction Error in model Output
%R - Residuals, V - Variance Accounted For, yp - Predicted Output
figNum=figNum+1; figure(figNum);clf;
subplot (2,1,1);
plot(V);%VAF (y-axis is percentage)
title('TI IRF VAF');
xlabel('Realization');
subplot (2,1,2);
plot(R); %Residuals (Plot shows that they are time varying)
title('TI Residuals');
xlabel('Time (s)');



%% Estimate a TV IRF using ensemble method 
disp('TV IRF Ensemble method');
tvIRFensemble = tvm; %time varying model class
set(tvIRFensemble,'tvIdentMethod','ensemble'); 
tvIRFensemble=nlident(tvIRFensemble, Z, I); %nlident method supports identification of: irf, polynom, and nl block models
figNum=figNum+1; figure(figNum);clf;

plot(tvIRFensemble);
figNum=figNum+1; figure(figNum);clf;

title('TV IRF ensemble estimate');
figNum=figNum+1; figure(figNum);clf;

tvResid(tvIRFensemble,Z);%function tvResid (V = tvResid(tvModel, Z))


%% Estimate a TV IRF using basis expansion 
disp('Estimate a TV IRF using basis expansion') 
tvIRFbasis = tvm;%time varying model class
set(tvIRFbasis,'tvIdentMethod','basisexpansion'); 
BF=polynom; %polynomial class for NLID toolbox
set(BF,'polyType','Bspline','polyOrder',10, 'polyRange',[ 0 10],'splineSD',1,'splineCenters',0.5:1:10);
nodeLocations = .5:1:10;
set(BF,'polyCoef',nodeLocations(:)); 
%Length of spline Centers must match up with length of nodeLocations
t=domain(Z);
BF1=basisfunction(BF,t);%Returns basis functions for polynomial P evaluated over x
plot(BF1); 
figNum=figNum+1; figure(figNum);clf;
tvIRFbasis=nlident(tvIRFbasis, Z(:,:,1:50), I,BF,'periodic','yes','method','Bayes');
%nlident method supports identification of: irf, polynom, and nl block models
figNum=figNum+1; figure(figNum);clf;

plot(tvIRFbasis);
title('TV IRF basic function  estimate');
figNum=figNum+1; figure(figNum);clf;
tvResid(tvIRFbasis,Z);%function tvResid (V = tvResid(tvModel, Z))
subplot (3,1,1); title('Residuals for TV IRF Basis Function Estimate');


%% TV Polynomial
disp('TV polynomial demo');
% Generate a TV Polynomial
p=polynom('polyType','power','polyOrder',2,'polyOrderMax',3) %polynomial class for NLID toolbox
PI=[];
coef=[ 0 1 1];
for i=1:1000,
    coef(2)=1+gain(i);
     pTemp=set(p,'polyCoef',coef);
     PI{i}=pTemp;
end
 TVpoly=tvm; %time varying model class
 set(TVpoly,'elements',PI,'tvIncr',X.domainIncr);
figNum=figNum+1; figure(figNum);clf;

 plot(TVpoly);
 title('TV Polyomial'); 
 % Simulate response to TV polynomial
Ypoly=nlsim(TVpoly,X); %simulate response of tvm to input data set
Zpoly=cat(2,X,Ypoly); %cat(2,A,B) is equal to [A;B]
figNum=figNum+1; figure(figNum);clf;
plot(Zpoly);
title('Simulated TV polynmial response'); 


% Identify a TV polynomial
TVpolyIdent=nlident(tvm,Zpoly,p); %identify
figNum=figNum+1; figure(figNum);clf;
plot(TVpolyIdent);
title('Estimated TV Polynomial'); 
%
simY=nlsim(TVpolyIdent,X); %simulate response of tvm to input data set
figNum=figNum+1; figure(figNum);clf;
subplot(3,1,1);
plot(Y);
title('Output');
subplot (3,1,2);
plot (simY);
title('Estimated Output');
subplot(3,1,3);
plot(Y-simY);
title('Residuals');
%% TV Hammerstein Identification - ensemble method
disp('TV Hammerstein  - In progres');
polyElements=TVpoly.elements;
irfElements=TVirf.elements;
nElements=length(polyElements);
NLBL=nlbl; %nonlinear - linear block model for NLID toolbox
%  Note that the length of the irf and the maximum order of the
%  nonlinearity are parameters of the nl class.

for i=1:1000
    set(NLBL,'elements',{polyElements{i} irfElements{i}});
    nlblElements{i}=NLBL;
end
TVnlbl=tvm; %time varying model class
set(TVnlbl,'elements',nlblElements,'tvStart',0,'tvIncr',.01);
y1=nlsim(TVirf,Ypoly); %simulate response of tvm to poly input data set

yNLBL=nlsim(TVnlbl,X); %simulate response of tvm to input data set

Znlbl=cat(2,X,yNLBL); %cat(2,A,B) is equal to [A;B]

TVnlblIdent=nlident(tvm,Zpoly,NLBL); %Identify
yPre=nlsim(TVnlblIdent, X); %simulate response of tvm to input data set

figNum=figNum+1; figure(figNum);clf;
plot(y1)

figNum=figNum+1; figure(figNum);clf;
plot(yNLBL)

figNum=figNum+1; figure(figNum);clf;
plot(Znlbl)

figNum=figNum+1; figure(figNum);clf;
plot(yPre)












   
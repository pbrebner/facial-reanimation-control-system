%% Module 2 Demonstration
%% QR Decompompostion
A=rand(5,3);
[Q,R]=qr(A);
disp('A');
disp(A);
disp('Q');
disp(Q)
disp('R')
disp(R)

%% Singular Value Decompoition

A=[1 2; 3 4; 5 6; 7 8];
A=abs(rand(4,4));

[U,S,V]=svd(A);
A
U
S
V
%% Reconstruct
A
Apre=U*S*V'

%% Orthogonal properties of U*U'

UUT=U*U'

VVT= V*V'

%% SECTION TITLE
% GAUSSIAN VARIABLES
%% randvar - nlid class supporting random variates
r=randvar;
disp(r)
% See what types are available
set(r,'randvarType');
set(r,'randvarType','Normal');disp(r);
% pdf supports randvar objects as well. 
clf;plot (pdf(r));


%% Generate a realization of 1000 points
% nlsim method of ranvar generates a realization of a randvar object
X=nlsim(r,[1:1000]')
set(X,'domainIncr',.001);
plot (X)
%% compare theoretical and observed PDFs
plot (pdf(X,'nBins',20))
px=pdf(r);  % Compute theoretical distribution
h=line(px);set(h,'color','r','linewidth',2)
legend ('data','theoretical');
%% Correlation Functions
X=nlsim(r,100000);Y=smo(X,10); 
%% Auto correlation fuction
figure(1);clf
c=cor(Y,'nSides',2,'nLags',25);
subplot (2,1,1);
plot (c);
title('First order correaltion function for linear signal '); 
subplot (2,1,2); 
c=cor(Y,'nSides',2,'nLags',25,'kernOrder',2);
plot (c);
title('Second order auto-correlation function for linear signal '); 
figMod(1)
%% Nonlinear 
figure(1);clf

Z=X.^2;
c=cor(Z,'nSides',2,'nLags',25);
subplot (2,1,1);
plot (c);
title('First order correaltion function for nonlinear signal'); 
subplot (2,1,2); 
c=cor(Z,'nSides',2,'nLags',25,'kernOrder',2);
plot (c);
title('Second order correlation function for nonlinear signal'); 
figMod(1)
%% Cross correlation 
figure(1); clf
Y=smo(X,10);
Z=cat(2,X,Y); 
c=cor(Z,'nSides',2,'nLags',25);
subplot (2,1,1);
plot (c);
title('First order cross-correaltion function for linear filter'); 
subplot (2,1,2); 
c=cor(Z,'nSides',2,'nLags',25,'kernOrder',2);
plot (c);
title('Second order correaltion function for linear filter '); 
%% Nonlinear filter
figure(1); clf
Y=smo(X,10);
Z=cat(2,X,Y); 
c=cor(Z,'nSides',2,'nLags',25);
subplot (2,1,1);
plot (c);
title('First order cross-correaltion function for linear filter'); 
subplot (2,1,2); 
c=cor(Z,'nSides',2,'nLags',25,'kernOrder',2);
plot (c);
title('Second order correalation function for linear filter '); 

%% Power Spectrum 
% SPectrum of a white signal
S=spect(X,'nFFT', 200);
figure(1); clf;
plot(S);

%% Low pass filter 
Xlow=smo(X,10);
S=spect(Xlow,'nFFT', 200);
figure(1); clf;
plot(S);


%% Power Series Polynomials
p=polynom;
set(p,'polyType','power');
pb = basisfunction(p);
pb.comment=[];
plot(pb);
streamer('Power series basis functions'); 
figMod(1);

%% Fit to a half wave rectifier
X=nlsim(r,1000);X.domainIncr=.001;
Y=max(0,X); 
Z=cat(2,X,Y);




%% Hermite  Series Polynomials
p=polynom;
set(p,'polyType','hermite');
pb = basisfunction(p);
pb.comment=[];
figure(1); clf;
plot(pb);
streamer('Hermite  series basis functions'); 
figMod(1);

%% Techebyshev Series Polynomials
p=polynom;
set(p,'polyType','tcheb');
pb = basisfunction(p);
pb.comment=[];
figure(1); clf;
plot(pb);
streamer('Tchebyshev series basis functions'); 
figMod(1);


%% 










% whisk
%load baseline_whisking4
%load ('baseline_whisking4','darb')
% Read data, create Nnldat object, detrend and decimate;
% 'a" #1-20 is a pre-surgical group of normal animals
% 'a' #21-28 are post-surgical (surgery on left side.. much weaker whisking
% after surgery)
% a2, 10 15, 22, and 25 seem to be good data
% 'b' is another group of pre-treatment normal animals (L and R sides)
% (1-20 all normal)
% b 3,15,18 seem to be pretty good data
% abes is an experiment looking for effect of brief electrical stimulation
% 'bes' on recovery.  Pre-operative data ('pre') is normal state,
% post-operative data (POD7-105) is data after transection and repair
% injury.  Both whisking and blinking data is available for these animals.
% animal 15 and 18 is a good example
% darb is an experiment with animals who received and did not receive
% darbopoeitin

%dataset = {'a'};
%iTrial = 29;  % Set trial to visualize / display
group = ['abes'];
subgroup = ['abes'];
animal = 17; 
trial = ['pre'];   
% pre or POD7...14,21,28,35,42,49,56,77,92...POD105  (POD = post-operative
% day) for abes
% pre or POD 6,7,8,8,10,11,12,13,14,15,16,17,18,19,20,21,28,36,62 for darb

%PROBLEM:  abes.abes18POD14 .... chatter?

func1 = ['Whisk'];  % Blink or Whisk
func2 = ['Whisk'];
side1 = ['L'];  % R or L
side2 = ['R'];

nFFT=200;
% Set frequency range of bandpass filter BPfilt
%fLow = 4;
%fHigh = 5;
%PB1 = fLow;
%PB2 = fHigh;
DF = 10;   % Set decimate factor
%Fs = 1000/DF;
%bpF = BPfilt(PB1,PB2,Fs);
%GainAvg = zeros(L,1);
%PhaseAvg = zeros(L,1);
%CohereAvg = zeros(L,1);

str1 = [group '.' subgroup, num2str(animal), '_', func1, trial, side1];
str2 = [group '.' subgroup, num2str(animal), '_', func2, trial, side2];
%strL = [dataset ; num2str(iTrial) ; 'L'];
%strR = [dataset ; num2str(iTrial) ; 'R'];
%fL=filtfilt(bpF,decimate(eval([strL{:}]),DF));
%fR=filtfilt(bpF,decimate(eval([strR{:}]),DF));
f1=decimate(eval(str1),DF);
f2=decimate(eval(str2),DF);
xStop=5000; 
%f1 = f1(1:xStop);
%f2 = f2(1:xStop);
%fL=decimate(eval([strL{:}]),DF);
%fR=decimate(eval([strR{:}]),DF);
s1 = [func1 '-' side1];
s2 = [func2 '-' side2];
W = nldat(cat(2,f1,f2),'domainIncr',.01, 'chanNames',{s1 s2});
W=ddt(W);
W=W-mean(W);
W=detrend(W); 
W=ddt(W);

 figure(1); clf
 W=smo(W,10);
 plot(W+[0 500])
 footer(['Animal: ' num2str(animal), '-', trial])
 
 
 % PDF and spectra
  figure(2)
 subplot (2,2,1); 
 plot (pdf(W(:,1))); 
 subplot(2,2,2);
 plot (pdf(W(:,2))); 
 subplot (2,2,3);
 plot (spect(W(:,1),'nFFT',nFFT));
  set(gca,'xlim',[0 20]); 
  subplot (2,2,4);
 plot (spect(W(:,2),'nFFT',nFFT));
  set(gca,'xlim',[0 20]); 
 footer(['Animal: ' num2str(animal), '-', trial])

 figure(3);  clf
 plot(W,'plotmode','xy') 
 vW=vaf(W);
cW= corrcoef(W)
titleStr= ['R = ' num2str(cW(1,2)) '; VAF = ' num2str(double(vW))]

p=polynom(W,'polyType','power','polyOrderMax',6);
hold on 
w1=double(W(:,1)); 
xAx=(linspace(min(w1),max(w1)))';
pSim=nlsim(p,xAx);
h=line(xAx,double(pSim));set(h,'color','r');
[R,vP]= nlid_resid(p,W,'plotFlag',false);

titleStr= ['R = ' num2str(cW(1,2)) '; static VAF = ' num2str(double(vW)) ';  polyVAF = ' num2str(double(vP))]

title(titleStr);
 footer(['Animal: ' num2str(animal), '-', trial])

% IRF 
figure(4);clf
i=irf(W,'nLags',32,'nSides',2,'irfPseudoInvMode','auto');
[R,V,yp]= nlid_resid(i,W,'plotFlag',false);
subplot (2,1,1);
plot(i)
subplot (2,1,2);
plot(spect(W(:,2)));
hold on;
h=line(spect(yp,'nFFT',nFFT));
title(['VAF = ' num2str(double(V))]);
set(h,'color','red');
set (gca,'xlim',[0 20]);
 footer(['Animal: ' num2str(animal), '-', trial])


% Frequency Response

figure(5);
F=fresp(W,'nFFT',nFFT);
plot(F(1:40,:));
[R,vF]= nlid_resid(F,W,'plotFlag',false);
 footer(['Animal: ' num2str(animal), '-', trial])

% orthogonal regeression
figure(6);clf
WF=(smo(W,10));
plot(WF,'plotmode','xy');
w=double(WF); 
[c,s,l]=pca(w);
slope=c(2,1)/c(1,1);
xGrid=linspace(min(w(:,1)),max(w(:,1)),50);
xFit= s(:,1)*c(:,1)';
h=line (xGrid,slope*xGrid); 




%{
G2 = gain (F );
xAx = domain(G2);
GdB =20*log10(double(G2));
P2 = phase(nldat(F ));
Pdeg =180*(unwrap(double(P2),[],1)/pi);
C2 = coherence(F);

figure(6);
subplot(3,1,1)
xAx = domain(G2);
GdB =20*log10(double(G2));
plot(xAx, GdB);
xlim([0 20]);
title(['Trial ',num2str(iTrial),' - Gain']); 
ylabel('Gain (dB)');
subplot(3,1,2)
plot(xAx, Pdeg);
xlim([0 20]);
ylabel('Phase (deg)');
title('Phase'); 
subplot(3,1,3)
plot(C2);
axis([0 20 0 1]);
ylabel('Coherence');

 %}

%%plot(eval([strL{:}]));



%{
%%  Segment Analysis
 nSamp=length(W);
 nSeg=50;
 nSampNew=nSamp/nSeg;
 clear S
 clear CC
 
S=nldat;
II=nldat;
RR=nldat; 
 for i=1:nSeg
     j1=1+(i-1)*nSampNew;
     j2=j1+nSampNew-1; 
     W1=(W(j1:j2,:));
     c=corrcoef(W1);
     s=std(W1);
     S(i,:)=s;      
     CC(i)=c(1,2);
     ii=irf(W1,'nLags',32,'nSides',2,'irfPseudoInvMode','auto');
[R,V]= nlid_resid(ii,W1,'plotFlag',false);
  F=fresp(W1,'nFFT',nFFT);
if i==1;
    II=ii;
    FF=F;
else
    II=cat(3,II,ii);
    FF=cat(3,FF,F);
end
VV(i)=V;
 end
  

%}edit synkinesis
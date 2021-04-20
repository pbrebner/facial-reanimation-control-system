load WhiskEMG.mat
f1 = EMG5(1:300000);
f2 = LWhisk5;
s1 = 'EMG';
s2 = 'WhiskDisplacement';
nFFT = 100;

W = nldat(cat(2,f1,f2),'domainIncr',.01, 'chanNames',{s1 s2});
W=ddt(W); % 1st differential - velocity
W=W-mean(W);  % Subtracting the mean velocity... is this necessary?
W=detrend(W);  % Unclear what detrend actually does... is it necessary?
W=ddt(W);   % 2nd differential .... acceleration.  Is it necessary?
%W=decimate(W,DF);   % smo makes a really makes a BIG difference in VAF
F=fresp(W,'nFFT',nFFT);
C = coherence(F);

% xAx = domain(C);
% idxfLow = find(xAx == fLow);
% idxfHigh = find(xAx == fHigh);
% CMax(n,1) = max(C(idxfLow:idxfHigh));
    
 figure(1); clf
 % W=smo(W,10);
 plot(W+[0 500])
% footer(['Trial: ' num2str(iTrial)])
 
 
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
% footer(['Trial: ' num2str(iTrial)])

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
% footer(['Trial: ' num2str(iTrial)])

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
set(h,'color','red');
set (gca,'xlim',[0 20]);
% footer(['Trial: ' num2str(iTrial)])


title(['VAF = ' num2str(double(V))]);
%footer(['Trial: ' num2str(iTrial)])
% Frequency Response

figure(5);
F=fresp(W,'nFFT',nFFT);
plot(F(1:40,:));
[R,vF]= nlid_resid(F,W,'plotFlag',false);
 %footer(['Trial: ' num2str(iTrial)])


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
title(['Trial - Gain']); 
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

figure(7)
subplot(2,1,1)
plot(f1);
subplot(2,1,2)
plot(f2);



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
 
 
 
 %{

%% Gain, Phase, Coherence

Gain = double(gain(F));



            
            
Phase = double(phase(F));
Coherence = double(coherence(F));
freq = abs(f);


f=nldat(F); 
f=F.dataSet;
f = f(:,1);

f1=gain(F); %Gain
f1=20*log10(double(f1)); 
fPhase=phase(F); 
f2=180*(unwrap(double(fPhase),[],1)/pi);
f3 = coherence(F );

figure
subplot(3,1,1)
plot(f,f1);
subplot(3,1,2)
plot(f,f2);
subplot(3,1,3)
plot(f,f3);
     
%}
%}

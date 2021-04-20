load video.mat

group = ['video'];
subgroup = ['trial'];
trial = 4;   %trial 4-8

func1 = ['Whisk'];  % Blink or Whisk
func2 = ['Blink'];
side1 = ['L'];  % R (normal) or L (operated)
side2 = ['L'];

nFFT=200;
DF = 1;   % Set decimate factor
BL = 100; % set length of blink in ms

str1 = [group '.' subgroup, num2str(trial), '_', func1, side1];
str2 = [group '.' subgroup, num2str(trial), '_', func2, side2];
f1=decimate(eval(str1),DF);
f2=decimate(eval(str2),DF);
f1 = diff(f1);
f2 = diff(f2);
idx = find(f2 > 7*std(f2));    % > 5 standard deviations  identifies blink

for n = 1:(length(idx)-1);
    if idx(n+1) - idx(n) <= 5
        idx(n+1) = 0;
    else
    end
 end
i2 = find(idx == 0); 
idx(i2) = [];

BL1 = BL +1;
F1 = zeros(BL1*length(idx),1);
F2 = zeros(BL1*length(idx),1);

for n = 1:length(idx)
    A = f1(idx(n)-BL/2:idx(n)+BL/2);
    B = f2(idx(n)-BL/2:idx(n)+BL/2);
    F1(n*BL1-BL:n*BL1) = A;
    F2(n*BL1-BL:n*BL1) = B;
end

s1 = [func1 '-' side1];
s2 = [func2 '-' side2];
W = nldat(cat(2,F1,F2),'domainIncr',.01, 'chanNames',{s1 s2});
W = ddt(W);
W = detrend(W);


 figure(1); clf
 %plot(W+[0 500])
 plot(W)
 footer(['Trial: ' num2str(trial)])
 
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
 footer(['Trial: ' num2str(trial)])

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
 footer(['Trial: ' num2str(trial)])

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
 footer(['Trial: ' num2str(trial)])


% Frequency Response

figure(5);
F=fresp(W,'nFFT',nFFT);
plot(F(1:40,:));
[R,vF]= nlid_resid(F,W,'plotFlag',false);
 footer(['Trial: ' num2str(trial)])


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
figure(7)
subplot(3,1,1)
p1 = decimate(eval(str1),DF);
p1 = p1(1:6000);
plot(p1, 'r');
%plot(eval([strR{:}]));
subplot(3,1,2)
p2 = decimate(eval(str2),DF);
p2 = p2(1:6000);
plot(p2, 'b');
subplot(3,1,3)
plot(p1, 'r');
hold on;
plot(p2, 'b');
footer(['Trial: ' num2str(trial)])
hold off;
%plot(eval([strL{:}]));



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
  

%}
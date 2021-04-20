load baseline_whisking
% Read data, create Nnldat object, detrend and decimate;
% 'a" #1-20 is a pre-surgical group of normal animals
% 'a' #21-28 are post-surgical (surgery on right side.. much weaker whisking
% after surgery)
% a2, 10 15, 22, and 25 seem to be good data
% 'b' is another group of pre-treatment normal animals (L and R sides)
% (1-20 all normal)
% b 3,15,18 seem to be pretty good data

dataset = {'a'};
%iTrial = 1;  % Set trial to visualize / display
%Trials = [1,3,15,18,5,10,13];  % good 'b' trials - normals
Trials = [2,3,4,10,12,13,14,15,16];  % good 'a' trials  - normals
L = length(Trials);
nFFT=200;
DF = 10;   % Set decimate factor

iTrial = 2;
strL = [dataset ; num2str(iTrial) ; 'L'];
strR = [dataset ; num2str(iTrial) ; 'R'];
figure(2)
subplot(2,1,1)
plot(eval([strR{:}]));
subplot(2,1,2)
plot(eval([strL{:}]));


rmsRatio = zeros(L,1);
PowerRatio = zeros(L,1);
static_polyfitVAF = zeros(L,1);
dynamicVAF = zeros(L,1);
cohere = zeros(101, L);

for n = 1:L;
iTrial = Trials(n);
strL = [dataset ; num2str(iTrial) ; 'L'];
strR = [dataset ; num2str(iTrial) ; 'R'];
%fL=filtfilt(bpF,decimate(eval([strL{:}]),DF));
%fR=filtfilt(bpF,decimate(eval([strR{:}]),DF));
fL=decimate(eval([strL{:}]),DF);
fR=decimate(eval([strR{:}]),DF);
W = nldat(cat(2,fL,fR),'domainIncr',.01, 'chanNames',{'LEFT' 'RIGHT'});
W=ddt(W);
W=W-mean(W);
W=detrend(W); 
W=ddt(W);
W=smo(W,10);    %  This GREATLY increases the VAF! (By ~ 10%)
rmsL = rms(fL);
rmsR = rms(fR);
rmsRatio(n) = rmsR/rmsL;
p=polynom(W,'polyType','power','polyOrderMax',6);
[R,vP]= nlid_resid(p,W,'plotFlag',false);
static_polyfitVAF(n) = double(vP);
i=irf(W,'nLags',32,'nSides',2,'irfPseudoInvMode','auto');
[R,V,yp]= nlid_resid(i,W,'plotFlag',false);
dynamicVAF(n) = double(V);
F=fresp(W,'nFFT',nFFT);
cohere(:,n) = coherence(F);
end

xAx = domain(gain(F));
cohereAVG = mean(cohere')';
figure(2)
plot(xAx, cohere);
hold on;
plot(xAx, cohereAVG, '--r', 'LineWidth',3);
ylabel('Coherence');
title('Coherence Plots'); 
axis([0 20 0 1]);
xlabel('Frequency (Hz)');


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
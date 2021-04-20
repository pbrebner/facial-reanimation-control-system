% Whisk_segments  - stationary analsyis for whisk data

figure(1); clf;
plot(W);

% Divide data  into segments
segLen=1000;
[nsamp,nchan,nreal]=size(W);
nReal=nsamp/segLen;
W1=reshapeChan(W,segLen,nReal);

figure(2); clf
% Variance
V=var(W1);
subplot (2,1,1);
stem(V');
ylabel('Variaince');
subplot(2,1,2);
cc=corrcoef(W1);
cc=squeeze(cc(1,2,:));
stem(cc);
ylabel('Correlation Coef');
xlabel('Segment Number'); 

%% IRF

I=irf(W1,'nLags',16,'nSides',2,'irfPseudoInvMode','auto');
[resid, v, yp]=nlid_resid(I, W1);
figure(3); clf
stem(v); 
xlabel('Segment Number');
ylabel('IRF %VAF');








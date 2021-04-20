%% TVModelDemo
% Demonstrationt of various aspects of TV identification for BMDE 502
%
% Use cell mode and run ech section separately.

T = 1/100; %Sampling time
duration= 10; % seconds
t = 0:T:duration-T; %Time vector
N= length(t);

%Input design 
[fil_n,fil_d]=butter(2,20/((1/T)/2),'low');  %2nd Order 20Hz low-pass filter.

Input = filtfilt (fil_n, fil_d,0.1*randn(N,1));
%input is low-pass filtered, white-Gaussing signal
d_Input = gradient(Input,T); %first derivative input
dd_Input = gradient(d_Input,T); %second derivative i
%%
%TI Model parameters 
K = 20; %Elasticity
B = 0.25; %viscosity 
I = 0.002; %inertia

%simulated ouput
Output = K*Input + B*d_Input + I*dd_Input;

%plots of input and ouput
figure(1);
subplot(3,1,1)
plot(t,Input,'LineWidth',2,'Color',[1 0 0])
title('Input')
ylabel('x(t)')
subplot(3,1,2)
plot(t,Output,'LineWidth',2,'Color',[0 0 1])
title('Output')
xlabel('Time(s)')
ylabel('y(t)')
subplot(3,1,3)
plot([0 10],[K K],'LineWidth',2,'Color','g')
title('K')
xlabel('Time(s)')
ylabel('K(t)')
streamer('Time Invariant K');


%%
%TV Model parameters 
K = [30*ones(N/4,1); linspace(30,10,N/2)';10*ones(N/4,1)]; % TV Elasticity
B = 0.25; %viscosity 
I = 0.002; %inertia

%simulated ouput
Output = K.*Input + B*d_Input + I*dd_Input;
%note that we use the notaion .* to signifiy item by item multiplication

%plots of input and ouput
figure(1);clf
subplot(3,1,1)
plot(t,Input,'LineWidth',2,'Color',[1 0 0])
title('Input')
ylabel('x(t)')
subplot(3,1,2)
plot(t,Output,'LineWidth',2,'Color',[0 0 1])
title('Output')
ylabel('y(t)')
subplot(3,1,3)
plot(t,K,'LineWidth',2,'Color',[0 1 0])
title('K'); 
streamer('Time-Varying Elasticity')
ylim([5,35])
ylabel('K(t)')
xlabel('Time(s)')
%%  
% We can use matlab's fitlm (fit linear model) function to find the model 
% parameters, other options are avaliable, such as '\' operator. They should provide 
% similar results. 

%estimating a TI system using TV data
tqEst=[];
kEst=[];
TI_model = fitlm([Input d_Input dd_Input],Output,'y ~ x1+x2+x3-1')
tqEst=TI_model.Fitted;
kEst=table2array(TI_model.Coefficients(1,1));

%estimate the models paraters using matlab's fitlm function 
%passing the command 'y ~ x1+x2+x3-1' tells the function that they are three
%inputs and that no constant value should be used

%plots of input and ouput
%plots of input and ouput
figure(2);
subplot(4,1,1)
plot(t,Output,'LineWidth',2,'Color',[1 0 0])
title('True Output')
ylabel('y(t)')
ylim([-5,5])

subplot(4,1,2)
plot(t,tqEst,'LineWidth',2,'Color','g')
title('TQ model Output')
ylabel('$$\hat{y}$$(t)','Interpreter','Latex');

ylim([-5,5])
subplot(4,1,3)
R=Output-tqEst;
plot(t,R,'LineWidth',2,'Color',[0 0 1])
VAF=100*(1-var(R)/var(Output));

title(['Residuals. %VAF=' num2str(VAF)]);
ylabel('y(t)-$$\hat{y}$$(t)','Interpreter','Latex');

ylabel('$$\hat{R}$$(t)','Interpreter','Latex');
ylim([-5,5])
subplot(4,1,4)
plot(t,K,'LineWidth',2,'Color','r')
h=line([0 10],[kEst kEst]); set(h,'linewidth',2,'color','g');
title('Time-Varying Elasticity');
ylim([5,35]);
ylabel('K(t)');
xlabel('Time(s)');
legend('True','Estimate');
streamer('Tine Invariant Identifcation');


%% Quasi-stationary  IdentificationIdentification 
% We can use matlab's fitlm (fit linear model) function to find the model 
% parameters, other options are avaliable, such as '\' operator. They should provide 
% similar results. 
%estimating a TI system using TV data
snr=0;
iStart=1;
TQ_model={};
tqEst=[];
kEst=[];
segLen=100;
nSeg=N/segLen;
N=1000;

% 
Noise = randn(N,1);
    power_noise = sum(Noise.^2);
    power_signal = sum((Output).^2);   
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);
    
OutputWithNoise=Output+Total_noise; 
    
for i=1:nSeg,
    iEnd=iStart+segLen-1;
 curModel = fitlm([Input(iStart:iEnd) d_Input(iStart:iEnd) dd_Input(iStart:iEnd)],OutputWithNoise(iStart:iEnd),'y ~ x1+x2+x3-1');
 tqEstCur=curModel.Fitted;
 kEstCur=zeros(segLen,1)+table2array(curModel.Coefficients(1,1));
 tqEst=cat(1,tqEst,tqEstCur);
 kEst=cat(1,kEst,kEstCur); 
iStart=iStart+segLen;
end

%estimate the models parameters using matlab's fitlm function 
%passing the command 'y ~ x1+x2+x3-1' tells the function that they are three
%inputs and that no constant value should be used

%plots of input and ouput
figure(1);
subplot(4,1,1)
plot(t,Output,'LineWidth',2,'Color',[1 0 0])
title('True Output')
ylabel('y(t)')
ylim([-5,5])
subplot(4,1,2)
plot(t,tqEst,'LineWidth',2,'Color','g')
title('TQ model Output')
xlabel('Time(s)')
ylabel('$$\hat{y}$$(t)','Interpreter','Latex')
ylim([-5,5])
subplot(4,1,3)
R=Output-tqEst;
plot(t,R,'LineWidth',2,'Color',[0 0 1])
title('Residuals')
xlabel('Time(s)')
ylabel('$$\hat{R}$$(t)','Interpreter','Latex')
ylim([-5,5])
subplot(4,1,4)
plot(t,K,'LineWidth',2,'Color','r')
h=line(t,kEst); set(h,'linewidth',2,'color','g');
title('Time-Varying Elasticity')
ylim([5,35])
ylabel('K(t)')
xlabel('Time(s)')
legend('True','Estimate');

streamer (['Quasistationary Identification SNR= ' num2str(snr)]);

figMod(3,'linewidth',2);


%% 
% The TI approximation provides the mean of the TV parameters. We need other 
% approach to track the changes in the parameter
% 
% *Ensemble estimation:*

%run your experiment 'hundreds' of times to understand what happens 
%at each point in time

%TV Model parameters 
K = [30*ones(N/4,1); linspace(30,10,N/2)';10*ones(N/4,1)]; % TV Elasticity
B = 0.25; %viscosity 
I = 0.002; %inertia

N_trials = 5;
INPUT = zeros(N,N_trials);
d_INPUT = zeros(N,N_trials);
dd_INPUT = zeros(N,N_trials);
OUTPUT = zeros(N,N_trials);
for i =1 : N_trials
    
    Input = filtfilt (fil_n, fil_d,0.1*randn(N,1));
    %input is low-pass filtered, white-Gaussing signal
    d_Input = gradient(Input,T); %first derivative input
    dd_Input = gradient(d_Input,T); %second derivative input
    
    %store this inforation in matrices 
    INPUT(:,i) = Input;
    d_INPUT(:,i) = d_Input;
    dd_INPUT(:,i) = dd_Input;
    
    %simulated ouput
    Output = K.*Input + B*d_Input + I*dd_Input;
    
    OUTPUT(:,i) = Output;
end
figure(1);
subplot (2,1,1);
waterfall(INPUT');
title('INPUT'); 
subplot(2,1,2);
waterfall(OUTPUT');
title('Output');
streamer('Ensemble of Input and Output Trials'); 

    
%% 
% Now that we have information for the system at each point in time, we 
% can estimate a linear model between the input and the output at each point in 
% time. This can be done by perfoming the parameter estimation at each point in 
% time

TV_params_est = zeros(N,3);
for j =1 : N %note that we are covering all the time up to N       
    Model_time_i =fitlm([INPUT(j,:)' d_INPUT(j,:)' dd_INPUT(j,:)'],OUTPUT(j,:)','y ~ x1+x2+x3-1');
    TV_params_est(j,:) = Model_time_i.Coefficients.Estimate;
end

figure(1);clf;
subplot(3,1,1);

plot(t,K,'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est(:,1),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('K')
ylim([5,35])
title(['Number of Trials = ',num2str(N_trials)])

subplot(3,1,2)
plot(t,B*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est(:,2),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('B')
ylim([0.2 0.3])

subplot(3,1,3)
plot(t,I*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est(:,3),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('I')
xlabel('Time (s)')
ylim([0.0015 0.0025])
streamer('Ensemble Estimates - No Noise'); 
%%
%Generate esnemble response with noise

N_trials = 500; %change the number of trials to visualize how this affects the results
INPUT = zeros(N,N_trials);
d_INPUT = zeros(N,N_trials);
dd_INPUT = zeros(N,N_trials);
OUTPUT = zeros(N,N_trials);
OUTPUT_noise = zeros(N,N_trials);
for i =1 : N_trials
    
    Input = filtfilt (fil_n, fil_d,0.1*randn(N,1));
    %input is low-pass filtered, white-Gaussing signal
    d_Input = gradient(Input,T); %first derivative input
    dd_Input = gradient(d_Input,T); %second derivative input
    
    %store this inforation in matrices 
    INPUT(:,i) = Input;
    d_INPUT(:,i) = d_Input;
    dd_INPUT(:,i) = dd_Input;
    
    %simulated ouput
    Output = K.*Input + B*d_Input + I*dd_Input;
    
    OUTPUT(:,i) = Output;
    
    Noise = randn(N,1);
    power_noise = sum(Noise.^2);
    power_signal = sum((Output).^2);
    snr=5;
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);
    
    OUTPUT_noise(:,i) = Output+Total_noise;
end

%
%results with 5 trials 
N_trials = 5;
TV_params_est_noise = zeros(N,3);
for j =1 : N %note that we are covering all the time up to N       
    Model_time_i =fitlm([INPUT(j,1:N_trials)' d_INPUT(j,1:N_trials)' ...
        dd_INPUT(j,1:N_trials)'],OUTPUT_noise(j,1:N_trials)','y ~ x1+x2+x3-1');
    TV_params_est_noise(j,:) = Model_time_i.Coefficients.Estimate;
end

 clf
subplot(3,1,1);
plot(t,K,'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,1),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('K')
ylim([0 40])
title(['Ensemble Identification with SNR=' num2str(snr) ' and number of Trials = ',num2str(N_trials)])
subplot(3,1,2)
plot(t,B*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,2),'LineWidth',2,'Color','g')
ylabel('B')
ylim([0.1 0.4])
subplot(3,1,3)
plot(t,I*ones(N,1),'LineWidth',2,'Color','r')
hold on
plot(t,TV_params_est_noise(:,3),'LineWidth',2,'Color','c')
ylabel('I')
xlabel('Time (s)')
ylim([0.000 0.004])
figMod(5,'lineWidth',2); 

%%
%results with 50 trials 
N_trials = 50;
TV_params_est_noise = zeros(N,3);
for j =1 : N %note that we are covering all the time up to N       
    Model_time_i =fitlm([INPUT(j,1:N_trials)' d_INPUT(j,1:N_trials)' ...
        dd_INPUT(j,1:N_trials)'],OUTPUT_noise(j,1:N_trials)','y ~ x1+x2+x3-1');
    TV_params_est_noise(j,:) = Model_time_i.Coefficients.Estimate;
end

figure(1); clf;
subplot(3,1,1);
plot(t,K,'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,1),'LineWidth',4,'Color',[0 0 1])
legend('True','Estimated')
ylabel('K')
ylim([0,40])
title(['Number of Trials = ',num2str(N_trials)])
subplot(3,1,2)
plot(t,B*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,2),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('B')
ylim([0.1 0.4])
subplot(3,1,3)
plot(t,I*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,3),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('I')
xlabel('Time (s)')
ylim([0.00 0.004])

%%
%results with 500 trials 
N_trials = 500;
TV_params_est_noise = zeros(N,3);
for j =1 : N %note that we are covering all the time up to N       
    Model_time_i =fitlm([INPUT(j,1:N_trials)' d_INPUT(j,1:N_trials)' ...
        dd_INPUT(j,1:N_trials)'],OUTPUT_noise(j,1:N_trials)','y ~ x1+x2+x3-1');
    TV_params_est_noise(j,:) = Model_time_i.Coefficients.Estimate;
end
figure(1); clf
subplot(3,1,1);
plot(t,K,'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,1),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('K')
ylim([0,40])
title(['Number of Trials = ',num2str(N_trials)])
subplot(3,1,2)
plot(t,B*ones(N,1),'LineWidth',2,'Color',[1 0 0])
hold on
plot(t,TV_params_est_noise(:,2),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('B')
ylim([0.1 0.4])
subplot(3,1,3)
plot(t,I*ones(N,1),'LineWidth',2,'Color',[1 0 0])

hold on
plot(t,TV_params_est_noise(:,3),'LineWidth',2,'Color',[0 0 1])
legend('True','Estimated')
ylabel('I')
xlabel('Time (s)')
ylim([0.00 0.004])
%% 
% Ensemble estimation is straighforward, but it needs a lot of data to produce 
% useful results. Getting all the data can be challenging and sometimes imposible. 
% Thus, other alternatives are needed. 
% 
% *Basis function Expansion:*

%create set of basis functions -- in this case we will use cubic B-splines
%for this example, we will use splines with fixed knots. Knot location can be 
%optimized using different approaches, that won't be discussed here

K = [30*ones(N/4,1); linspace(30,10,N/2)';10*ones(N/4,1)]; % TV Elasticity

%create basis function expansions 
q=linspace(0,1,N);
N_basis  = 10;
centers = linspace(0,1,N_basis);  %equaly spaced B-splines  
number=length(centers);
Basis=zeros(length(q),number);
sd = 0.1;
figure(1); clf
for i=1:number
    Basis(:,i)=(1/(2*pi*sd)^(1/2))*exp(-(0.5/(sd^2))*((q-centers(i)).*(q-centers(i))));
     subplot(2,5,i); 
     plot(Basis(:,i),'LineWidth',2);
     title(['Basis function ' num2str(i)]);
end
streamer('Cubic B-spline BasisFunctions'); 
%%
figure;clf
subplot(3,5,[1:5])
plot(t,K,'LineWidth',2,'Color',[0 0 1])
title('TV parameter')
ylabel('K(t)')
ylim([5 35])
for i=1:10,
subplot(3,5,5+i)
plot(t,Basis(:,i),'LineWidth',2,'Color',[1 0 0])
title(['Basis #', num2str(i)])
end

figMod(1,'linewidth',2); 

%%

figure(1); clf

params = Basis\K;
fit_to_K = Basis*params;
subplot (2,1,1);
stem(params);
xlabel('Basis Function Number');
ylabel('Scale Factor'); 
title('Basis Function Parameters for K'); 

subplot (2,1,2);

plot(t,K,'LineWidth',2,'Color',[0 0 1])
hold on
plot(t,fit_to_K,'LineWidth',2,'Color',[1 0 0])
title('K and Its Basis Function Approximation')
ylabel('K(t)')
xlabel('Time (s)')
legend('True', 'Approximation')
figMod(1,'linewidth',2); 

%% Scale basis functions
figure(1); clf
for i=1:10,
subplot(3,5,5+i)
plot(t,fit_to_K(:,i),'LineWidth',2,'Color',[1 0 0])  % Doesn't seem to work
title(['Basis #', num2str(i)])
end

%% 
% Now we use this idea to estimate the model parameters from data 
% 
% will use a single trial and  noise with various SNR
snr=0; %change the SNR value to see how noise affects estimated parameters
%model simulation 
%TV Model parameters 
K = [30*ones(N/4,1); linspace(30,10,N/2)';10*ones(N/4,1)]; % TV Elasticity
B = 0.25; %viscosity 
I = 0.002; %inertia

N_trials = 1; %change the number of trials to visualize how this affects the results
INPUT = zeros(N,N_trials);
d_INPUT = zeros(N,N_trials);
dd_INPUT = zeros(N,N_trials);
OUTPUT = zeros(N,N_trials);
OUTPUT_noise = zeros(N,N_trials);

snrList=[40 20 0];
for iSNR=1:3,
    snr=snrList(iSNR); 

for i =1 : N_trials
    
    Input = filtfilt (fil_n, fil_d,0.1*randn(N,1));
    %input is low-pass filtered, white-Gaussing signal
    d_Input = gradient(Input,T); %first derivative input
    dd_Input = gradient(d_Input,T); %second derivative input
    
    %store this inforation in matrices 
    INPUT(:,i) = Input;
    d_INPUT(:,i) = d_Input;
    dd_INPUT(:,i) = dd_Input;
    
    %simulated ouput
    Output = K.*Input + B*d_Input + I*dd_Input;
    
    OUTPUT(:,i) = Output;
    
    Noise = randn(N,1);
    power_noise = sum(Noise.^2);
    power_signal = sum((Output).^2);
    
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);
    
    OUTPUT_noise(:,i) = Output+Total_noise;
end


%new inputs
X_new = zeros(N,N_basis);
for i = 1 : N_basis
    X_new(:,i)= Input.*Basis(:,i);
end

model_des= 'y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12-1';

TInew_Model =fitlm([X_new d_INPUT dd_INPUT],...
    OUTPUT_noise,model_des);
TInew_params_est = TInew_Model.Coefficients.Estimate;

Est_K = Basis*TInew_params_est(1:N_basis);
Est_B = TInew_params_est(N_basis+1);
Est_I = TInew_params_est(N_basis+2);

figure (iSNR);clf
subplot(2,1,1)
plot(t,OUTPUT,'LineWidth',2,'Color',[0 0 1])
hold on
plot(t,TInew_Model.Fitted,'LineWidth',1,'Color',[1 0 0])
V=vaf(OUTPUT,TInew_Model.Fitted)
ylabel('y(t)')
title(['Noise Free Output and  Estimated Model Predictions (SNR =' , num2str(snr), 'dB) VAF=' num2str(V)])
legend('Noise free measurements','Estimated Model prediction')
ylim([-5 5])
subplot(2,1,2)
plot(t,K,'LineWidth',2,'Color',[0 0 1])
hold on
plot(t,Est_K,'LineWidth',2,'Color',[1 0 0])
title('Elastic parameter')
ylabel('K(t)')
xlabel('Time (s)')
legend('True', 'Estimation')
end

%% 
% We can combine Ensemble learning + basis function expansion to produce 
% a method that is more robust that each method individually
% 
% 

%we will use a single trial and  noise with various SNR
T = 1/100; %Sampling time
duration= 10; % seconds
t = 0:T:duration-T; %Time vector
N= length(t);
%model simulation 
%TV Model parameters 
K = [30*ones(N/4,1); linspace(30,10,N/2)';10*ones(N/4,1)]; % TV Elasticity
B = 0.25; %viscosity 
I = 0.002; %inertia

N_trials = 10; %change the number of trials to visualize how this affects the results
INPUT = zeros(N,N_trials);
d_INPUT = zeros(N,N_trials);
dd_INPUT = zeros(N,N_trials);
OUTPUT = zeros(N,N_trials);
OUTPUT_noise = zeros(N,N_trials);

snr=0; %change the SNR value to see how noise affects estimated parameters

for i =1 : N_trials
    
    Input = filtfilt (fil_n, fil_d,0.1*randn(N,1));
    %input is low-pass filtered, white-Gaussing signal
    d_Input = gradient(Input,T); %first derivative input
    dd_Input = gradient(d_Input,T); %second derivative input
    
    %store this inforation in matrices 
    INPUT(:,i) = Input;
    d_INPUT(:,i) = d_Input;
    dd_INPUT(:,i) = dd_Input;
    
    %simulated ouput
    Output = K.*Input + B*d_Input + I*dd_Input;
    
    OUTPUT(:,i) = Output;
    
    Noise = randn(N,1);
    power_noise = sum(Noise.^2);
    power_signal = sum((Output).^2);
    
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);
    
    OUTPUT_noise(:,i) = Output+Total_noise;
end


%new inputs
X_new = zeros(N,N_basis);
XX_new = [];
for k = 1 : N_trials
    temp = INPUT(:,k);
    for i = 1 : N_basis
        X_new(:,i)= temp.*Basis(:,i);
    end
    
    XX_new = [XX_new;X_new];
end

%params_TI_ensemle = pinv([XX_new reshape(d_INPUT,N*N_trials,1)...
%    reshape(dd_INPUT,N*N_trials,1)])*reshape(OUTPUT_noise,N*N_trials,1);
TInew_Model_ensemble =fitlm([XX_new reshape(d_INPUT,N*N_trials,1)...
    reshape(dd_INPUT,N*N_trials,1)],...
    reshape(OUTPUT_noise,N*N_trials,1),model_des);
TInew_params_est_ensemble = TInew_Model_ensemble.Coefficients.Estimate;

Est_K_ensemble = Basis*TInew_params_est_ensemble (1:N_basis);
Est_B_ensemble = TInew_params_est_ensemble (N_basis+1);
Est_I_ensemble = TInew_params_est_ensemble (N_basis+2);

figure (1);clf
subplot(2,1,1)
plot(t,OUTPUT(:,1),'LineWidth',2,'Color',[0 0 1])
hold on
plot(t,TInew_Model_ensemble.Fitted(1:N),'LineWidth',1,'Color',[1 0 0])
ylabel('y(t)')
V=vaf(OUTPUT(:,1),TInew_Model_ensemble.Fitted(1:N))

title(['Nose free output vs. Estimated Model Predictions (SNR =' , num2str(snr), 'dB) NTrials=' num2str(N_trials) ' VAF='  num2str(V)])
legend('Noise free measurement','Estimated Model prediction')
ylim([-5 5])
subplot(2,1,2)
plot(t,K,'LineWidth',2,'Color',[0 0 1])
hold on
plot(t,Est_K_ensemble,'LineWidth',2,'Color',[1 0 0])
title('TV parameter')
ylabel('K(t)')
xlabel('Time (s)')
legend('True', 'Estimation')
figMod(1,'lineWidth',2); 

%% 
% 
% 
%
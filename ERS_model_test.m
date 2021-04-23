%ERS Model Test (MODEL VALIDATION)

%Validate the estimated model by using it to predict the movements
%generated in response to an input realization which is different than that
%used to identify the model

%INSTRUCTIONS: RUN ERS_model FIRST and THEN RUN THIS

%% Set initial Parameters
noise_snr = [];
set_output_noise_power = 0;
output_noise_power = [];
figNum = 50;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%Desired Displacement Signal Type (PRBS or "Physiological")
PRBS_movement = true;
%physiological_movement = true;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = true;
N = PRBS_movement_time/10;      %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                      %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 10;            %PRBS Amplitude

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds) 
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%Set which model you want to use
if compare_two_models == true
    NHK = NHK2;
end

%Set number of Validation Trials (To get a mean and std VAF)
num_trials = 1;
validation_accuracy = [];

%% Generate Desired Displacement Signals for Validation

for trial = 1:num_trials
    
    if PRBS_movement == true
        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;

        A = 0;                                      %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %First interval is at max PRBS Amplitude
                else
                    R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
                end
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude;                      %Otherwise set as Constant Amplitude
        end

        Range = [0,0.001]; %Specify that the single-channel PRBS value switches between -2 and 2

        %Specify the clock period of the signal as 1 sample. 
        %That is, the signal value can change at each time step. 
        %For PRBS signals, the clock period is specified in Band = [0 B], 
        %where B is the inverse of the required clock period
        %(Must be less than 1)
        Band = [0 0.01];

        %Generate a nonperiodic PRBS of length 100 samples.
        u = idinput(time*1000+1,'prbs',Band,Range);

        %Create an iddata object from the generated signal. 
        %For this example, specify the sample time as 1 second.
        u = iddata([],u,0.001);

        U = (u.InputData)';
        desired_displacement = A.*U;

        %Power Pectrum of Desired Displacement
        desired_displacement_zero = desired_displacement - mean(desired_displacement);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

    else

        % "Physiological" Movements
        t_total = 0:0.001:physiological_movement_time;
        time = physiological_movement_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 1.8;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);

        AR = makedist('Uniform','lower',0,'upper',physiological_movement_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];

        for j = 1 : nf    
            t  = 0 : 0.001 : t_interval;         % Length of Time Interval
            if j == 1
                Freq = FrequenciesRandom_max;
                A = physiological_movement_max_amplitude;
            else
                Freq = random(FrequenciesRandom,1,1);
                A = random(AmplitudesRandom,1,1);
            end

            if chance_of_zero == true
                nums = randi([0 1], 1, 1);
            else
                nums = 0;
            end

            if nums == 0
                g = 1/Freq;
                D = (1:g:t_interval)';     % pulse delay times
                data = (A*pulstran(t,D,@rectpuls,W))';
                movement_frequency = Freq;
                Pulses_per_interval = t_interval/g;
                data(end) = [];
            else
                data = zeros(20000,1);
                movement_frequency = 0;
                Pulses_per_interval = 0;
            end

            %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
            desired_displacement = [desired_displacement; data];
            Freq_test = [Freq_test movement_frequency];
            Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
        end

        Pulses_per_interval_total = sum(Pulses_per_interval_test);
        Freq_test_average = sum(Freq_test)/length(Freq_test);

        %Power Pectrum of Desired Displacement
        desired_displacement_zero = desired_displacement - mean(desired_displacement);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

        desired_displacement = desired_displacement';

    end

    %% Create Neural Input for ERS Simulation
    %Create Frequency and Amplitude Parameters of the Neural Input based on
    %Desired Displacement Amplitude
    Amplitude = desired_displacement*100;    %mV
    Frequency = desired_displacement*14000;  %Hz

    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model

    %Power Pectrum of Neural Input
    neural_zero = neural - mean(neural);
    [Pxx2,f2] = pwelch(neural_zero,Nfft,[],Nfft,Fs);


    %% Execute ERS Simulation Simulink Model

    %Set Output Noise (set as zero)
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink Model
    out = sim('ERS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %% Get Output Signals from ERS Simulation Model
    
    %EMG
    emg_simulink = out.EMGout;
    emg_pdf = emg_simulink;
    emg_pdf(emg_pdf==0) = []; %Removing the zero values
    
    %Muscle Force
    force_simulink = out.EMG_Model_Force;
    
    %Healthy Displacement Output and Time
    output_displacement_simulink = out.EMG_Model_Displacement;
    t_simulink = out.tout;


    %% Input/Output for Model Validation
    
    Zcur = [emg_simulink,output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

    %% Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    %% Model Validation
    %Test the Model to predict movement from EMG
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

    %Plot Results of Identification
    figure(figNum)
    %figNum = figNum+1;
    [R, V, yp] = nlid_resid(NHK,Zcur);

    validation_accuracy = [validation_accuracy V];
    
end

%Calculate validation mean and Std
validation_accuracy_mean = mean(validation_accuracy)
validation_accuracy_std = std(validation_accuracy)

%% Plot the Validation results of last validation trial

%Plots just the superimposed part with %VAF in the title
figNum = figNum+1;
figure(figNum);
figNum = figNum+1;
pred = double(yp);
plot(t_total,pred);
hold on
plot(t_total, output_displacement_simulink)
ax = gca;
ax.FontSize = 14;
hold off
title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 24)
xlabel('Time (s)', 'Fontsize', 18)
ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 18)
legend('Predicted', 'Observed','Fontsize',14)

%Plot Residuals
figure(figNum)
figNum = figNum+1;
subplot(3,1,1)
plot(R)
ax = gca;
ax.FontSize = 12;
title('Residuals of ERS Hammerstein Model','Fontsize',20)
xlabel('Time (s)','Fontsize',18)
ylabel('Displacement (m)','Fontsize',18)

subplot(3,1,2)
p = pdf(R);
plot(p)
ax = gca;
ax.FontSize = 12;
xlabel('Displacement (m)','Fontsize',18)
ylabel('Density','Fontsize',18)
title('Residual Distribution','Fontsize',20)

R_zero = R - mean(R);
S = spect(R_zero);
S_frequency = 0.0556:0.0556:0.0556*length(S);
% subplot(3,1,3)
% semilogy(S_frequency(:,1:50),S(1:50,:),'LineWidth',1.5);
% ax = gca;
% ax.FontSize = 12;
% title('Power Spectrum of Residuals','Fontsize',20);
% ylabel('PSD (log)','Fontsize',18); 
% xlabel('Frequency (Hz)','Fontsize',18);
% grid on

[PxxR,fR] = pwelch(double(R_zero),[],[],[],Fs);
subplot(3,1,3)
%semilogy(fR(1:110,:),PxxR(1:110,:),'LineWidth',1.5);
semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
ax = gca;
ax.FontSize = 12;
title('(c) Power Spectrum of Residuals','Fontsize',20);
ylabel('PSD (log)','Fontsize',18); 
xlabel('Frequency (Hz)','Fontsize',18);
grid on

figure(figNum)
figNum = figNum+1;
subplot(2,2,[1 2])
plot(t_total,pred);
hold on
plot(t_total, output_displacement_simulink)
ax = gca;
ax.FontSize = 12;
hold off
title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 18)
xlabel('Time (s)', 'Fontsize', 16)
ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 16)
legend('Predicted', 'Observed','Fontsize', 14)

subplot(2,2,3)
plot(p)
ax = gca;
ax.FontSize = 12;
xlabel('Displacement (m)','Fontsize',16)
ylabel('Density','Fontsize',16)
title('(b) Residual Distribution','Fontsize',18)
grid on

subplot(2,2,4)
semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
ax = gca;
ax.FontSize = 12;
title('(c) Power Spectrum of Residuals','Fontsize',18);
ylabel('PSD (log)','Fontsize',16); 
xlabel('Frequency (Hz)','Fontsize',16);
grid on

figure(figNum)
figNum = figNum+1;
subplot(2,2,[1 2])
plot(t_total,pred);
hold on
plot(t_total, output_displacement_simulink)
ax = gca;
ax.FontSize = 12;
hold off
title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 18)
xlabel('Time (s)', 'Fontsize', 16)
ylabel('Displacement, Pos_H(t) (m)', 'Fontsize', 16)
legend('Predicted', 'Observed','Fontsize', 14)

subplot(2,2,3)
plot(p)
ax = gca;
ax.FontSize = 12;
xlabel('Displacement (m)','Fontsize',16)
ylabel('Density','Fontsize',16)
title('(b) Residual Distribution','Fontsize',18)
grid on

subplot(2,2,4)
R_double = double(R_zero);
L = length(R_double)-1;
maxlag = L*0.10;
[res_corr,lags] = xcorr(R_double,maxlag,'normalized');
%[res_corr,lags] = xcov(R,maxlag,'coeff');
lags = lags/Fs;
subplot(2,2,4),plot(lags,res_corr,'LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
title('(c) Auto-Correlation Residuals','Fontsize',18)
xlabel('Time (s)','Fontsize',16)
ylabel('Correlation','Fontsize',16)
grid on

figure(figNum)
figNum = figNum+1;
Rx = double(R);
plot(Rx, output_displacement_simulink)
ax = gca;
ax.FontSize = 14;
title('Residuals vs Displacement Amplitude','Fontsize',24);
ylabel('Displacement Amplitude (m)','Fontsize',18); 
xlabel('Residuals (m)','Fontsize',18);
grid on
    




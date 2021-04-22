% EMG Model Experiment (Accuracy vs Ouput Noise)

%Accuracy vs Noise with Constant Record Length
%Increases the output noise(measurement error) in simulink model

%Identifies one set of models for each signal type (PRBS and Physiological)
%with increasing levels of noise for each model (Identifies 2x
%'noise_level_iters' models)
%Validates each identified model with 'num_trials' number of physiological
%signals
%Plots identification accuracy and monte carlo validation accuracy


%UPDATED TO COMPARE WITH CLEAN SIGNAL

clc
clear all

tStart = tic;

%% Set initial Parameters
noise_level_iters = 18;                             %Number of times the noise power is multiplied
set_output_noise_initial = 1e-15;                   %Initial Output Noise Power
set_output_noise_power = set_output_noise_initial;
noise_multiplier = 5;                               %Noise Power Multiplier
set_seed = 23341;

%For Power Spectrums or FFT
Fs = 1000; 
Nfft = 10000;

%PRBS Stimulus
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10; %mm

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = 18;                                        %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

accuracy_identification = [];
accuracy_identification_all = [];
accuracy_validation = [];
accuracy_validation_all = [];

noise_snr_clean = [];
noise_snr = [];
noise_snr_all_clean = [];
noise_snr_all = [];
output_noise_power = [];
output_noise_power_all = [];

NHK_all = [];
models_all = [];
Zcur_all = [];
emg_all = [];

%%
num_trials = 30;
signals = 2;

%%
for signal = 1:signals
    
    if signal == 1
        PRBS_movement = false;
    else
        PRBS_movement = true;
    end
    
    if PRBS_movement == true
        
        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;
        
        A = 0;                    %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %First interval is at max amplitude
                else
                    R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and 10
                end
                
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude;              %Constant Amplitude
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

        %Plot, and examine the generated signal.
        %plot(u);
        %title('Non-Periodic Signal')

        U = (u.InputData)';
        desired_displacement = A.*U;
        
    else
        
        %(Physiologically Based Movment - Square pulses with random amplitude and random frequency)
        t_total = 0:0.001:physiological_movement_time;
        time = physiological_movement_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 1.8;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);

        AR = makedist('Uniform','lower',0,'upper',physiological_movement_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement=[0];
        Freq_test = [];
        Pulses_per_interval_test = [];

        for j = 1 : nf    
            t  = 0 : 0.001 : t_interval;         % Time Samples

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
        [Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        desired_displacement = desired_displacement';
        
    end
    
    %%
    %Create Frequency and Amplitude Parameters based on Analog Signal
    Amplitude = desired_displacement*100;  %mV
    Frequency = desired_displacement*14000;  %Hz
    
    %%
    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model

    %Power Pectrum of Neural Input
    [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    %% Generate Clean Signals for Evaluation
    %Execute Simulink Model with No Noise
    
    set_output_noise_power_clean = 0;
    set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power_clean')
    %output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('EMG_Model_Simulink',time);
    
    %Get Output Signals from Simulink
    emg_simulink_clean = out.EMGout;
    force_simulink_clean = out.EMG_Model_Force;
    
    output_displacement_simulink_clean = out.EMG_Model_Displacement;
    t_simulink_clean = out.tout;

    Zcur_clean = [emg_simulink_clean,output_displacement_simulink_clean];
    Zcur_clean = nldat(Zcur_clean,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
    
    %% Add Noise to the Model
    for noise_level = 1:noise_level_iters
        
        %Execute Simulink Model with Set Output Noise

        %Set Output Noise
        set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink Model
        out = sim('EMG_Model_Simulink',time);

        set_output_noise_power = noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %%
        %Get Output Signals from Simulink
        emg_simulink = out.EMGout;

        emg_pdf = emg_simulink;
        emg_pdf(emg_pdf==0) = []; %Removing the zero values

        %Power Pectrum of EMG
        [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        force_simulink = out.EMG_Model_Force;

        output_displacement_simulink = out.EMG_Model_Displacement;
        t_simulink = out.tout;

        Zcur = [emg_simulink,output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

        %%
        %Calculate Signal to Noise Ratio
        output_noise_simulink = out.Output_Noise;
        
        %SNR between clean signal and noise
        signal_to_noise_clean = snr(output_displacement_simulink_clean, output_noise_simulink);
        noise_snr_clean = [noise_snr_clean signal_to_noise_clean];
        
        %SNR between clean signal + noise and noise
        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];
        
        %% Model Identification
        %Hammerstein System Identification (EMG-Movement Model)
        %HK Method
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        NHK=nlbl;
        set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

        % Set umber of lags in IRF
        I=NHK{1,2};
        set(I,'nLags',2400); %(accuracy increases if nlags is increased)
        NHK{1,2}=I;
        
        %Identify the model with the unclean output
        NHK=nlident(NHK,Zcur);
        
        % Evaluate by comparing with the clean Output
        figure(signal);
        [R, V, yp] = nlid_resid(NHK,Zcur_clean);

        accuracy_identification = [accuracy_identification V];
        
        NHK_all = [NHK_all NHK];
        NHK = [];
        Zcur_all = [Zcur_all Zcur];
        emg_all = [emg_all emg_simulink];
        
    end
    
    set_output_noise_power = set_output_noise_initial;
        
    accuracy_identification_all(signal,:) = accuracy_identification;
    noise_snr_all_clean(signal,:) = noise_snr_clean;
    noise_snr_all(signal,:) = noise_snr;
    output_noise_power_all(signal,:) = output_noise_power;
    
    models_all = [models_all;NHK_all];
    
    accuracy_identification = [];
    noise_snr_clean = [];
    noise_snr = [];
    output_noise_power = [];
    NHK_all = [];
    
end

%% Model Validation with Physiological Signals
%Set initial Parameters
set_output_noise_power_validation = 0;

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = 18;                  %number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%%
for trial = 1:num_trials
    
    %(Physiologically Based Movment - Square pulses with random amplitude and random frequency)
    t_total = 0:0.001:physiological_movement_time;             %Total time
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
        t  = 0 : 0.001 : t_interval;         % Time Samples

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
    [Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

    desired_displacement = desired_displacement';
    
    %%
    %Create Frequency and Amplitude Parameters based on Analog Signal
    Amplitude = desired_displacement*100;  %mV
    Frequency = desired_displacement*14000;  %Hz

    %%
    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model

    %Power Pectrum of Neural Input
    [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    %%
    %Execute Simulink Model with Set Output Noise

    %Set Output Noise
    set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power_validation')
    %output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink Model
    out = sim('EMG_Model_Simulink',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;
    
    %%
    %Get Output Signals from Simulink
    emg_simulink = out.EMGout;
    
    emg_pdf = emg_simulink;
    emg_pdf(emg_pdf==0) = []; %Removing the zero values
    
    %Power Pectrum of EMG
    [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    force_simulink = out.EMG_Model_Force;
    
    output_displacement_simulink = out.EMG_Model_Displacement;
    t_simulink = out.tout;

    Zcur = [emg_simulink,output_displacement_simulink];

    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
    
    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    %noise_snr = [noise_snr signal_to_noise]
    
    %%
    %Model Validation
    
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    for signal = 1:signals
        
        for noise_level = 1:noise_level_iters
            
            figure(signal)
            [R, V, yp] = nlid_resid(models_all(signal, noise_level),Zcur);
            
            accuracy_validation(trial,noise_level,signal) = V;
            
        end
        
    end
    
end

%% Plot Accuracy (Identification & Validation) vs Noise
figNum = 8000;

accuracy_identification_all = max(0, accuracy_identification_all);

figure(figNum)
figNum = figNum+1;
plot(noise_snr_all(2,:),accuracy_identification_all(2,:),'LineWidth',3);
hold on
plot(noise_snr_all(1,:),accuracy_identification_all(1,:),'LineWidth',3);
hold off
ax = gca;
ax.FontSize = 15;
title('Identification Accuracy vs Output Noise','Fontsize',24);
legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

figure(figNum)
figNum = figNum+1;
plot(noise_snr_all_clean(2,:),accuracy_identification_all(2,:),'LineWidth',3);
hold on
plot(noise_snr_all_clean(1,:),accuracy_identification_all(1,:),'LineWidth',3);
hold off
ax = gca;
ax.FontSize = 15;
title('Identification Accuracy vs Output Noise','Fontsize',24);
legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

accuracy_validation = max(0, accuracy_validation);
accuracy_validation_mean = mean(accuracy_validation);
accuracy_validation_std = std(accuracy_validation);

figure(figNum)
figNum = figNum+1;
hold on
plot(noise_snr_all(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(noise_snr_all(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3);

plot(noise_snr_all(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all(2,:) fliplr(noise_snr_all(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(noise_snr_all(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all(1,:) fliplr(noise_snr_all(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('Validation Accuracy vs Output Noise','Fontsize',24);
legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

figure(figNum)
figNum = figNum+1;
hold on
plot(noise_snr_all_clean(2,[1:6 8:end]),accuracy_validation_mean(:,[1:6 8:end],2),'LineWidth',3);
plot(noise_snr_all_clean(1,[1:10 12:end]),accuracy_validation_mean(:,[1:10 12:end],1),'LineWidth',3);

plot(noise_snr_all_clean(2,[1:6 8:end]),min(accuracy_validation_mean(:,[1:6 8:end],2)+accuracy_validation_std(:,[1:6 8:end],2),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(2,[1:6 8:end]),max(accuracy_validation_mean(:,[1:6 8:end],2)-accuracy_validation_std(:,[1:6 8:end],2),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(2,[1:6 8:end]) fliplr(noise_snr_all_clean(2,[1:6 8:end]))], [min(accuracy_validation_mean(:,[1:6 8:end],2)+accuracy_validation_std(:,[1:6 8:end],2),100) fliplr(max(accuracy_validation_mean(:,[1:6 8:end],2)-accuracy_validation_std(:,[1:6 8:end],2),0))], 'b','FaceAlpha',0.2)

plot(noise_snr_all_clean(1,[1:10 12:end]),min(accuracy_validation_mean(:,[1:10 12:end],1)+accuracy_validation_std(:,[1:10 12:end],1),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(1,[1:10 12:end]),max(accuracy_validation_mean(:,[1:10 12:end],1)-accuracy_validation_std(:,[1:10 12:end],1),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(1,[1:10 12:end]) fliplr(noise_snr_all_clean(1,[1:10 12:end]))], [min(accuracy_validation_mean(:,[1:10 12:end],1)+accuracy_validation_std(:,[1:10 12:end],1),100) fliplr(max(accuracy_validation_mean(:,[1:10 12:end],1)-accuracy_validation_std(:,[1:10 12:end],1),0))], 'y','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('Validation Accuracy vs Output Noise','Fontsize',24);
legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on


tEnd = toc(tStart)/60


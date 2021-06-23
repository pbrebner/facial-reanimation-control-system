%% Effects of Noise on the ERS, SRS, and FRCS

%Noise is added to the output of the models to quantify the effects on
%accuracy

%Runs completely through the ERS noise, followed by the SRS noise, and then
%finally the models identified in those forst two sections are used to
%build the FRCS

%Plots the noise vs accuracy at the end of each section

%% ERS Model Test for Accuracy vs Ouput Noise (Measurement Error)

%Accuracy vs Noise with Constant Record Length
%Iteratively increases the output noise(measurement error) in simulink
%model, identifies models with this noisy signal, validates the models,
%and plots the results

%Identifies one set of models for each signal type (PRBS and Physiological)
%with increasing levels of noise for each model (Identifies 2x 'noise_level_iters' models)
%Validates each identified model with 'num_trials' number of physiological
%signals to get a VAF mean and std
%Plots identification accuracy vs noise and validation accuracy vs noise

%UPDATED TO COMPARE WITH CLEAN SIGNAL

%When running the script, you need to provide the following input:
% 1. Number of Validation Trials?
%       Number of validation trials used to calculate the VAF mean and std

clc
clear all

%% User Input Prompts

prompt1 = 'Number of Validation Trials? 1-50 [30]: ';
str1 = input(prompt1);
if str1<1 | str1>50
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 30;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
noise_level_iters = 18;                             %Number of times the noise power is multiplied
set_output_noise_initial = 1e-15;                   %Initial Output Noise Power
set_output_noise_power = set_output_noise_initial;  %Noise power passed into ERS simulation
noise_multiplier = 5;                               %Noise Power Multiplier
set_seed = 23341;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = false;        %PRBS can either be constant amplitude or variable amplitude
N = PRBS_movement_time/10;         %Number of times the amplitude randomly changes (For Variable Amplitude Only)     
M = 10000;                         %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 10;               %PRBS Amplitude (mm)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;    %Physiological Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds)
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

ERS_accuracy_identification = [];
ERS_accuracy_identification_all = [];
ERS_accuracy_validation = [];
ERS_accuracy_validation_all = [];

ERS_noise_snr_clean = [];
ERS_noise_snr = [];
ERS_noise_snr_all_clean = [];
ERS_noise_snr_all = [];
ERS_output_noise_power = [];
ERS_output_noise_power_all = [];

%Initialize Models
NHK_all = [];
ERS_models_all = [];
ERS_Zcur_ident = [];
ERS_Zcur_clean_all = [];
ERS_emg_all = [];

%% Num trials and Signal Types
num_trials = str1;
%signals = 2;
signals = 1; %Only the Phys

%% Generate the Desired Displacement for Model Identification

% Generates a Physiological signal first, identifies models with increasing
% noise levels, evaluates the predicted output of these models with a clean output,
% stores these models, and then repeats for PRBS signal
for signal = 1:signals
    
    if signal == 1
        PRBS_movement = false;
    else
        PRBS_movement = true;
    end
    
    if PRBS_movement == true
        
        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;
        
        A = 0;                                      %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %First interval is at max PRBS amplitude
                else
                    R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS amplitude
                end
                
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude;                     %Else set as Constant PRBS Amplitude
        end

        Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

        %Specify the clock period of the signal as 1 sample. 
        %That is, the signal value can change at each time step. 
        %For PRBS signals, the clock period is specified in Band = [0 B], 
        %where B is the inverse of the required clock period
        %(Must be less than 1)
        Band = [0 0.01];

        %Generate a nonperiodic PRBS of length time.
        u = idinput(time*1000+1,'prbs',Band,Range);

        %Create an iddata object from the generated signal. 
        %For this example, specify the sample time as 0.001 second.
        u = iddata([],u,0.001);

        U = (u.InputData)';
        desired_displacement = A.*U;
        
    else
        
        %"Physiological" Movement
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
            t  = 0 : 0.001 : t_interval;         % Time Intervals

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
                D = (1:g:t_interval)';          % pulse delay times
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

        desired_displacement = desired_displacement';
        
    end
    
    %% Create Neural Input for ERS Simulation
    
    %Create Frequency and Amplitude Parameters of Neural Input based on
    %Desired Displacement Amplitude
    Amplitude = desired_displacement*100;    %mV
    Frequency = desired_displacement*14000;  %Hz
    
    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model
    
    %% Generate Clean Signals for Evaluation (To compare with Noisy Signals)
    
    %Execute Simulink Model with No Noise
    set_output_noise_power_clean = 0;
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power_clean')

    %Run Simulink;
    out = sim('ERS_simulation',time);
    
    %Get Clean Outputs from Simulink (EMG, Force, Healthy Displacement)
    emg_simulink_clean = out.ERS_Simulation_EMG;
    force_simulink_clean = out.ERS_Simulation_Force;
    output_displacement_simulink_clean = out.ERS_Simulation_Displacement;
    t_simulink_clean = out.tout;
    
    %Clean Input/Output
    Zcur_clean = [emg_simulink_clean,output_displacement_simulink_clean];
    Zcur_clean = nldat(Zcur_clean,'domainIncr',0.001,'comment','Output EMG, Clean Output Displacement','chanNames', {'EMG (V)' 'Clean Displacement (m)'});
    
    %% Iteratively Add Noise to the Model
    
    for noise_level = 1:noise_level_iters
        
        %Execute Simulink Model with Set Output Noise
        set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
        ERS_output_noise_power = [ERS_output_noise_power set_output_noise_power];

        %Run Simulink Model
        out = sim('ERS_simulation',time);

        set_output_noise_power = noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %% Get Noisy Outputs from ERS simulation (Simulink)
        
        %EMG, Muscle Force, and Healthy Displacement Output
        emg_simulink = out.ERS_Simulation_EMG;
        force_simulink = out.ERS_Simulation_Force;
        output_displacement_simulink = out.ERS_Simulation_Displacement;
        t_simulink = out.tout;
        
        %% Noisy Input/Output for Model Identification
        
        Zcur = [emg_simulink,output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

        %% Calculate Signal to Noise Ratio (SNR)
        
        %Get output noise from ERS simulation (Simulink)
        output_noise_simulink = out.Output_Noise;
        
        %SNR between clean signal and noise
        signal_to_noise_clean = snr(output_displacement_simulink_clean, output_noise_simulink);
        ERS_noise_snr_clean = [ERS_noise_snr_clean signal_to_noise_clean];
        
        %SNR between clean signal + noise and noise
        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        ERS_noise_snr = [ERS_noise_snr signal_to_noise];
        
        %% Model Identification
        
        %Hammerstein System Identification (HK Method)
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        NHK=nlbl;
        set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

        % Set Number of lags in IRF
        I=NHK{1,2};
        set(I,'nLags',2400);
        NHK{1,2}=I;
        
        %Identify the model with the Noisy Output
        NHK=nlident(NHK,Zcur);
        
        %Evaluate by comparing with the Clean Output
        figure(signal);
        [R, V, yp] = nlid_resid(NHK,Zcur_clean);

        ERS_accuracy_identification = [ERS_accuracy_identification V];
        
        %Store the Models Identified for each noise level
        NHK_all = [NHK_all NHK];
        NHK = [];
        ERS_Zcur_ident = [ERS_Zcur_ident Zcur];
        ERS_emg_all = [ERS_emg_all emg_simulink];
        
    end
    
    %Reset the set noise power to the inital noise power (For next Signal
    %Type)
    set_output_noise_power = set_output_noise_initial;
    
    %Store the accuracy, SNR, and identified models for this signal type
    %over all noise levels
    ERS_accuracy_identification_all(signal,:) = ERS_accuracy_identification;
    ERS_noise_snr_all_clean(signal,:) = ERS_noise_snr_clean;
    ERS_noise_snr_all(signal,:) = ERS_noise_snr;
    ERS_output_noise_power_all(signal,:) = ERS_output_noise_power;
    ERS_models_all = [ERS_models_all;NHK_all];
    ERS_Zcur_clean_all = [ERS_Zcur_clean_all Zcur_clean];
    
    %Reset for next signal type
    ERS_accuracy_identification = [];
    ERS_noise_snr_clean = [];
    ERS_noise_snr = [];
    ERS_output_noise_power = [];
    NHK_all = [];
    
end

%% Set Initial Parameters for Model Validation with Physiological Signals

%Noise Parameters for Validation
set_output_noise_power_validation = 0;

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds)
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%% Generate "Physiological" Desired Displacement Signals for Model Validation

for trial = 1:num_trials
    
    % "Physiological" Movement
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
        t  = 0 : 0.001 : t_interval;         % Time Intervals

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

    desired_displacement = desired_displacement';
    
    %% Create Neural Input for ERS Simulation
    %Create Frequency and Amplitude Parameters of Neural Input based on
    %Desired DIsplacement Amplitude
    Amplitude = desired_displacement*100;  %mV
    Frequency = desired_displacement*14000;  %Hz

    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model
    
    %% Execute ERS simulation (Simulink)

    %Set Output Noise (validation noise of zero)
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power_validation')

    %Run Simulink Model
    out = sim('ERS_simulation',time);
    
    %% Get Output Signals from Simulink
    
    %EMG, Muscle Force, and Healthy Displacement Output
    emg_simulink = out.ERS_Simulation_EMG;
    force_simulink = out.ERS_Simulation_Force;
    output_displacement_simulink = out.ERS_Simulation_Displacement;
    t_simulink = out.tout;

    %% Input/Output for Model Validation
    
    Zcur = [emg_simulink,output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
    
    %% Model Validation
    
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    %Run through each identified model (for different noise levels) of each signal
    %type
    for signal = 1:signals
        
        for noise_level = 1:noise_level_iters
            
            figure(signal)
            [R, V, yp] = nlid_resid(ERS_models_all(signal, noise_level),Zcur);
            
            ERS_accuracy_validation(trial,noise_level,signal) = V;
            
        end
        
    end
    
end

%% Plot Accuracy (Identification & Validation) vs Output Noise
figNum = 400;

ERS_accuracy_identification_all = max(0, ERS_accuracy_identification_all);

% figure(figNum)
% figNum = figNum+1;
% plot(noise_snr_all(2,:),accuracy_identification_all(2,:),'LineWidth',3);
% hold on
% plot(noise_snr_all(1,:),accuracy_identification_all(1,:),'LineWidth',3);
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

figure(figNum)
figNum = figNum+1;
% plot(ERS_noise_snr_all_clean(2,:),ERS_accuracy_identification_all(2,:),'LineWidth',3);
hold on
plot(ERS_noise_snr_all_clean(1,:),ERS_accuracy_identification_all(1,:),'LineWidth',3);
hold off
ax = gca;
ax.FontSize = 15;
title('ERS Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

ERS_accuracy_validation = max(0, ERS_accuracy_validation);
ERS_accuracy_validation_mean = mean(ERS_accuracy_validation);
ERS_accuracy_validation_std = std(ERS_accuracy_validation);

% figure(figNum)
% figNum = figNum+1;
% hold on
% plot(noise_snr_all(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
% plot(noise_snr_all(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3);
% 
% plot(noise_snr_all(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
% plot(noise_snr_all(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
% patch([noise_snr_all(2,:) fliplr(noise_snr_all(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)
% 
% plot(noise_snr_all(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
% plot(noise_snr_all(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
% patch([noise_snr_all(1,:) fliplr(noise_snr_all(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Validation Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

figure(figNum)
figNum = figNum+1;
hold on
% plot(ERS_noise_snr_all_clean(2,:),ERS_accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(ERS_noise_snr_all_clean(1,:),ERS_accuracy_validation_mean(:,:,1),'LineWidth',3);

% plot(ERS_noise_snr_all_clean(2,:),min(ERS_accuracy_validation_mean(:,:,2)+ERS_accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
% plot(ERS_noise_snr_all_clean(2,:),max(ERS_accuracy_validation_mean(:,:,2)-ERS_accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
% patch([ERS_noise_snr_all_clean(2,:) fliplr(ERS_noise_snr_all_clean(2,:))], [min(ERS_accuracy_validation_mean(:,:,2)+ERS_accuracy_validation_std(:,:,2),100) fliplr(max(ERS_accuracy_validation_mean(:,:,2)-ERS_accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(ERS_noise_snr_all_clean(1,:),min(ERS_accuracy_validation_mean(:,:,1)+ERS_accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(ERS_noise_snr_all_clean(1,:),max(ERS_accuracy_validation_mean(:,:,1)-ERS_accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([ERS_noise_snr_all_clean(1,:) fliplr(ERS_noise_snr_all_clean(1,:))], [min(ERS_accuracy_validation_mean(:,:,1)+ERS_accuracy_validation_std(:,:,1),100) fliplr(max(ERS_accuracy_validation_mean(:,:,1)-ERS_accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('ERS Validation Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on


%% SRS Model Test for Accuracy vs Ouput Noise (Measurement Error)

%Accuracy vs Noise with Constant Record Length
%Iteratively increases the output noise(measurement error) in simulink
%model, identifies models with this noisy signal, validates the models,
%and plots the results

%Identifies one set of models for each signal type (PRBS and Physiological)
%with increasing levels of noise for each model (Identifies 2x 'noise_level_iters' models)
%Validates each identified model with 'num_trials' number of physiological
%signals to get a VAF mean and std
%Plots identification accuracy vs noise and validation accuracy vs noise

%UPDATED TO COMPARE WITH CLEAN SIGNAL

%When running the script, you need to provide the following input:
% 1. Type of Model Structure? LNL/Hammerstein/Wiener/IRF
%       The Model structure used for the Identified SRS (Default is Wiener)
% 2. Number of Validation Trials?
%       Number of validation trials used to calculate the VAF mean and std


%% User Input Prompts

prompt1 = 'Type of Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Wiener]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'LNL') & ~strcmp(str1,'Hamm') & ~strcmp(str1,'Wiener') & ~strcmp(str1,'IRF') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Wiener';
end

prompt2 = 'Number of Validation Trials? 1-50 [30]: ';
str2 = input(prompt2);
if str2<1 | str2>50
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 30;
end

%% Set initial Parameters

%Noise Parameters
noise_level_iters = 18;                             %Number of times the noise power is multiplied
set_output_noise_initial = 1e-15;                   %Initial Output Noise Power
set_output_noise_power = set_output_noise_initial;  %Noise power passed into ERS simulation
noise_multiplier = 5;                               %Noise Power Multiplier
set_seed = 23341;

%PRBS Signal Parameters
PRBS_stimulus_time = 180;
variable_amplitude = true;      %PRBS can either be constant amplitude or variable amplitude
N = PRBS_stimulus_time/10;      %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                      %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 20;            %PRBS Amplitude (mm)

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;    %Physiological Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                                       %Width of signal pulse (seconds)
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%Model Type
LNL_model = false;
Hammerstein_model = false;
Weiner_model = false;
Linear_IRF_model = false;

if strcmp(str1,'LNL')
    LNL_model = true;
elseif strcmp(str1,'Hamm')
    Hammerstein_model = true;
elseif strcmp(str1,'Wiener')
    Weiner_model = true;
elseif strcmp(str1,'IRF')
    Linear_IRF_model = true;
end

SRS_accuracy_identification = [];
SRS_accuracy_identification_all = [];
SRS_accuracy_validation = [];

SRS_noise_snr_clean = [];
SRS_noise_snr = [];
SRS_noise_snr_all_clean = [];
SRS_noise_snr_all = [];
SRS_output_noise_power = [];
SRS_output_noise_power_all = [];

%Initialize Models
LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];

SRS_models_all = [];

SRS_Zcur_ident = [];
SRS_Zcur_clean_all = [];

%% Num trials and Signal Types
num_trials = str2;
%signals = 2;
signals = 1; %Only the Phys

%% Generate the Desired Displacement and Amplitude Modulation for Model Identification

% Generates a Physiological signal first, identifies models with increasing
% noise levels, evaluates the predicted output of these models with a clean output,
% stores these models, and then repeats for PRBS signal
for signal = 1:signals
    
    if signal == 1
        PRBS_stimulus = false;
    else
        PRBS_stimulus = true;
    end
    
    if PRBS_stimulus == true

        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;

        A = [0];                                    %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %First interval is at max amplitude
                else
                    R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS amplitude
                end
                
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude;                     %Else set as Constant Amplitude
        end

        Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

        %Specify the clock period of the signal as 1 sample. 
        %That is, the signal value can change at each time step. 
        %For PRBS signals, the clock period is specified in Band = [0 B], 
        %where B is the inverse of the required clock period
        %(Must be less than 1)
        Band = [0 0.01];

        %Generate a nonperiodic PRBS of length time.
        u = idinput(time*1000+1,'prbs',Band,Range);

        %Create an iddata object from the generated signal. 
        %For this example, specify the sample time as 0.001 second.
        u = iddata([],u,0.001);

        U = (u.InputData)';
        desired_displacement = A.*U;

        stim_frequency = 50;
        stim_amplitude = desired_displacement*170;

        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    else

        t_total = 0:0.001:physiological_stimulus_time;
        time = physiological_stimulus_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 2.1;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);

        AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];

        for j = 1 : nf    
            t  = 0 : 0.001 : t_interval;         % Time Intervals

            if j == 1
                Freq = FrequenciesRandom_max;
                A = physiological_stimulus_max_amplitude;
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
                stim_frequency = Freq;
                Pulses_per_interval = t_interval/g;
                data(end) = [];
            else
                data = zeros(20000,1);
                stim_frequency = 0;
                Pulses_per_interval = 0;
            end

            %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
            desired_displacement = [desired_displacement; data];
            Freq_test = [Freq_test stim_frequency];
            Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
        end

        Pulses_per_interval_total = sum(Pulses_per_interval_test);
        Freq_test_average = sum(Freq_test)/length(Freq_test);

        desired_displacement = desired_displacement';

        stim_amplitude = desired_displacement*170;  %mV
        stim_frequency = 50;

        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    end
    
    %% Generate Clean Signals for Evaluation
    
    %Execute Simulink Model with No Noise
    set_output_noise_power_clean = 0;
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power_clean')

    %Run Simulink;
    out = sim('SRS_simulation',time);
    
    %Get Clean Output Signals from Simulink
    %Muscle Force
    force_simulink_clean = out.SRS_Simulation_Force;
    
    %Clean Output Paralyzed Displacement
    output_displacement_simulink_clean = out.SRS_Simulation_Displacement;
    t_simulink_clean = out.tout;

    %Clean Input/Output
    Zcur_clean = [amplitude_modulation',output_displacement_simulink_clean];
    Zcur_clean = nldat(Zcur_clean,'domainIncr',0.001,'comment','Input Amplitude Modulation, Clean Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Clean Displacement (m)'});
    
    %% Add Noise to the Model
    
    for noise_level = 1:noise_level_iters
        
        %Execute Simulink Model with Set Output Noise
        set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
        SRS_output_noise_power = [SRS_output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('SRS_simulation',time);

        set_output_noise_power = noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
    
        %% Get Output Signals from SRS Simulation

        %Muscle Force and Output Paralyzed Displacement
        force_simulink = out.SRS_Simulation_Force;
        output_displacement_simulink = out.SRS_Simulation_Displacement;
        t_simulink = out.tout;

        %% Input/Output for Model Identification
        
        Zcur = [amplitude_modulation',output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});
 
        %% Calculate Signal to Noise Ratio
        
        %Get Output Noise from SRS Simulation
        output_noise_simulink = out.Output_Noise;
        
        %SNR between clean signal and noise
        signal_to_noise_clean = snr(output_displacement_simulink_clean, output_noise_simulink);
        SRS_noise_snr_clean = [SRS_noise_snr_clean signal_to_noise_clean];
        
        %SNR between clean signal + noise and noise
        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        SRS_noise_snr = [SRS_noise_snr signal_to_noise];
        
        %% Model Identification (LNL, Hammerstein, Wiener, or Linear IRF)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            LNL=lnlbl;
            set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
            'hkTolerance', 0.1, 'nhkMaxIts', 10, 'nhkMaxInner', 5);

            LNL=nlident(LNL,Zcur);
            
            %Evaluate by comparing with the clean Output
            figure(signal)
            [R, V, yp] = nlid_resid(LNL,Zcur_clean);

            SRS_accuracy_identification = [SRS_accuracy_identification V];

            LNL_all = [LNL_all LNL];
            LNL = [];
            SRS_Zcur_ident = [SRS_Zcur_ident Zcur];
            
        elseif Hammerstein_model == true
        
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Hammerstein=nlbl;
            set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
            
            I=Hammerstein{1,2};
            set(I,'nLags',2400,'nSides',1);
            Hammerstein{1,2}=I;

            Hammerstein=nlident(Hammerstein,Zcur);
            
            %Evaluate by comparing with the clean Output
            figure(signal)
            [R, V, yp] = nlid_resid(Hammerstein,Zcur_clean);

            SRS_accuracy_identification = [SRS_accuracy_identification V];

            Hammerstein_all = [Hammerstein_all Hammerstein];
            Hammerstein = [];
            SRS_Zcur_ident = [SRS_Zcur_ident Zcur];

        elseif Weiner_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Weiner = lnbl; %Wiener
            set(Weiner,'idMethod','hk');
            
            I = Weiner{1,1};
            set(I,'nLags',850,'nSides',1); % Set Number of lags and Sides in IRF
            Weiner{1,1} = I;
            
            %Identify the model with the unclean output
            Weiner=nlident(Weiner,Zcur);
            
            % Evaluate by comparing with the clean Output
            figure(signal)
            [R, V, yp] = nlid_resid(Weiner,Zcur_clean);

            SRS_accuracy_identification = [SRS_accuracy_identification V];

            Weiner_all = [Weiner_all Weiner];
            Weiner = [];
            SRS_Zcur_ident = [SRS_Zcur_ident Zcur];

        elseif Linear_IRF_model == true
            
            %Identify a two-sided IRF Model
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);
            
            % Evaluate by comparing with the clean Output
            figure(signal)
            [R, V, yp] = nlid_resid(IRF_model,Zcur_clean);

            SRS_accuracy_identification = [SRS_accuracy_identification V];

            if signal == 1
                IRF_model_phys{noise_level} = IRF_model;
            elseif signal == 2
                IRF_model_PRBS{noise_level} = IRF_model;
            end
            
            IRF_model = [];
            SRS_Zcur_ident = [SRS_Zcur_ident Zcur];

        end
            
    end
    
    set_output_noise_power = set_output_noise_initial;
        
    SRS_accuracy_identification_all(signal,:) = SRS_accuracy_identification;
    SRS_noise_snr_all_clean(signal,:) = SRS_noise_snr_clean;
    SRS_noise_snr_all(signal,:) = SRS_noise_snr;
    SRS_output_noise_power_all(signal,:) = SRS_output_noise_power;
    
    if LNL_model == true
        SRS_models_all = [SRS_models_all;LNL_all];
    elseif Hammerstein_model == true
        SRS_models_all = [SRS_models_all;Hammerstein_all];
    elseif Weiner_model == true
        SRS_models_all = [SRS_models_all;Weiner_all];
    end
    
    SRS_Zcur_clean_all = [SRS_Zcur_clean_all Zcur_clean];
    
    SRS_accuracy_identification = [];
    SRS_noise_snr_clean = [];
    SRS_noise_snr = [];
    SRS_output_noise_power = [];
    LNL_all = [];
    Hammerstein_all = [];
    Weiner_all = [];

end

%% Set Initial Parameters for Model Validation with Physiological Signals

%Noise Parameters for Validation
set_output_noise_power_validation = 0;

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%% Generate Desired Displacement and Amplitude Modulation Signals for Model Validation

for trial = 1:num_trials
    
    t_total = 0:0.001:physiological_stimulus_time;
    time = physiological_stimulus_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 2.1;
    FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
    freq_distribution = random(FrequenciesRandom,10000,1);

    AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
    AmplitudesRandom = AR;
    amp_distribution = random(AmplitudesRandom,10000,1);

    desired_displacement= 0;
    Freq_test = [];
    Pulses_per_interval_test = [];

    for j = 1 : nf    
        t  = 0 : 0.001 : t_interval;         % Time Intervals

        if j == 1
            Freq = FrequenciesRandom_max;
            A = physiological_stimulus_max_amplitude;
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
            stim_frequency = Freq;
            Pulses_per_interval = t_interval/g;
            data(end) = [];
        else
            data = zeros(20000,1);
            stim_frequency = 0;
            Pulses_per_interval = 0;
        end

        %data = A*square(2*pi*f/Fs*t)';        % Generate Square Wave
        desired_displacement = [desired_displacement; data];
        Freq_test = [Freq_test stim_frequency];
        Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
    end

    Pulses_per_interval_total = sum(Pulses_per_interval_test);
    Freq_test_average = sum(Freq_test)/length(Freq_test);
  
    desired_displacement = desired_displacement';

    stim_amplitude = desired_displacement*170;  %mV
    stim_frequency = 50;
    
    %Theoretical Electrical Stimulus
    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    amplitude_modulation = stim_amplitude;
    amplitude_modulation_simulink = [t_total' amplitude_modulation'];
    
    %% Execute SRS Simulation (Simulink Model)
    
    %Set Output Noise as zero
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power_validation')

    %Run Simulink;
    out = sim('SRS_simulation',time);

    %% Get Output Signals from SRS Simulation (Simulink Model)

    %Muscle Force and Output Paralyzed Displacement
    force_simulink = out.SRS_Simulation_Force;
    output_displacement_simulink = out.SRS_Simulation_Displacement;
    t_simulink = out.tout;

    %% Input/Output for Model Validation
    
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

    %% Calculate Signal to Noise Ratio
    
    %Get Output Noise from SRS Simulation
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    
    %% Model Validation (LNL, Hammerstein, Wiener, or Linear IRF)
    
    if LNL_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 

                figure(signal)
                [R, V, yp] = nlid_resid(SRS_models_all(signal,noise_level),Zcur);
                SRS_accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
        
    elseif Hammerstein_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 

                figure(signal)
                [R, V, yp] = nlid_resid(SRS_models_all(signal,noise_level),Zcur);
                SRS_accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
    
    elseif Weiner_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 

                figure(signal)
                [R, V, yp] = nlid_resid(SRS_models_all(signal,noise_level),Zcur);
                SRS_accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
    
    elseif Linear_IRF_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 
                
                if signal == 1
                    figure(signal)
                    [R, V, yp] = nlid_resid(IRF_model_phys{noise_level},Zcur);
                    SRS_accuracy_validation(trial,noise_level,signal) = V;
                elseif signal == 2
                    figure(signal)
                    [R, V, yp] = nlid_resid(IRF_model_PRBS{noise_level},Zcur);
                    SRS_accuracy_validation(trial,noise_level,signal) = V;
                end
                
            end
            
        end

    end
    
end

%% Plot Accuracy (Identification & Validation) vs Noise
figNum = 402;

SRS_accuracy_identification_all = max(0, SRS_accuracy_identification_all);

% figure(figNum)
% figNum = figNum+1;
% plot(noise_snr_all(2,:),accuracy_identification_all(2,:),'LineWidth',3);
% hold on
% plot(noise_snr_all(1,:),accuracy_identification_all(1,:),'LineWidth',3);
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

figure(figNum)
figNum = figNum+1;
% plot(SRS_noise_snr_all_clean(2,:),SRS_accuracy_identification_all(2,:),'LineWidth',3);
hold on
plot(SRS_noise_snr_all_clean(1,:),SRS_accuracy_identification_all(1,:),'LineWidth',3);
hold off
ax = gca;
ax.FontSize = 15;
title('SRS Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

SRS_accuracy_validation = max(0, SRS_accuracy_validation);
SRS_accuracy_validation_mean = mean(SRS_accuracy_validation);
SRS_accuracy_validation_std = std(SRS_accuracy_validation);

% figure(figNum)
% figNum = figNum+1;
% hold on
% plot(noise_snr_all(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
% plot(noise_snr_all(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3);
% 
% plot(noise_snr_all(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
% plot(noise_snr_all(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
% patch([noise_snr_all(2,:) fliplr(noise_snr_all(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)
% 
% plot(noise_snr_all(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
% plot(noise_snr_all(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
% patch([noise_snr_all(1,:) fliplr(noise_snr_all(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Validation Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

figure(figNum)
figNum = figNum+1;
hold on
% plot(SRS_noise_snr_all_clean(2,:),SRS_accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),SRS_accuracy_validation_mean(:,:,1),'LineWidth',3);

% plot(SRS_noise_snr_all_clean(2,:),min(SRS_accuracy_validation_mean(:,:,2)+SRS_accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
% plot(SRS_noise_snr_all_clean(2,:),max(SRS_accuracy_validation_mean(:,:,2)-SRS_accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
% patch([SRS_noise_snr_all_clean(2,:) fliplr(SRS_noise_snr_all_clean(2,:))], [min(SRS_accuracy_validation_mean(:,:,2)+SRS_accuracy_validation_std(:,:,2),100) fliplr(max(SRS_accuracy_validation_mean(:,:,2)-SRS_accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(SRS_accuracy_validation_mean(:,:,1)+SRS_accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(SRS_accuracy_validation_mean(:,:,1)-SRS_accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(SRS_accuracy_validation_mean(:,:,1)+SRS_accuracy_validation_std(:,:,1),100) fliplr(max(SRS_accuracy_validation_mean(:,:,1)-SRS_accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('SRS Validation Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on


%% Inverse SRS

prompt2 = 'Inverse SRS Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Hamm]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'LNL') & ~strcmp(str2,'Hamm') & ~strcmp(str2,'Wiener') & ~strcmp(str2,'IRF') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Hamm';
end
SRS_inverse_model_structure = str2;

prompt3 = 'Number of Validation Trials? 1-30 [1]: ';
str3 = input(prompt3);
if str3<1 | str3>30
    disp('Invalid Input')
    return
elseif isempty(str3)
    str3 = 1;
end

%% Set Initial Parameters

figNum = 100;
Fs = 1000;

%Simulated Input (PRBS) Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;      %Can either be constant or variable amplitude
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20;            %PRBS Signal Amplitude (mm)

%Number of Signals (Identification and Validation)
num_signals = str3+1;
%num_signals = str3+noise_level_iters; %If we want to use a new signal for
%each inverse SRS identification

%Inverse Model Structure
SRS_inverse_LNL = false;
SRS_inverse_Hamm = false;
SRS_inverse_Weiner = false;
SRS_inverse_IRF = false;

if strcmp(str2,'LNL')
    SRS_inverse_LNL = true;
elseif strcmp(str2,'Hamm')
    SRS_inverse_Hamm = true;
elseif strcmp(str2,'Wiener')
    SRS_inverse_Weiner = true;
elseif strcmp(str2,'IRF')
    SRS_inverse_IRF = true;
end

inverse_SRS_identification_accuracy = [];
SRS_inverse_all = [];

%% Simulated Input (PRBS Input)

simulated_input = [];

for signal = 1:num_signals

    t_total = 0:0.001:PRBS_stimulus_time;
    time = PRBS_stimulus_time;

    A = [0];                                    %Intialize amplitude
    if variable_amplitude == true   
        for k = 1:N
            if k == 1
                R = PRBS_amplitude;
            else
                R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
            end

            for j = 1:M
                A = [A R];
            end
        end
    else
        A = PRBS_amplitude;              %Else set as Constant Amplitude
    end

    Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

    %Specify the clock period of the signal as 1 sample. 
    %That is, the signal value can change at each time step. 
    %For PRBS signals, the clock period is specified in Band = [0 B], 
    %where B is the inverse of the required clock period
    %(Must be less than 1)
    Band = [0 0.01];

    %Generate a nonperiodic PRBS of length time.
    u = idinput(time*1000+1,'prbs',Band,Range);

    %Create an iddata object from the generated signal. 
    %For this example, specify the sample time as 0.001 second.
    u = iddata([],u,0.001);

    U = (u.InputData)';
    desired_displacement = A.*U;

    stim_frequency = 50;
    stim_amplitude = desired_displacement*170;

    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    simulated_input(:,signal) = stim_amplitude;

end

%% Simulate the reponse of the SRS model to the Input Realizations

%Runs through each SRS model and identifies all inverse SRS models
%Uses the same simulated input/output each model
for noise_level = 1:noise_level_iters

    %Simulated Outputs of the SRS Model
    simulated_output = [];
    simulated_outputs = [];

    %Runs through all input realizations
    for signal = 1:num_signals

        simulated_output = nlsim(SRS_models_all(1,noise_level),simulated_input(:,signal));
        set(simulated_output, 'domainIncr',0.001);

        simulated_output = max(simulated_output,0);
        simulated_outputs(:,signal) = simulated_output.dataSet;

    end

    %% Identify a model between a simulated output (as input) and a simulated input (as output)

    %Input/Output for Inverse SRS Identification
    Zcur_simulated = [simulated_outputs(1000:end,1) simulated_input(1000:end,1)];
    Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

%     %plot the Input/Output
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,1,1)
%     plot(t_total(1,1:end-999),simulated_input(1000:end,1))
%     ax = gca;
%     ax.FontSize = 15;
%     title('Simulated Amplitude Modulation, A(t)','Fontsize',20)
%     ylabel('Amplitude (V)','Fontsize',18)
%     grid on
% 
%     subplot(2,1,2)
%     plot(t_total(1,1:end-999),simulated_outputs(1000:end,1))
%     ax = gca;
%     ax.FontSize = 15;
%     title('Simulated Paralyzed Displacement, Pos_P(t)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',18)
%     xlabel('Time (s)','Fontsize',18)
%     grid on

    %Outputs to other input realizations 
    validation_output = [];

    for signal = 2:num_signals

        Validation_Output = nldat(simulated_outputs(1000:end,signal));
        set(Validation_Output,'domainIncr',0.001);

        validation_output(:,signal-1) = Validation_Output.dataSet;

    end

    %% Identify the Inverse SRS

    %Number of Lags of Inverse SRS IRF(s)
    nLags = 400;

    if SRS_inverse_LNL == true

        SRS_inverse = lnlbl;  %LNL Model
        set(SRS_inverse,'idMethod','hk','hkTolerance', 0.1,...
            'nhkMaxIts', 4, 'nhkMaxInner', 4);

        I1 = irf;
        set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
        P = polynom;
        I3 = irf;
        set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
        SRS_inverse.elements = {I1 P I3};

        SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

    elseif SRS_inverse_Hamm == true

        SRS_inverse = nlbl;  %Hammerstein
        set(SRS_inverse,'idMethod','hk','displayFlag',true,'threshNSE',.001);
        I2 = irf;
        set(I2,'nLags',nLags, 'nSides', 2,'domainIncr',0.001); % Set number of lags and Sides in IRF
        SRS_inverse{1,2} = I2;

        SRS_inverse = nlident(SRS_inverse,Zcur_simulated);
        
        SRS_inverse_IRF = SRS_inverse{1,2};
        set(SRS_inverse_IRF, 'domainIncr', 1.0e-3);
        SRS_inverse{1,2} = SRS_inverse_IRF;

    elseif SRS_inverse_Weiner == true

        SRS_inverse = lnbl; %Wiener
        set(SRS_inverse,'idMethod','hk');
        I1 = irf;
        set(I1,'nLags',nLags,'nSides',2,'domainIncr',0.001); % Set Number of lags and Sides in IRF
        SRS_inverse{1,1} = I1;

        SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

    elseif SRS_inverse_IRF == true

        %Two-Sided IRF
        SRS_inverse = irf(Zcur_simulated,'nLags',nLags,'nSides',2,'domainIncr',0.001);

    end

    figure(figNum)
    figNum = figNum+1;
    [R, V, yp] = nlid_resid(SRS_inverse,Zcur_simulated);
    
    inverse_SRS_identification_accuracy = [inverse_SRS_identification_accuracy V];
    
    %% Inverse SRS Validation

    control_system_test_output_double = [];
    inverse_SRS_validation_accuracy = [];

    for signal = 1:num_signals-1

        control_system_test_output = nlsim(SRS_inverse,validation_output(:,signal));
        set(control_system_test_output,'domainIncr',0.001);
        control_system_test_output_double(:,signal) = control_system_test_output.dataSet;

        inverse_SRS_validation_accuracy(signal,noise_level) = vaf(simulated_input(1999:end-1000,signal+1),control_system_test_output_double(1000:end-1000,signal));

    end
    
    SRS_inverse_all = [SRS_inverse_all;SRS_inverse];
    SRS_inverse = [];
    
end

%% Plot
figNum = 404;

inverse_SRS_identification_accuracy = max(0, inverse_SRS_identification_accuracy);

figure(figNum)
figNum = figNum+1;
%plot(noise_snr_all_clean(2,:),accuracy_identification_all(2,:),'LineWidth',3);
hold on
plot(SRS_noise_snr_all_clean(1,:),inverse_SRS_identification_accuracy(1,:),'LineWidth',3);
hold off
ax = gca;
ax.FontSize = 15;
title('Inverse SRS Identification Accuracy vs Output Noise','Fontsize',24);
%legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

%Mean and Standard Deviation of Inverse SRS Validation Accuracy
inverse_SRS_validation_accuracy = max(0, inverse_SRS_validation_accuracy);
inverse_SRS_validation_accuracy_mean = mean(inverse_SRS_validation_accuracy);
inverse_SRS_validation_accuracy_std = std(inverse_SRS_validation_accuracy);

%disp(['Validation Accuracy Mean: ' num2str(round(validation_accuracy_mean,1)) '%'])
%disp(['Validation Accuracy Std: ' num2str(round(validation_accuracy_std,1)) '%'])

figure(figNum)
figNum = figNum+1;
hold on
%plot(noise_snr_all_clean(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),inverse_SRS_validation_accuracy_mean(1,:),'LineWidth',3);

%plot(noise_snr_all_clean(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
%plot(noise_snr_all_clean(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
%patch([noise_snr_all_clean(2,:) fliplr(noise_snr_all_clean(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(inverse_SRS_validation_accuracy_mean(1,:)+inverse_SRS_validation_accuracy_std(1,:),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(inverse_SRS_validation_accuracy_mean(1,:)-inverse_SRS_validation_accuracy_std(1,:),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(inverse_SRS_validation_accuracy_mean(1,:)+inverse_SRS_validation_accuracy_std(1,:),100) fliplr(max(inverse_SRS_validation_accuracy_mean(1,:)-inverse_SRS_validation_accuracy_std(1,:),0))], 'b','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('Inverse SRS Validation Accuracy vs Output Noise','Fontsize',24);
%legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

%% FRCS Execute and Evaluation

prompt1 = 'Type of FRCS Input? PRBS/Physiological(Phys) [Phys]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'PRBS') & ~strcmp(str1,'Phys') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Phys';
end

if strcmp(str1,'PRBS')
    PRBS_signal = true;
elseif strcmp(str1,'Phys')
    PRBS_signal = false;
end 
FRCS_input_type = str1;

prompt2 = 'Number of Validation Trials? 1-100 [1]: ';
str2 = input(prompt2);
if str2<1 | str2>100
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 1;
end

% prompt3 = 'Number of FRCS Versions? 1-20 [1]: ';
% str3 = input(prompt3);
% if str3<1 | str3>20
%     disp('Invalid Input')
%     return
% elseif isempty(str3)
%     str3 = 1;
% end

%% Set initial Parameters

Fs = 1000; 
Nfft = 10000;
figNum = 406;

%ERS and SRS Model
% if strcmp(SRS_inverse_input_type,'PRBS')
%     EMG_response_system = ERS_all(1);
% elseif strcmp(SRS_inverse_input_type,'Phys')
%     EMG_response_system = ERS_all(2);
% end

EMG_response_system = ERS_models_all;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10;                            %PRBS Amplitude (mm)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_stimulus_max_amplitude = 0.01;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%Number of Validation Trials and FRCS Versions
num_val_trials = str2;
%num_FRCS_models = str3;

FRCS_validation_accuracy = [];

%% Generate the Desired Displacement Signal for ERS Simulation

for noise_level = 1:noise_level_iters
    
    %Tests the FRCS on multiple signal realizations
    %Only plots the results of the first realization
    for trial = 1:num_val_trials

        %This is used to get an EMG Signal for Input to the FRCS
        if PRBS_signal == true

            t_total = 0:0.001:PRBS_movement_time;
            time = PRBS_movement_time;

            A = 0;                                      %Intialize amplitude
            if variable_amplitude == true   
                for k = 1:N
                    if k == 1
                        R = PRBS_amplitude;
                    else
                        R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
                    end
                    for j = 1:M
                        A = [A R];
                    end
                end
            else
                A = PRBS_amplitude;              %Else set as Constant Amplitude
            end

            Range = [0,0.001]; %Specify what the single-channel PRBS value switches between

            %Specify the clock period of the signal as 1 sample. 
            %That is, the signal value can change at each time step. 
            %For PRBS signals, the clock period is specified in Band = [0 B], 
            %where B is the inverse of the required clock period
            %(Must be less than 1)
            Band = [0 0.01];

            %Generate a nonperiodic PRBS of length time
            u = idinput(time*1000+1,'prbs',Band,Range);

            %Create an iddata object from the generated signal. 
            %For this example, specify the sample time as 0.001 second.
            u = iddata([],u,0.001);

            U = (u.InputData)';
            desired_displacement = A.*U;

        else

            t_total = 0:0.001:physiological_movement_time;
            time = physiological_movement_time;

            FR = makedist('Normal','mu',fr,'sigma',sig);
            FrequenciesRandom_max = 1.8;
            FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
            freq_distribution = random(FrequenciesRandom,10000,1);

            AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
            AmplitudesRandom = AR;
            amp_distribution = random(AmplitudesRandom,10000,1);

            desired_displacement= 0;
            Freq_test = [];
            Pulses_per_interval_test = [];

            for j = 1 : nf    
                t  = 0 : 0.001 : t_interval;         % Time Intervals

                if j == 1
                    Freq = FrequenciesRandom_max;
                    A = physiological_stimulus_max_amplitude;
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
                    stim_frequency = Freq;
                    Pulses_per_interval = t_interval/g;
                    data(end) = [];
                else
                    data = zeros(20000,1);
                    stim_frequency = 0;
                    Pulses_per_interval = 0;
                end

                desired_displacement = [desired_displacement; data];
                Freq_test = [Freq_test stim_frequency];
                Pulses_per_interval_test = [Pulses_per_interval_test Pulses_per_interval];
            end

            Pulses_per_interval_total = sum(Pulses_per_interval_test);
            Freq_test_average = sum(Freq_test)/length(Freq_test);

            desired_displacement = desired_displacement';

        end

        %% Generate Neural Input Command Signal for ERS Simulation

        %Neural Input Parameters are based on desired displacement amplitude
        Amplitude = desired_displacement*100;  %mV
        Frequency = desired_displacement*14000;  %Hz

        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';
        neural_simulink = [t_total' neural];

        %% Generate EMG Signal use EMG Simulation (Simulink Model)

        %Run Simulink;
        out = sim('ERS_simulation',time);

        %Get Output from Simulink Model (EMG)
        EMG_simulink = out.ERS_Simulation_EMG;
        t_simulink = out.tout;

        EMG_simulink = nldat(EMG_simulink);
        set(EMG_simulink, 'domainIncr',0.001);

        EMG_double = EMG_simulink.dataSet;

        %% Simulate FRCS Response to EMG Input

        %Set domain increments of EMG response model
%         EMG_response_system_IRF = EMG_response_system{1,2};
%         set(EMG_response_system_IRF, 'domainIncr', 1.0e-3);
%         EMG_response_system{1,2} = EMG_response_system_IRF;

        %Simulate Response of ERS to EMG Input
        control_system_displacement = nlsim(ERS_models_all(1,noise_level),EMG_simulink);
        set(control_system_displacement,'domainIncr',0.001);

        %Simulate Response of Inverse SRS to Healthy Displacement (Output of ERS)
        control_system_output = nlsim(SRS_inverse_all(1,noise_level),control_system_displacement);
        set(control_system_output,'domainIncr',0.001);

        control_system_displacement_double = control_system_displacement.dataSet;
        control_system_output_double = control_system_output.dataSet;

        %% Plot Inputs and Output of Facial Reanimation Control System

    %     if trial == 1
    %         figure(figNum)
    %         figNum = figNum+1;
    %         sgtitle('Inputs/Outputs of FRCS','Fontsize',14)
    % 
    %         subplot(3,2,1)
    %         plot(t_total,neural)
    %         title('Neural Command Input for ERS Simulation','Fontsize',14)
    %         xlabel('Time (s)','Fontsize',14)
    %         ylabel('Amplitude (% of MUs)', 'Fontsize', 14)
    %         grid on
    % 
    %         subplot(3,2,2)
    %         plot(t_total,EMG_double)
    %         title('EMG Input from ERS Simulation','Fontsize',14)
    %         xlabel('Time (s)','Fontsize',14)
    %         ylabel('Amplitude (V)', 'Fontsize', 14)
    %         grid on
    % 
    %         subplot(3,2,[3 4])
    %         plot(t_total(1,1:end-1999),control_system_displacement_double(1000:end-1000,1));
    %         title('ERS Output (Healthy Side Predicted Displacement)', 'Fontsize', 14)
    %         ylabel('Displacement (m)', 'Fontsize', 14)
    %         grid on
    % 
    %         subplot(3,2,[5 6])
    %         plot(t_total(1,1:end-1999),control_system_output_double(1000:end-1000,:));
    %         title('FRCS Output (Amplitude Modulation)', 'Fontsize', 14)
    %         xlabel('Time (s)','Fontsize',14)
    %         ylabel('Amplitude', 'Fontsize', 14)
    %         grid on
    %     end

        %% Run the Control System Output (Amplitude Modulation) through SRS Simulation

        predicted_healthy_displacement = control_system_displacement_double(1000:end-1000,1);

        %Execute SRS Simulation with amplitude modulation input
        amplitude_modulation = control_system_output_double;
        amplitude_modulation_simulink = [(t_total(1,1:end-1999))' amplitude_modulation(1000:end-1000,:)];

        %Set Output Noise as zero
        set_output_noise_power = 0;
        set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')

        %Run Simulink;
        time = length(t_total(1,1:end-2000))/1000;
        out = sim('SRS_simulation',time);

        %Get Output signals from SRS Simulation
        electrical_stimulus = out.SRS_Simulation_Stimulus;
        predicted_paralyzed_displacement = out.SRS_Simulation_Displacement;
        t_simulink = out.tout;

        Variance = vaf(predicted_healthy_displacement,predicted_paralyzed_displacement); 
        FRCS_validation_accuracy(trial,noise_level) = Variance;

        %To save the displacement data
        FRCS_desired_displacement_all(:,trial,noise_level) = desired_displacement;
        FRCS_predicted_healthy_displacement_all(:,trial,noise_level) = predicted_healthy_displacement;
        FRCS_predicted_paralyzed_displacement_all(:,trial,noise_level) = predicted_paralyzed_displacement;

    %     if trial == 1
    %         %plot the stimulus that is created based on the amplitude modulation
    %         figure(figNum)
    %         figNum = figNum+1;
    %         plot(t_total(1,1:end-1999),electrical_stimulus)
    %         title('Electrical Stimulus Signal based on FRCS Amplitude Modulation Output')
    %         ylabel('Voltage (V)')
    %         xlabel('Time (s)')
    %         grid on
    % 
    %         %Compare the predicted healthy displacement and predicted paralyzed
    %         %displacement with the desired displacement
    %         figure(figNum)
    %         figNum = figNum+1;
    %         subplot(4,1,1)
    %         plot(t_total(1,1:end-1999),desired_displacement(:,1000:end-1000))
    %         ax = gca;
    %         ax.FontSize = 10;
    %         title('(a) Desired Displacement', 'Fontsize', 12)
    %         ylabel('Displacement (m)', 'Fontsize', 12)
    %         grid on
    % 
    %         subplot(4,1,2)
    %         plot(t_total(1,1:end-1999),predicted_healthy_displacement)
    %         ax = gca;
    %         ax.FontSize = 10;
    %         title('(b) Predicted Healthy Displacement, Pos_H(t)', 'Fontsize', 12)
    %         ylabel('Displacement (m)', 'Fontsize', 12)
    %         grid on
    % 
    %         subplot(4,1,3)
    %         plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
    %         ax = gca;
    %         ax.FontSize = 10;
    %         title('(c) Predicted Paralyzed Displacement, Pos_P(t)', 'Fontsize', 12)
    %         ylabel('Displacement (m)', 'Fontsize', 12)
    %         grid on
    % 
    %         subplot(4,1,4)
    %         hold on
    %         plot(t_total(1,1:end-1999),predicted_healthy_displacement)
    %         plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
    %         hold off
    %         ax = gca;
    %         ax.FontSize = 10;
    %         title(['(d) Superimposed, VAF = ' num2str(round(Variance,1)) '%'], 'Fontsize', 12)
    %         xlabel('Time (s)', 'Fontsize', 12)
    %         ylabel('Displacement (m)', 'Fontsize', 12)
    %         legend('Pos_H(t)', 'Pos_P(t)','Fontsize',10)
    %         grid on
    % 
    %         %PDF and Spectrum of Residuals (between predicted healthy displacement and
    %         %predicted paralyzed displacement)
    %         residuals = predicted_healthy_displacement - predicted_paralyzed_displacement;
    % 
    %         Nfft = 10000;
    %         residuals_zero = residuals - mean(residuals);
    %         [Prr,f] = pwelch(residuals_zero,Nfft,[],Nfft,Fs);
    % 
    %         figure(figNum)
    %         figNum = figNum+1;
    %         sgtitle('Residuals between Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)
    % 
    %         subplot(2,2,[1 2])
    %         plot(t_total(1,1:end-3998),residuals(1000:end-1000,:))
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(a) Residuals','Fontsize',16)
    %         ylabel('Displacement (m)','Fontsize',16)
    %         xlabel('Time (s)','Fontsize',16)
    %         grid on
    % 
    %         subplot(2,2,3)
    %         histogram(residuals(1000:end-1000,:))
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(b) Residual Distribution','Fontsize',16)
    %         xlabel('Displacement (m)','Fontsize',16)
    %         ylabel('Density','Fontsize',16)
    %         grid on
    % 
    %         subplot(2,2,4)
    %         semilogy(f(1:200,:),Prr(1:200,:))
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(c) Spectrum of Residuals','Fontsize',16)
    %         ylabel('PSD (log)','Fontsize',16);
    %         xlabel('Frequency (Hz)','Fontsize',16);
    %         grid on;
    % 
    %         %Spectrum of Predicted Healthy Displacement and Predicted Paralyzed
    %         %Displacement
    %         predicted_healthy_displacement_zero = predicted_healthy_displacement - mean(predicted_healthy_displacement);
    %         predicted_paralyzed_displacement_zero = predicted_paralyzed_displacement - mean(predicted_paralyzed_displacement);
    % 
    %         Nfft = 10000;
    %         [Pxx,~] = pwelch(predicted_healthy_displacement_zero,Nfft,[],Nfft,Fs);
    %         [Pyy,f] = pwelch(predicted_paralyzed_displacement_zero,Nfft,[],Nfft,Fs);
    % 
    %         figure(figNum)
    %         figNum = figNum+1;
    %         subplot(2,2,1),plot(t_total(1,1:end-1999),predicted_healthy_displacement)
    %         title('Predicted Healthy Displacement, Pos_H(t)','Fontsize',14)
    %         xlabel('Time (s)','Fontsize',16)
    %         ylabel('Displacement (m)','Fontsize',16)
    %         grid on
    % 
    %         subplot(2,2,2),plot(t_total(1,1:end-1999),predicted_paralyzed_displacement)
    %         title('Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)
    %         xlabel('Time (s)','Fontsize',16)
    %         ylabel('Displacement (m)','Fontsize',16)
    %         grid on
    % 
    %         subplot(2,2,3),semilogy(f(1:200,:),Pxx(1:200,:));
    %         title('Predicted Healthy Displacement Spectrum','Fontsize',14);
    %         ylabel('PSD (log scale)','Fontsize',16);
    %         xlabel('Frequency (Hz)','Fontsize',16);
    %         grid on;
    % 
    %         subplot(2,2,4),semilogy(f(1:200,:),Pyy(1:200,:));
    %         title('Predicted Paralyzed Displacement Spectrum','Fontsize',14);
    %         ylabel('PSD','Fontsize',16);
    %         xlabel('Frequency (Hz)','Fontsize',16);
    %         grid on;
    % 
    %         %Auto-Correlation and Cross Correlation of Predicted Healthy Displacement
    %         %and Predicted Paralyzed Displacement
    %         figure(figNum)
    %         figNum = figNum+1;
    %         sgtitle('Correlation of Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)
    % 
    %         maxlag = floor(length(predicted_healthy_displacement_zero)*0.05);
    %         [Rxx,lags] = xcorr(predicted_healthy_displacement_zero,maxlag,'coeff');
    %         lags = lags/Fs;
    %         subplot(3,1,1),plot(lags,Rxx)
    %         title('(a) Auto-Correlation of Predicted Healthy Displacement','Fontsize',16)
    %         ylabel('Correlation','Fontsize',16)
    %         grid on
    % 
    %         [Ryy,lags] = xcorr(predicted_paralyzed_displacement_zero,maxlag,'coeff');
    %         lags = lags/Fs;
    %         subplot(3,1,2),plot(lags,Ryy)
    %         title('(b) Auto-Correlation of Predicted Paralyzed Displacement','Fontsize',16)
    %         ylabel('Correlation','Fontsize',16)
    %         grid on
    % 
    %         [Rxy,lags] = xcorr(predicted_paralyzed_displacement_zero,predicted_healthy_displacement_zero,maxlag,'coeff');
    %         lags = lags/Fs;
    %         subplot(3,1,3),plot(lags,Rxy)
    %         title('(c) Cross-Correlation','Fontsize',16)
    %         xlabel('Time (s)','Fontsize',16)
    %         ylabel('Correlation','Fontsize',16)
    %         grid on
    % 
    %         % Frequency Response between Predicted Healthy Side Displacement and
    %         % Predicted Paralyzed Displacement
    %         %Transfer Function
    %         Nfft = 2000;
    %         txy = tfestimate(predicted_healthy_displacement,predicted_paralyzed_displacement,Nfft,[],Nfft,Fs);
    %         G = abs(txy);
    %         G = 20*log10(G);
    %         P = angle(txy);
    %         P = P*(180/pi);
    %         [Cxy, f] = mscohere(predicted_healthy_displacement,predicted_paralyzed_displacement,Nfft,[],Nfft,Fs);
    % 
    %         figure(figNum)
    %         figNum = figNum+1;
    %         sgtitle('Frequency Response of Predicted Healthy Displacement, Pos_H(t), and Predicted Paralyzed Displacement, Pos_P(t)','Fontsize',14)
    % 
    %         subplot(3,1,1),plot(f(1:20,:),G(1:20,:));
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(a) Gain','Fontsize',16);
    %         ylabel('Gain (dB)','Fontsize',16);
    %         grid on;
    % 
    %         subplot(3,1,2),plot(f(1:20,:),P(1:20,:));
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(b) Phase Shift','Fontsize',16);
    %         ylabel('Phase (degrees)','Fontsize',16);
    %         grid on;
    % 
    %         subplot(3,1,3),plot(f(1:20,:),Cxy(1:20,:));
    %         ax = gca;
    %         ax.FontSize = 14;
    %         title('(c) Coherence','Fontsize',16);
    %         ylabel('Coherence','Fontsize',16);
    %         xlabel('Frequency (Hz)','Fontsize',16);
    %         grid on;
    % 
    %         figure(figNum)
    %         figNum = figNum+1;
    %         plot(f(1:20,:),Cxy(1:20,:));
    %         ax = gca;
    %         ax.FontSize = 16;
    %         title('Coherence between Pos_H(t) and Pos_P(t)','Fontsize',20);
    %         ylabel('Coherence','Fontsize',20);
    %         xlabel('Frequency (Hz)','Fontsize',20);
    %         grid on;
    % 
    %     end

    end
    
end

%% Plot
figNum = 408;

FRCS_validation_accuracy_mean = [];
FRCS_validation_accuracy_std = [];

%Mean and Standard Deviation of Inverse SRS Validation Accuracy
FRCS_validation_accuracy = max(0, FRCS_validation_accuracy);
FRCS_validation_accuracy_mean = mean(FRCS_validation_accuracy,1);
FRCS_validation_accuracy_std = std(FRCS_validation_accuracy,0,1);

%disp(['Validation Accuracy Mean: ' num2str(round(validation_accuracy_mean,1)) '%'])
%disp(['Validation Accuracy Std: ' num2str(round(validation_accuracy_std,1)) '%'])

figure(figNum)
figNum = figNum+1;
hold on
%plot(noise_snr_all_clean(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),FRCS_validation_accuracy_mean(1,:),'LineWidth',3);

%plot(noise_snr_all_clean(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
%plot(noise_snr_all_clean(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
%patch([noise_snr_all_clean(2,:) fliplr(noise_snr_all_clean(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(FRCS_validation_accuracy_mean(1,:)+FRCS_validation_accuracy_std(1,:),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(FRCS_validation_accuracy_mean(1,:)-FRCS_validation_accuracy_std(1,:),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(FRCS_validation_accuracy_mean(1,:)+FRCS_validation_accuracy_std(1,:),100) fliplr(max(FRCS_validation_accuracy_mean(1,:)-FRCS_validation_accuracy_std(1,:),0))], 'b','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('FRCS Accuracy vs Output Noise','Fontsize',24);
%legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

%% Plot all
figNum = 410;

figure(figNum)
figNum = figNum+1;
hold on
plot(ERS_noise_snr_all_clean(1,:),ERS_accuracy_validation_mean(:,:,1),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),SRS_accuracy_validation_mean(:,:,1),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),inverse_SRS_validation_accuracy_mean(1,:),'LineWidth',3);
plot(SRS_noise_snr_all_clean(1,:),FRCS_validation_accuracy_mean(1,:),'LineWidth',3);

plot(ERS_noise_snr_all_clean(1,:),min(ERS_accuracy_validation_mean(:,:,1)+ERS_accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(ERS_noise_snr_all_clean(1,:),max(ERS_accuracy_validation_mean(:,:,1)-ERS_accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([ERS_noise_snr_all_clean(1,:) fliplr(ERS_noise_snr_all_clean(1,:))], [min(ERS_accuracy_validation_mean(:,:,1)+ERS_accuracy_validation_std(:,:,1),100) fliplr(max(ERS_accuracy_validation_mean(:,:,1)-ERS_accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(SRS_accuracy_validation_mean(:,:,1)+SRS_accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(SRS_accuracy_validation_mean(:,:,1)-SRS_accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(SRS_accuracy_validation_mean(:,:,1)+SRS_accuracy_validation_std(:,:,1),100) fliplr(max(SRS_accuracy_validation_mean(:,:,1)-SRS_accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(inverse_SRS_validation_accuracy_mean(1,:)+inverse_SRS_validation_accuracy_std(1,:),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(inverse_SRS_validation_accuracy_mean(1,:)-inverse_SRS_validation_accuracy_std(1,:),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(inverse_SRS_validation_accuracy_mean(1,:)+inverse_SRS_validation_accuracy_std(1,:),100) fliplr(max(inverse_SRS_validation_accuracy_mean(1,:)-inverse_SRS_validation_accuracy_std(1,:),0))], 'b','FaceAlpha',0.2)

plot(SRS_noise_snr_all_clean(1,:),min(FRCS_validation_accuracy_mean(1,:)+FRCS_validation_accuracy_std(1,:),100),'LineStyle','none','LineWidth',2)
plot(SRS_noise_snr_all_clean(1,:),max(FRCS_validation_accuracy_mean(1,:)-FRCS_validation_accuracy_std(1,:),0),'LineStyle','none','LineWidth',2)
patch([SRS_noise_snr_all_clean(1,:) fliplr(SRS_noise_snr_all_clean(1,:))], [min(FRCS_validation_accuracy_mean(1,:)+FRCS_validation_accuracy_std(1,:),100) fliplr(max(FRCS_validation_accuracy_mean(1,:)-FRCS_validation_accuracy_std(1,:),0))], 'b','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('FRCS Accuracy vs Output Noise','Fontsize',24);
legend('ERS','SRS','Inverse SRS','FRCS','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

tEnd = toc(tStart)/60
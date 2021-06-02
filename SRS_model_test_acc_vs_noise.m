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

clc
clear all

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

tStart = tic;

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

accuracy_identification = [];
accuracy_identification_all = [];
accuracy_validation = [];

noise_snr_clean = [];
noise_snr = [];
noise_snr_all_clean = [];
noise_snr_all = [];
output_noise_power = [];
output_noise_power_all = [];

%Initialize Models
LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];

models_all = [];

Zcur_ident = [];
Zcur_clean_all = [];

%% Num trials and Signal Types
num_trials = str2;
signals = 2;

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
        output_noise_power = [output_noise_power set_output_noise_power];

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
        noise_snr_clean = [noise_snr_clean signal_to_noise_clean];
        
        %SNR between clean signal + noise and noise
        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];
        
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

            accuracy_identification = [accuracy_identification V];

            LNL_all = [LNL_all LNL];
            LNL = [];
            Zcur_ident = [Zcur_ident Zcur];
            
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

            accuracy_identification = [accuracy_identification V];

            Hammerstein_all = [Hammerstein_all Hammerstein];
            Hammerstein = [];
            Zcur_ident = [Zcur_ident Zcur];

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

            accuracy_identification = [accuracy_identification V];

            Weiner_all = [Weiner_all Weiner];
            Weiner = [];
            Zcur_ident = [Zcur_ident Zcur];

        elseif Linear_IRF_model == true
            
            %Identify a two-sided IRF Model
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);
            
            % Evaluate by comparing with the clean Output
            figure(signal)
            [R, V, yp] = nlid_resid(IRF_model,Zcur_clean);

            accuracy_identification = [accuracy_identification V];

            if signal == 1
                IRF_model_phys{noise_level} = IRF_model;
            elseif signal == 2
                IRF_model_PRBS{noise_level} = IRF_model;
            end
            
            IRF_model = [];
            Zcur_ident = [Zcur_ident Zcur];

        end
            
    end
    
    set_output_noise_power = set_output_noise_initial;
        
    accuracy_identification_all(signal,:) = accuracy_identification;
    noise_snr_all_clean(signal,:) = noise_snr_clean;
    noise_snr_all(signal,:) = noise_snr;
    output_noise_power_all(signal,:) = output_noise_power;
    
    if LNL_model == true
        models_all = [models_all;LNL_all];
    elseif Hammerstein_model == true
        models_all = [models_all;Hammerstein_all];
    elseif Weiner_model == true
        models_all = [models_all;Weiner_all];
    end
    
    Zcur_clean_all = [Zcur_clean_all Zcur_clean];
    
    accuracy_identification = [];
    noise_snr_clean = [];
    noise_snr = [];
    output_noise_power = [];
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
                [R, V, yp] = nlid_resid(models_all(signal,noise_level),Zcur);
                accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
        
    elseif Hammerstein_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 

                figure(signal)
                [R, V, yp] = nlid_resid(models_all(signal,noise_level),Zcur);
                accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
    
    elseif Weiner_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 

                figure(signal)
                [R, V, yp] = nlid_resid(models_all(signal,noise_level),Zcur);
                accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end
    
    elseif Linear_IRF_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters 
                
                if signal == 1
                    figure(signal)
                    [R, V, yp] = nlid_resid(IRF_model_phys{noise_level},Zcur);
                    accuracy_validation(trial,noise_level,signal) = V;
                elseif signal == 2
                    figure(signal)
                    [R, V, yp] = nlid_resid(IRF_model_PRBS{noise_level},Zcur);
                    accuracy_validation(trial,noise_level,signal) = V;
                end
                
            end
            
        end

    end
    
end

%% Plot Accuracy (Identification & Validation) vs Noise
figNum = 400;

accuracy_identification_all = max(0, accuracy_identification_all);

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
plot(noise_snr_all_clean(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
plot(noise_snr_all_clean(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3);

plot(noise_snr_all_clean(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(2,:) fliplr(noise_snr_all_clean(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(noise_snr_all_clean(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(1,:) fliplr(noise_snr_all_clean(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)
hold off
ax = gca;
ax.FontSize = 15;
title('Validation Accuracy vs Output Noise','Fontsize',24);
legend('PRBS','Physiological','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on

tEnd = toc(tStart)/60
%% SRS Model Test for Number of Movement Pulses of "Physiological" Signal used for Identification

%Evaluates the Model based on the number of movement pulses of the 
%"Physiological" Input used for identification. Plots the identifcation
%accuracy vs number of movement pulses and validation accuracy vs number of
%movement pulses

%The script generates a set of physiological desired displacement signals, 
%starting at 1 movement pulse and iteratively increases the number of movements 
%by 1 for subsequent signals. A specified number of Models are identified for 
%each of these signals. A set of these models are then validated with 
%typical "physiological" movements for a specified number of trials to
%calculate the validation VAF mean and std.

%When running the script, you need to provide the following input:
% 1. Type of Model Structure? LNL/Hammerstein/Wiener/IRF
%       The Model structure used for the Identified SRS (Default is Wiener)
% 2. Maximum number of Movement Pulses?
%       The signals used for identification will start at 1 movement pulse
%       and iteratively increase by 1 until the maximum number of movement
%       pulse is reached
% 3. Number of Identification Trials?
%       Number of model identification trials to calculate VAF mean and
%       std for model identification. Increasing this value has a large 
%       impact on how long the script needs to run 
% 4. Number of Validation Trials?
%       Number of  model validation trials to calculate VAF mean and std

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

prompt2 = 'Maximum Number of Movement Pulses? 1-50 [30]: ';
str2 = input(prompt2);
if str2<1 | str2>50
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 30;
end

prompt3 = 'Number of Identification Trials? 1-20 [1]: ';
str3 = input(prompt3);
if str3<1 | str3>20
    disp('Invalid Input')
    return
elseif isempty(str3)
    str3 = 1;
end

prompt4 = 'Number of Validation Trials? 1-50 [30]: ';
str4 = input(prompt4);
if str4<1 | str4>50
    disp('Invalid Input')
    return
elseif isempty(str4)
    str4 = 30;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                                       %Width of signal pulse (seconds)
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%Type of Model
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

%Initialize Models
SRS_models_temp = [];
SRS_models = [];

Zcur_all_temp = [];
Zcur_all = [];

identification_accuracy = [];
validation_accuracy = [];

%% Set Number of Movement Pulses

use_fr = false;         %Use frequency instead of number of pulses for Physiological Input
if use_fr == true
    fr = 0:0.1:0.3;
    sig = 0.6;
end

num_pulses = 1:1:str2;
num_trials_ident = str3;
num_trials_val = str4;
variable_time = false;

if use_fr == true
    num_models_for_ident = length(fr);
else
    num_models_for_ident = length(num_pulses);
end

validation_frequency = [];
validation_pulses_total = [];

%% Generate the Desired Displacement and Amplitude Modulation Signals for Model Identification

for trial = 1:num_trials_ident
    
    for model = 1:num_models_for_ident

        if variable_time == true 
            t_total = 0:0.001:((physiological_stimulus_time/num_models_for_ident)*model);
            time = ((physiological_stimulus_time/num_models_for_ident)*model);
        else
            t_total = 0:0.001:physiological_stimulus_time;
            time = physiological_stimulus_time;
        end

        if use_fr == true
            FR = makedist('Normal','mu',fr(model),'sigma',sig);
            FrequenciesRandom_max = 2.1;
            FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
            freq_distribution = random(FrequenciesRandom,10000,1);
        end

        AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);

        desired_displacement= 0;
        Freq_test = [];
        Pulses_per_interval_test = [];

        if use_fr == true

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

        else
            t  = 0 : 0.001 : time;

            A = physiological_stimulus_max_amplitude;
            g = time/num_pulses(model);

            DS = makedist('Uniform','lower',0,'upper',g);
            DelayStart = DS;
            delay_start_distribution = random(DelayStart,10000,1);

            d = random(DelayStart,1,1);

            Delay = d;

            for hh = 1:model-1

                %normal distribution from 0 - g
                DR = makedist('Uniform','lower',0,'upper',g);
                DelaysRandom = DR;
                delay_random_distribution = random(DelaysRandom,10000,1);

                random_delay = random(DelaysRandom,1,1);

                Delay = [Delay hh*g+random_delay];
                %Delay = [Delay Delay(end)+random_delay];

            end

            %Delay = (d:g:time)';

            data = (A*pulstran(t,Delay,@rectpuls,W))';
            desired_displacement = data;

            Pulses_per_interval_total = num_pulses(model);
            Freq_test_average = [];

        end

        desired_displacement = desired_displacement';

        stim_amplitude = desired_displacement*170;  %mV
        stim_frequency = 50;
        
        %Theoretical Electrical Stimulus
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

        %% Execute the SRS Simulation
        
        %Set Output noise as zero
        set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('SRS_simulation',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;

        %% Get Output Signals from SRS Simulation (Simulink Model)

        %Muscle Force
        force_simulink = out.SRS_Simulation_Force;

        %Input Stimulus and Output Paralyzed Displacement
        input_stimulus = out.SRS_Simulation_Stimulus;
        output_displacement_simulink = out.SRS_Simulation_Displacement;
        t_simulink = out.tout;

        %% Input/Output for Model Identification
        
        Zcur = [amplitude_modulation',output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

        %% Calculate Signal to Noise Ratio
        
        %Get Output Noise from SRS Simulation
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];

        %% Model Identification (LNL, Hammerstein, Wiener, Linear IRF)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            LNL=lnlbl;
            set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
            'hkTolerance', 0.1, 'nhkMaxIts', 10, 'nhkMaxInner', 5);

            LNL=nlident(LNL,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(LNL,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            SRS_models_temp = [SRS_models_temp LNL];
            LNL = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];
        
        elseif Hammerstein_model == true
        
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Hammerstein=nlbl;
            set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
            
            I=Hammerstein{1,2};
            set(I,'nLags',2400,'nSides',1); %(accuracy increases if nlags is increased)
            Hammerstein{1,2}=I;

            Hammerstein=nlident(Hammerstein,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(Hammerstein,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            SRS_models_temp = [SRS_models_temp Hammerstein];
            Hammerstein = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        elseif Weiner_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Weiner = lnbl; %Wiener
            set(Weiner,'idMethod','hk');
            
            I = Weiner{1,1};
            set(I,'nLags',850,'nSides',1); % Set Number of lags and Sides in IRF
            Weiner{1,1} = I;

            Weiner=nlident(Weiner,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(Weiner,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end

            SRS_models_temp = [SRS_models_temp Weiner];
            Weiner = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        elseif Linear_IRF_model == true
            
            %Identify a two-sided IRF Model
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);

            figure(model)
            [R, V, yp] = nlid_resid(IRF_model,Zcur);

            identification_accuracy(trial,model) = V;

            if use_fr == true
                validation_frequency(trial,model) = Freq_test_average;
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            else
                validation_pulses_total(trial,model) = Pulses_per_interval_total;
            end
            
            SRS_models_temp{model+((trial-1)*num_models_for_ident)} = IRF_model;

            IRF_model = [];
            Zcur_all_temp = [Zcur_all_temp Zcur];

        end

    end
    
    if Linear_IRF_model == true
        SRS_models = SRS_models_temp(:,:);
    else
        SRS_models = [SRS_models;SRS_models_temp];
        SRS_models_temp = [];
    end

    Zcur_all = [Zcur_all;Zcur_all_temp];
    Zcur_all_temp = [];

end

%% Set Initial Parameters for Model Validation

%Noise Parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                                       %Width of signal pulse (seconds)
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%% Generate Desired Displacement and Amplitude Modulation Signals for Model Validation

for trial = 1:num_trials_val

    t_total = 0:0.001:physiological_stimulus_time;
    time = physiological_stimulus_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 2.21;
    FrequenciesRandom = truncate(FR,0,1.8);
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
    
    %Set Output Noise as Zero
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('SRS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %% Get Output Signals from SRS Simulation (Simulink Model)
    
    %Muscle Force
    force_simulink = out.SRS_Simulation_Force;

    %Input Stimulus and Output Paralyzed Displacement
    input_stimulus = out.SRS_Simulation_Stimulus;
    output_displacement_simulink = out.SRS_Simulation_Displacement;
    t_simulink = out.tout;

    %% Input/Output for Model Validation
    
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});
    
    %% Calculate Signal to Noise Ratio
    
    %Get Output Noise from SRS Simulation
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];
    
    %% Model Validation (LNL, Hammerstein, Wiener, or Linear IRF)

    if LNL_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(SRS_models(1,model),Zcur);
            validation_accuracy(trial,model) = V;
        end
        
    elseif Hammerstein_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(SRS_models(1,model),Zcur);
            validation_accuracy(trial,model) = V;
        end
        
    elseif Weiner_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(SRS_models(1,model),Zcur);
            validation_accuracy(trial,model) = V;
        end
        
    elseif Linear_IRF_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models_for_ident
            
            figure(model)
            [R, V, yp] = nlid_resid(SRS_models{model},Zcur);
            validation_accuracy(trial,model) = V;
        end
    end   
end

%% Plot Accuracy (%VAF) vs Number of Movement Pulses in Identification Signal

figNum = 200;

if use_fr == true
    figure(figNum)
    figNum = figNum+1;
    plot(fr,validation_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 26;
    title('Expected Average Frequency vs Validation Accuracy','Fontsize', 36)
    xlabel('Frequency (Hz)','Fontsize', 32)
    ylabel('Accuracy (%VAF)','Fontsize', 32)
    grid on

    figure(figNum)
    figNum = figNum+1;
    plot(validation_frequency,validation_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 26;
    title('Observed Average Frequency vs Validation Accuracy','Fontsize', 36)
    xlabel('Frequency (Hz)','Fontsize', 32)
    ylabel('Accuracy (%VAF)','Fontsize', 32)
    grid on
    
else
    
    identification_accuracy = max(0, identification_accuracy);
    identification_accuracy_mean = mean(identification_accuracy,1);
    identification_accuracy_var = var(identification_accuracy,0,1);
    identification_accuracy_std = std(identification_accuracy,0,1);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(validation_pulses_total(1,:),min(identification_accuracy_mean+identification_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(validation_pulses_total(1,:),max(identification_accuracy_mean-identification_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([validation_pulses_total(1,:) fliplr(validation_pulses_total(1,:))], [min(identification_accuracy_mean+identification_accuracy_std,100) fliplr(max(identification_accuracy_mean-identification_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(validation_pulses_total(1,:),identification_accuracy_mean,'LineWidth', 3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Total Pulses vs Identifcation Accuracy','Fontsize', 24) 
    xlabel('Number of Movement Pulses in Identification Signal','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
    validation_accuracy = max(0, validation_accuracy);
    validation_accuracy_mean = mean(validation_accuracy,1);
    validation_accuracy_var = var(validation_accuracy,0,1);
    validation_accuracy_std = std(validation_accuracy,0,1);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(validation_pulses_total(1,:),min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(validation_pulses_total(1,:),max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([validation_pulses_total(1,:) fliplr(validation_pulses_total(1,:))], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(validation_pulses_total(1,:),validation_accuracy_mean,'LineWidth',3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Total Pulses vs Validation Accuracy','Fontsize', 24) 
    xlabel('Number of Movement Pulses in Identification Signal','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
end

tEnd = toc(tStart)/60


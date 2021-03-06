%% SRS Model Test for Identification Signal Amplitude and Record Length

%Runs each test seperately
%Can Evaluate the Identified Models based on the:
% 1. Amplitude of the PRBS Input used for Identification
% 2. Record Length of PRBS Input used for Identification

%Identifies models with one set of PRBS signals, each with a different max 
%amplitude or record length (depending on test selected), and then validates
%with multiple realizations of the physiological signals (Since PRBS is 
%deterministic, the identifcation results do not vary much so multiple 
%identification trials to get a identification VAF mean and std are unnecessary)

%When running the script, you need to provide the following input:
% 1. Type of Model Test?
%       Choose which model test you want to run, amplitude vs model accuracy
%       or record length vs model accuracy
% 2. Type of Model Structure? LNL/Hammerstein/Wiener/IRF
%       The Model structure used for the Identified SRS (Default is Wiener)
% 3. Number of Validation Trials?
%       Number of validation trials used to calculate the VAF mean and std

clc
clear all

%% User Input Prompts

prompt1 = 'Type of Model Test? Amplitude(Amp)/Record Length(Rec) [Amp]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'Amp') & ~strcmp(str1,'Rec') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Amp';
end

prompt2 = 'Type of Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Wiener]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'LNL') & ~strcmp(str2,'Hamm') & ~strcmp(str2,'Wiener') & ~strcmp(str2,'IRF') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Wiener';
end

prompt3 = 'Number of Validation Trials? 1-50 [30]: ';
str3 = input(prompt3);
if str3<1 | str3>50
    disp('Invalid Input')
    return
elseif isempty(str3)
    str3 = 30;
end

tStart = tic;

%Can change record length or amplitude of signals used for identification
if strcmp(str1,'Amp')
    variable_time = false;
    variable_signal = true;
elseif strcmp(str1,'Rec')
    variable_time = true;
    variable_signal = false;
end

%% Set initial Parameters

%Noise Parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];

%PRBS Signal Parameters (Dependent on the test selected)
if variable_time == true
    %PRBS_stimulus_time = [1:1:15 20:5:25 30:15:90];
    PRBS_stimulus_time = [1 2 4 5 7 9:1:15 20:5:25 30:15:90];
else
    PRBS_stimulus_time = 180;
end
variable_amplitude = false;
N = PRBS_stimulus_time/10;
M = 10000;
if variable_signal == true
    PRBS_amplitude = 5.5:0.5:20;
else
    PRBS_amplitude = 20;
end

%Model Type
LNL_model = false;
Hammerstein_model = false;
Weiner_model = false;
Linear_IRF_model = false;

if strcmp(str2,'LNL')
    LNL_model = true;
elseif strcmp(str2,'Hamm')
    Hammerstein_model = true;
elseif strcmp(str2,'Wiener')
    Weiner_model = true;
elseif strcmp(str2,'IRF')
    Linear_IRF_model = true;
end

%Initialize Models
LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];
IRF_model_all = [];

Zcur_ident = [];

identification_accuracy = [];
validation_accuracy = [];

%% Set number of models/trials
num_record_lengths = length(PRBS_stimulus_time);
num_models_for_ident = length(PRBS_amplitude);

num_trials = str3;

%% Generate Desired Displacement and Amplitude Modulation for Model Identification

for record_length = 1:num_record_lengths
    
    for model = 1:num_models_for_ident
        
        t_total = 0:0.001:PRBS_stimulus_time(record_length);
        time = PRBS_stimulus_time(record_length);

        A = 0;                                             %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude(model);             %First interval is at max amplitude
                else
                    R = rand(1,1)*PRBS_amplitude(model);   %Randomly generate a number between 0 and PRBS Amplitude
                end

                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude(model);                      %Else set as Constant Amplitude
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
        
        %Theoretical Electrical Stimulus
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];
        
        %% Execute SRS Simulation (Simulink Model)
        
        %Set Output Noise as zero
        set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('SRS_simulation',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %% Get Output from SRS Simulation (Simulink Model)

        %Muscle Force,Input Stimulus and Output Paralyzed Displacement
        force_simulink = out.SRS_Simulation_Force;
        input_stimulus = out.SRS_Simulation_Stimulus;
        output_displacement_simulink = out.SRS_Simulation_Displacement;
        t_simulink = out.tout;

        %% Input/Output for Model Identification
        
        Zcur = [amplitude_modulation',output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});
        
        %% Model Identification (LNL, Hammerstein, Wiener, or Linear IRF)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            LNL=lnlbl;
            set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
            'hkTolerance', 0.1, 'nhkMaxIts', 10, 'nhkMaxInner', 5);

            LNL=nlident(LNL,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(LNL,Zcur);

            identification_accuracy = [identification_accuracy V];

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

            figure(model)
            [R, V, yp] = nlid_resid(Hammerstein,Zcur);

            identification_accuracy = [identification_accuracy V];

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

            Weiner=nlident(Weiner,Zcur);

            figure(model)
            [R, V, yp] = nlid_resid(Weiner,Zcur);

            identification_accuracy = [identification_accuracy V];

            Weiner_all = [Weiner_all Weiner];
            Weiner = [];
            Zcur_ident = [Zcur_ident Zcur];

        elseif Linear_IRF_model == true
            
            % Identify a two-sided IRF Model
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);

            figure(model)
            [R, V, yp] = nlid_resid(IRF_model,Zcur);

            identification_accuracy = [identification_accuracy V];
            
            if variable_time == true
                IRF_model_all{record_length} = IRF_model;
            elseif variable_signal == true
                IRF_model_all{model} = IRF_model;
            end
            
            IRF_model = [];
            Zcur_ident = [Zcur_ident Zcur];
            
        end  
        
    end
    
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
W = 0.45;                                       %Width of Movement Pulse (seconds)
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

Zcur_val = [];

%% Generate Desired Displacement and Amplitude Modulation Signal for Model Validation

for trial = 1:num_trials

    t_total = 0:0.001:physiological_stimulus_time;
    time = physiological_stimulus_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 2.21;
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
    
    %Set output noise as zero
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('SRS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %% Get Output Signals from Simulink

    %Muscle Force, Input Stimulus and Output Paralyzed Displacement
    force_simulink = out.SRS_Simulation_Force;
    input_stimulus = out.SRS_Simulation_Stimulus;
    output_displacement_simulink = out.SRS_Simulation_Displacement;
    t_simulink = out.tout;

    %% Input/Output for Model Validation
    
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});
    
    %% Calculate Signal to Noise Ratio
    
    %Get Output noise from SRS Simulation
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];
    
    %% Model Validation (LNL, Hammerstein, Wiener, or Linear IRF)
    
    if variable_time == true
        num_models = num_record_lengths;
    elseif variable_signal == true
        num_models = num_models_for_ident;
    end

    if LNL_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models

            figure(2)
            [R, V, yp] = nlid_resid(LNL_all(model),Zcur);
            validation_accuracy(trial,model) = V;

        end
    
    elseif Hammerstein_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models

            figure(2)
            [R, V, yp] = nlid_resid(Hammerstein_all(model),Zcur);
            validation_accuracy(trial,model) = V;

        end
        
    elseif Weiner_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models

            figure(2)
            [R, V, yp] = nlid_resid(Weiner_all(model),Zcur);
            validation_accuracy(trial,model) = V;

        end
        
    elseif Linear_IRF_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        for model = 1:num_models

            figure(2)
            [R, V, yp] = nlid_resid(IRF_model_all{model},Zcur);
            validation_accuracy(trial,model) = V;

        end

    end
    
    Zcur_val = [Zcur_val Zcur];
    
end

%% Plot Accuracy (%VAF) vs Identifcation Signal Amplitude or Record Length
figNum = 300;

if variable_signal == true
    figure(figNum)
    figNum = figNum+1;
    identification_accuracy = max(0, identification_accuracy);
    identification_accuracy_mean = mean(identification_accuracy);
    identification_accuracy_var = var(identification_accuracy);
    identification_accuracy_std = std(identification_accuracy);
    plot(PRBS_amplitude*0.001,identification_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 14;
    title('Amplitude vs Identification Accuracy','Fontsize', 24) 
    xlabel('Identification Signal Amplitude (m)','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
    validation_accuracy = max(validation_accuracy, 0);
    figure(figNum)
    figNum = figNum+1;
    validation_accuracy_mean = mean(validation_accuracy);
    validation_accuracy_var = var(validation_accuracy);
    validation_accuracy_std = std(validation_accuracy);
    hold on
    plot(PRBS_amplitude*0.001,min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(PRBS_amplitude*0.001,max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([PRBS_amplitude*0.001 fliplr(PRBS_amplitude*0.001)], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(PRBS_amplitude*0.001,validation_accuracy_mean,'LineWidth', 3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Amplitude vs Validation Accuracy','Fontsize', 24) 
    xlabel('Identification Signal Amplitude (m)','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
end

if variable_time == true
    
    figure(figNum)
    figNum = figNum+1;
    identification_accuracy = max(0, identification_accuracy);
    identification_accuracy_mean = mean(identification_accuracy);
    identification_accuracy_var = var(identification_accuracy);
    identification_accuracy_std = std(identification_accuracy);
    plot(PRBS_stimulus_time,identification_accuracy,'LineWidth', 3)
    ax = gca;
    ax.FontSize = 14;
    title('Identification Signal Record Length vs Identification Accuracy','Fontsize', 24) 
    xlabel('Record Length (s)','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
    
    validation_accuracy = max(validation_accuracy, 0);
    validation_accuracy_mean = mean(validation_accuracy);
    validation_accuracy_var = var(validation_accuracy);
    validation_accuracy_std = std(validation_accuracy);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(PRBS_stimulus_time,min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(PRBS_stimulus_time,max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([PRBS_stimulus_time fliplr(PRBS_stimulus_time)], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(PRBS_stimulus_time,validation_accuracy_mean,'LineWidth',3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Identification Signal Record Length vs Validation Accuracy','Fontsize', 24) 
    xlabel('Record Length (s)','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
end

tEnd = toc(tStart)/60

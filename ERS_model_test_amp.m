%% ERS Model Test for Identification Signal Amplitude and Record Length

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
% 2. Number of Validation Trials?
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

prompt2 = 'Number of Validation Trials? 1-50 [30]: ';
str2 = input(prompt2);
if str2<1 | str2>50
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 30;
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

%% Set initial parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];

%PRBS Signal Parameters (Dependent on the test selected)
if variable_time == true
    PRBS_movement_time = [1:1:15 20:5:25 30:15:90];
else
    PRBS_movement_time = 180;
end
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
if variable_signal == true
    PRBS_amplitude = 4.5:0.5:9.5;
else
    PRBS_amplitude = 10;
end

NHK_all = [];
Zcur_all = [];

identification_accuracy = [];
validation_accuracy = [];

%% Set number of models/trials
num_record_lengths = length(PRBS_movement_time);
num_models_for_ident = length(PRBS_amplitude);

num_trials = str2;

%% Generate Desired Displacement for Model Identification

for record_length = 1:num_record_lengths
    
    for model = 1:num_models_for_ident
        
        t_total = 0:0.001:PRBS_movement_time(record_length);
        time = PRBS_movement_time(record_length);
        
        A = 0;                                             %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude(model);             %First interval is at max PRBS amplitude
                else
                    R = rand(1,1)*PRBS_amplitude(model);   %Randomly generate a number between 0 and PRBS amplitude
                end
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude(model);                      %Constant Amplitude
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
        
        %% Create Neural Input for ERS Simulation
        
        %Create Frequency and Amplitude Parameters of Neural Input based on
        %Desired Displacement
        Amplitude = desired_displacement*100;    %mV
        Frequency = desired_displacement*14000;  %Hz

        %Generate Neural Command Signal with Frequency and Amplitdue Parameters
        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

        neural_simulink = [t_total' neural]; %Input for the Simulink Model
        
        %% Execute ERS Simulation
        
        %Set Output Noise (as zero)
        set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('ERS_simulation',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %% Get Output Signals from ERS simulation
        
        %EMG, Muscle Force, and Output Displacement
        emg_simulink = out.EMGout;
        force_simulink = out.EMG_Model_Force;
        output_displacement_simulink = out.EMG_Model_Displacement;
        t_simulink = out.tout;
        
        %% Input/Output for Model Identification
        
        Zcur = [emg_simulink,output_displacement_simulink];
        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
        
        %% Calculate Signal to Noise Ratio
        
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];
        
        %% Model Identification
        
        %Hammerstein System Identification (HK Method)
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        NHK=nlbl;
        set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

        % Set Number of lags in IRF
        I=NHK{1,2};
        set(I,'nLags',min(round(length(Zcur)/2.5),2400));
        NHK{1,2}=I;
        
        %Identify Model
        NHK=nlident(NHK,Zcur);

        figure(model);
        [R, V, yp] = nlid_resid(NHK,Zcur);

        identification_accuracy = [identification_accuracy V];
        
        NHK_all = [NHK_all NHK];
        NHK = [];
        Zcur_all = [Zcur_all Zcur];
        
    end
     
end

%% Set Initial Parameters for Model Validation

noise_snr = [];
set_output_noise_power = 0;
output_noise_power = [];

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;    %Physiological Signal max amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds)
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%% Generate Desired Displacement for Model Validation

for trial = 1:num_trials
    
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
    %Desired Displacement
    Amplitude = desired_displacement*100;    %mV
    Frequency = desired_displacement*14000;  %Hz

    %Generate Neural Command Signal with Frequency and Amplitdue Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model
    
    %% Execute ERS Simulation

    %Set Output Noise
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink Model
    out = sim('ERS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;
    
    %% Get Output Signals from Simulink
    
    %EMG, Muscle Force and Healthy Displacement Output
    emg_simulink = out.EMGout;
    force_simulink = out.EMG_Model_Force;
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
    
    num_models = max(num_models_for_ident, num_record_lengths);
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    %Validates the current Physioloigcal Signal on all identified models
    for model = 1:num_models
        
        figure(10000)
        [R, V, yp] = nlid_resid(NHK_all(model),Zcur);
        
        validation_accuracy(trial,model) = V;
        
    end
    
end

%% Plot Accuracy (%VAF)
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
    % errorbar(validation_pulses_total,validation_accuracy_mean,validation_accuracy_std,'CapSize',14,'LineWidth',3)
    hold on
%     plot(PRBS_stimulus_time*0.001,min(validation_accuracy_mean+validation_accuracy_std,100),'--r','LineWidth',2)
%     plot(PRBS_stimulus_time*0.001,max(validation_accuracy_mean-validation_accuracy_std,0),'--r','LineWidth',2)
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
    
%     figure(figNum)
%     figNum = figNum+1;
%     plot(const_amp*0.001,identification_accuracy,'LineWidth', 3)
%     hold on
%     plot(const_amp*0.001,validation_accuracy,'LineWidth', 3)
%     ax = gca;
%     ax.FontSize = 26;
%     hold off
%     title('Amplitude vs Accuracy','Fontsize', 36) 
%     xlabel('Identification Signal Amplitude (m)','Fontsize', 32)
%     ylabel('Accuracy (%VAF)','Fontsize', 32)
%     legend('Identification Accuracy', 'Validation Accuracy','Fontsize', 26,'Location','southeast')
%     grid on
    
end

if variable_time == true
    
    figure(figNum)
    figNum = figNum+1;
    identification_accuracy = max(0, identification_accuracy);
    identification_accuracy_mean = mean(identification_accuracy);
    identification_accuracy_var = var(identification_accuracy);
    identification_accuracy_std = std(identification_accuracy);
    plot(PRBS_movement_time,identification_accuracy,'LineWidth', 3)
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
    
    PRBS_movement_time = PRBS_movement_time(1,[1:4 6:end]);
    validation_accuracy_mean = validation_accuracy_mean(1,[1:4 6:end]);
    validation_accuracy_std = validation_accuracy_std(1,[1:4 6:end]);
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(PRBS_movement_time,min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
    plot(PRBS_movement_time,max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
    patch([PRBS_movement_time fliplr(PRBS_movement_time)], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)
    plot(PRBS_movement_time,validation_accuracy_mean,'LineWidth',3)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('Identification Signal Record Length vs Validation Accuracy','Fontsize', 24) 
    xlabel('Record Length (s)','Fontsize', 18)
    ylabel('Accuracy (%VAF)','Fontsize', 18)
    grid on
end


tEnd = toc(tStart)/60

        
        
        

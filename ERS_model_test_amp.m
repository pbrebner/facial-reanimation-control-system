%EMG Model Test

%Validate the estimated model by using it to predict the movements
%generated in response to an input realization which is different than that
%used to identify the model

%Tests the Model based on the Amplitude of the PRBS Input used for
%Identification

%Identifies models with one set of PRBS signals and then validates with multiple
%different physiiloigcal signals (Since const amp PRBS is deterministic, the
%identifcation results are always the same for the same signal amplitude)

clc
clear all

tStart = tic;

%Can change record length or amplitude of signals used for identification
variable_time = true;
variable_signal = false;

%%
%Set initla parameters
set_output_noise_power = 0;
%noise_multiplier = 10;
noise_snr = [];
output_noise_power = [];
% figNum = 1;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%PRBS Stimulus
%PRBS_movement_time = 180;
PRBS_movement_time = [1:1:15 20:5:25 30:15:90];
variable_amplitude = false;
N = PRBS_movement_time/10;
M = 10000;
PRBS_amplitude = 10; %mm
%PRBS_amplitude = 4.5:0.5:9.5;

NHK_all = [];
Zcur_all = [];
%emg_all = [];

identification_accuracy = [];
validation_accuracy = [];

%%
num_record_lengths = length(PRBS_movement_time);
num_models_for_ident = length(PRBS_amplitude);

num_trials = 30;

%%
for record_length = 1:num_record_lengths
    
    for model = 1:num_models_for_ident
        
        t_total = 0:0.001:PRBS_movement_time(record_length);
        time = PRBS_movement_time(record_length);
        
        A = 0;                    %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude(model);             %First interval is at max amplitude
                else
                    R = rand(1,1)*PRBS_amplitude(model);   %Randomly generate a number between 0 and 10
                end
                for j = 1:M
                    A = [A R];
                end
            end
        else
            A = PRBS_amplitude(model);              %Constant Amplitude
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
        
        %%
        %Create Frequency and Amplitude Parameters based on Analog Signal
        Amplitude = desired_displacement*100;  %mV
        Frequency = desired_displacement*14000;  %Hz

        %%
        %Generate Neural Command Signal with Frequency and Amplitdue Parameters
        neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(t_total,neural)
    %     ax = gca;
    %     ax.FontSize = 26;
    %     xlabel('Time (s)','Fontsize',32)
    %     ylabel('Amplitude (% of MUs)','Fontsize',32);
    %     title('Neural Command Input','Fontsize',36)
    %     grid on

        neural_simulink = [t_total' neural]; %Input for the Simulink Model

        %Power Pectrum of Neural Input
%         [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);
        
        %%
        %Execute Simulink Model with Set Output Noise
        
        %Set Output Noise
        set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('EMG_Model_Simulink',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %%
        %Get Output Signals from Simulink
        emg_simulink = out.EMGout;

        %PDF of EMG
    %     figure(figNum)
    %     figNum = figNum+1;
    %     emg_pdf = emg_simulink;
    %     emg_pdf(emg_pdf==0) = []; %Removing the zero values
    %     histogram(emg_pdf,'Normalization','pdf')
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Volts (V)','Fontsize',18)
    %     ylabel('Density','Fontsize',18)
    %     title('EMG Distribution','Fontsize',24)
    %     grid on

        %Power Pectrum of EMG
%         [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        force_simulink = out.EMG_Model_Force;

        %Power Pectrum of FOrce
%         [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        output_displacement_simulink = out.EMG_Model_Displacement;
        t_simulink = out.tout;

        Zcur = [emg_simulink,output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

%         figure(figNum)
%         figNum = figNum+1;
%         subplot(2,1,1)
%         plot(t_total,emg_simulink)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('EMG (V)','Fontsize',18);
%         title('(a) Output EMG','Fontsize',24)
%         grid on
%         
%         subplot(2,1,2)
%         plot(t_total,output_displacement_simulink)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Displacement (m)','Fontsize',18);
%         title('(b) Output Displacement','Fontsize',24)
%         grid on
        
        %%
        %Calculate Signal to Noise Ratio
        output_noise_simulink = out.Output_Noise;

        signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
        noise_snr = [noise_snr signal_to_noise];
        
        %%
        %Hammerstein System Identification (EMG-Movement Model)
        %HK Method
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        NHK=nlbl;
        set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

        % Set umber of lags in IRF
        I=NHK{1,2};
        set(I,'nLags',min(round(length(Zcur)/2.5),2400)); %(accuracy increases if nlags is increased)
        NHK{1,2}=I;

        NHK=nlident(NHK,Zcur);

        figure(model);
        [R, V, yp] = nlid_resid(NHK,Zcur);

        identification_accuracy = [identification_accuracy V];
        
        NHK_all = [NHK_all NHK];
        NHK = [];
        Zcur_all = [Zcur_all Zcur];
        %emg_all = [emg_all emg_simulink];
        
    end
     
end

%% Model Validation
% Set Initial Parameters
noise_snr = [];
set_output_noise_power = 0;
%noise_multiplier = 10;
output_noise_power = [];

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = physiological_movement_time/10;                  %number of random signal changes
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

%     figure(figNum)
%     figNum = figNum+1;
%     plot(t_total,desired_displacement)
%     xlabel('Time (s)','Fontsize',18)
%     ylabel('Desired Displacement (m)','Fontsize',18)
%     title('Analog Signal of Desired Movement','Fontsize',24)
%     grid on

    %Power Pectrum of Desired Displacement
%     [Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

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
%     [Pxx2,f2] = pwelch(neural,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    %%
    %Execute Simulink Model with Set Output Noise

    %Set Output Noise
    set_param('EMG_Model_Simulink/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

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
%     [Pxx_emg,f_emg] = pwelch(emg_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    force_simulink = out.EMG_Model_Force;
    
    output_displacement_simulink = out.EMG_Model_Displacement;
    t_simulink = out.tout;

    Zcur = [emg_simulink,output_displacement_simulink];

    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});

    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];
    
    %%
    %Model Validation
    num_models = max(num_models_for_ident, num_record_lengths);
    
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    
    for model = 1:num_models
        
        figure(2)
        [R, V, yp] = nlid_resid(NHK_all(model),Zcur);
        
        validation_accuracy(trial,model) = V;
        
    end
    
end

%% Plot Accuracy (%VAF)
figNum = 600;

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

        
        
        

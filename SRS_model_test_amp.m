%Paralyzed Model Test

%Validate the estimated model by using it to predict the movements
%generated in response to an input realization which is different than that
%used to identify the model

%Tests the Model based on the Amplitude of the PRBS Input used for
%Identification

%Identifies models with one set of PRBS signals and then validates with multiple
%different physiological signals (Since const amp PRBS is deterministic, the
%identifcation results are always the same for the same signal amplitude)

clc
clear all

tStart = tic;

%Can change record length or amplitude of signals used for identification
variable_time = true;
variable_signal = false;

%%
%Set initial Parameters
%set_output_noise_power = 1e-10;        %Output Noise Power
set_output_noise_power = 0;
%noise_multiplier = 10;
noise_snr = [];
output_noise_power = [];
% figNum = 1;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

% PRBS_stimulus = true;
% physiological_stimulus = true;

%PRBS Stimulus
%PRBS_stimulus_time = 180;
PRBS_stimulus_time = [1:1:2 4:1:5 7 9:1:15 20:5:25 30:15:90];
variable_amplitude = false;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm
%PRBS_amplitude = 4.5:0.5:9.5;

%Set Physiological Signal Parameters
% physiological_stimulus_time = 180;
% % physiological_stimulus_amplitude = 6.5;
% fr = 0.1;                 %Frequency distribution mean (Hz)
% sig = 0.8;                %Std of Frequency Distribution (Hz)
% W = 0.45;                  
% nf = physiological_stimulus_time/10;                  %number of random signal changes
% t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
% chance_of_zero = false;

LNL_model = false;
Hammerstein_model = false;
Weiner_model = true;
Linear_IRF_model = true;

LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];
IRF_model_all = [];
NLN_all = [];
Zcur_all = [];

identification_accuracy = [];
validation_accuracy = [];

%%
num_record_lengths = length(PRBS_stimulus_time);
num_models_for_ident = length(PRBS_amplitude);

num_trials = 30;

%%
for record_length = 1:num_record_lengths
    
    for model = 1:num_models_for_ident
        
        t_total = 0:0.001:PRBS_stimulus_time(record_length);
        time = PRBS_stimulus_time(record_length);

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


        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
        Pxx1 = Pxx1';

        f1 = (Fs*(0:(L/2))/L)';
    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(f1(1:5000,:),Pxx1(1:5000,:)) 
    %     title('FFT of Desired Displacement')
    %     xlabel('f (Hz)')
    %     ylabel('|Pxx1(f)|')

        %Power Pectrum of Desired Displacement
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stim_frequency = 50;
        stim_amplitude = desired_displacement*170;

        %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

        %FFT of Input Stimulus
        L = length(input_stimulus);

        Y = fft(input_stimulus);
        P2 = abs(Y/L);
        Pxx2 = P2(1:L/2+1);
        Pxx2(2:end-1) = 2*Pxx2(2:end-1);
        Pxx2 = Pxx2';

        f2 = (Fs*(0:(L/2))/L)';
    %     figure(figNum)
    %     figNum = figNum +1;
    %     plot(f2(1:12000,:),Pxx2(1:12000,:)) 
    %     title('FFT of Input Stimulus')
    %     xlabel('f (Hz)')
    %     ylabel('|Pxx2(f)|')

        %Power Pectrum of Neural Input
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stimulus_simulink = [t_total' input_stimulus'];
        
        %%
        %Execute Simulink Model with Set Output Noise

        set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('Paralyzed_Model_Simulink',time);

        %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
        
        %%
        %Get Output Signals from Simulink

        %Muscle Force
        force_simulink = out.Paralyzed_Model_Force;

        %Power Pectrum of Force
%         [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);


        %%
        %Input Stimulus and Output Displacement
        input_stimulus = out.Paralyzed_Model_Stimulus;
        output_displacement_simulink = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        %Zcur = [input_stimulus',output_displacement_simulink];
        Zcur = [amplitude_modulation',output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

        
        %%
        %Model Identification (LNL Model or NLN Model)

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
            Zcur_all = [Zcur_all Zcur];
            
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

            identification_accuracy = [identification_accuracy V];

            Hammerstein_all = [Hammerstein_all Hammerstein];
            Hammerstein = [];
            Zcur_all = [Zcur_all Zcur];

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
            Zcur_all = [Zcur_all Zcur];

        elseif Linear_IRF_model == true
            % Identify a two-sided IRF Model

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);

            figure(model)
            [R, V, yp] = nlid_resid(IRF_model,Zcur);

            identification_accuracy = [identification_accuracy V];

            IRF_model_all = [IRF_model_all IRF_model];
            IRF_model = [];
            Zcur_all = [Zcur_all Zcur];
            
        else
            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
            Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
            %'domainIncr',0.001,'comment','Input Stimulus, Output Displacement','chanNames', {'Stimulus (V)' 'Displacement (m)'});

            set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');

            NLN = nlhw(Zcur_NLN,[2 2 2]);

            [y,fit,ic] = compare(Zcur_NLN,NLN);

            pred = y.y(:,1);
            V = vaf(output_displacement_simulink,pred);
            identification_accuracy = [identification_accuracy V]; 

            R = output_displacement_simulink - pred;

            %Plots just the superimposed part with %VAF in the title
            figure(model);
            plot(t_total,pred);
            hold on
            plot(t_total, output_displacement_simulink)
            ax = gca;
            ax.FontSize = 18;
            hold off
            title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
            xlabel('Time (s)', 'Fontsize', 24)
            ylabel('Displacement (m)', 'Fontsize', 24)
            legend('Predicted', 'Observed', 'Fontsize', 20)
            grid on

            NLN_all = [NLN_all NLN];
            NLN = [];
            Zcur_all = [Zcur_all Zcur];
            
        end  
        
    end
    
end

%% Model Validation
%Set initial Parameters
%set_output_noise_power = 1e-10;        %Output Noise Power
set_output_noise_power = 0;
%noise_multiplier = 10;
noise_snr = [];
output_noise_power = [];

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.8;                %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;                  %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

%%
for trial = 1:num_trials

    t_total = 0:0.001:physiological_stimulus_time;
    time = physiological_stimulus_time;

    FR = makedist('Normal','mu',fr,'sigma',sig);
    FrequenciesRandom_max = 2.21;
    FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
    freq_distribution = random(FrequenciesRandom,10000,1);
%     figure(figNum)
%     figNum = figNum+1;
%     histogram(freq_distribution,100)
%     title('Frequency Distribution')
%     xlabel('Frequency (Hz)')

    AR = makedist('Uniform','lower',0,'upper',physiological_stimulus_max_amplitude);  %Full Amplitude Range
    AmplitudesRandom = AR;
    amp_distribution = random(AmplitudesRandom,10000,1);
    %figure(figNum)
    %figNum = figNum+1;
    %histogram(amp_distribution,100)
    %title('Amplitude Distribution')
    %xlabel('Amplitude (mm)')

    %PWR = makedist('Normal','mu',0.001,'sigma',0.005);
    %PulseWidthRandom = truncate(PWR,0.001,0.010);   % 1ms - 10 ms
    %pw_distribution = random(PulseWidthRandom,10000,1);
    %figure(107)
    %histogram(pw_distribution,100)
    %title('Pulse Width Distribution')

    desired_displacement= 0;
    Freq_test = [];
    Pulses_per_interval_test = [];

    for j = 1 : nf    
        t  = 0 : 0.001 : t_interval;         % Time Samples
        
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

    
    %FFT of Desired Displacement
    L = length(desired_displacement);
    
    Y = fft(desired_displacement);
    P1 = abs(Y/L);
    Pxx1 = P1(1:L/2+1);
    Pxx1(2:end-1) = 2*Pxx1(2:end-1);
%     Pxx1 = Pxx1';
    
    f1 = (Fs*(0:(L/2))/L)';
%     figure(figNum)
%     figNum = figNum+1;
%     plot(f1(1:5000,:),Pxx1(1:5000,:)) 
%     title('FFT of Desired Displacement')
%     xlabel('f (Hz)')
%     ylabel('|Pxx1(f)|')
    
    %Power Pectrum of Desired Displacement
    %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);


    desired_displacement = desired_displacement';

    stim_amplitude = desired_displacement*170;  %mV
    stim_frequency = 50;

    %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
    
    amplitude_modulation = stim_amplitude;
    amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    %FFT of Input Stimulus
    L = length(input_stimulus);
    
    Y = fft(input_stimulus);
    P2 = abs(Y/L);
    Pxx2 = P2(1:L/2+1);
    Pxx2(2:end-1) = 2*Pxx2(2:end-1);
    Pxx2 = Pxx2';
    
    f2 = (Fs*(0:(L/2))/L)';
%     figure(figNum)
%     figNum = figNum +1;
%     plot(f2(1:15000,:),Pxx2(1:15000,:)) 
%     title('FFT of Input Stimulus')
%     xlabel('f (Hz)')
%     ylabel('|Pxx2(f)|')
    
    %Power Pectrum of Stimulus Input
    %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

    stimulus_simulink = [t_total' input_stimulus'];
    
    %%
    %Execute Simulink Model with Set Output Noise

    set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('Paralyzed_Model_Simulink',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %%
    %Get Output Signals from Simulink

    %Muscle Force
    force_simulink = out.Paralyzed_Model_Force;

    %Power Pectrum of Force
%     [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);


    %%
    %Input Stimulus and Output Displacement
    input_stimulus = out.Paralyzed_Model_Stimulus;
    output_displacement_simulink = out.Paralyzed_Model_Displacement;
    t_simulink = out.tout;

    %Zcur = [input_stimulus',output_displacement_simulink];
    Zcur = [amplitude_modulation',output_displacement_simulink];

    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

    
    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    
    %%
    %Model Validation
    num_models = max(num_models_for_ident, num_record_lengths);

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
            [R, V, yp] = nlid_resid(IRF_model_all(model),Zcur);

            validation_accuracy(trial,model) = V;

        end

    else

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
        set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');
        
        for model = 1:num_models

            [y,fit,ic] = compare(Zcur_NLN,NLN_all(model));

            pred = y.y(:,1);
            V = vaf(output_displacement_simulink,pred);
            validation_accuracy(trial,model) = V; 

            R = output_displacement_simulink - pred;

            %Plots just the superimposed part with %VAF in the title
            figure(2);
            plot(t_total,pred);
            hold on
            plot(t_total, output_displacement_simulink)
            ax = gca;
            ax.FontSize = 18;
            hold off
            title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
            xlabel('Time (s)', 'Fontsize', 24)
            ylabel('Displacement (m)', 'Fontsize', 24)
            legend('Predicted', 'Observed', 'Fontsize', 20)
            grid on

        end

    end
    
end

%% Plot Accuracy (%VAF)
figNum = 6000;

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

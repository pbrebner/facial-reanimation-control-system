% Paralyzed Model Experiment (Accuracy vs Ouput Noise)

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

%%
%Set initial Parameters
noise_level_iters = 18;
set_output_noise_initial = 1e-15;        %Output Noise Power
set_output_noise_power = set_output_noise_initial;
noise_multiplier = 5;
set_seed = 23341;

%For Power Spectrums or FFT
Fs = 1000; 
Nfft = 10000;

% PRBS_stimulus = true;
% physiological_stimulus = true;

%PRBS Stimulus
PRBS_stimulus_time = 180;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm

%Set Physiological Signal Parameters
physiological_stimulus_time = 180;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.8;                %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;                  %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

LNL_model = false;
Hammerstein_model = false;
Weiner_model = true;
Linear_IRF_model = true;

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

LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];
IRF_model_all = [];
NLN_all = [];
models_all = [];
Zcur_all = [];

%%
num_trials = 30;
signals = 2;

%%
for signal = 1:signals
    
    if signal == 1
        PRBS_stimulus = false;
    else
        PRBS_stimulus = true;
    end
    
    if PRBS_stimulus == true

        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;

        A = [0];                    %Intialize amplitude
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

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
        Pxx1 = Pxx1';

        f1 = (Fs*(0:(L/2))/L)';
%         figure(figNum)
%         figNum = figNum+1;
%         plot(f1(1:5000,:),Pxx1(1:5000,:)) 
%         title('FFT of Desired Displacement')
%         xlabel('f (Hz)')
%         ylabel('|Pxx1(f)|')

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
%         figure(figNum)
%         figNum = figNum +1;
%         plot(f2(1:12000,:),Pxx2(1:12000,:)) 
%         title('FFT of Input Stimulus')
%         xlabel('f (Hz)')
%         ylabel('|Pxx2(f)|')

        %Power Pectrum of Neural Input
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stimulus_simulink = [t_total' input_stimulus'];

    else

        t_total = 0:0.001:physiological_stimulus_time;
        time = physiological_stimulus_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 2.1;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);
%         figure(figNum)
%         figNum = figNum+1;
%         histogram(freq_distribution,100)
%         title('Frequency Distribution')
%         xlabel('Frequency (Hz)')

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

%         figure(signal)
%         plot(t_total,desired_displacement)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Desired Displacement (m)','Fontsize',18)
%         title('Analog Signal of Desired Displacement','Fontsize',24)
%         grid on

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
    %     Pxx1 = Pxx1';

        f1 = (Fs*(0:(L/2))/L)';
%         figure(figNum)
%         figNum = figNum+1;
%         plot(f1(1:5000,:),Pxx1(1:5000,:)) 
%         title('FFT of Desired Displacement')
%         xlabel('f (Hz)')
%         ylabel('|Pxx1(f)|')

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
%         figure(figNum)
%         figNum = figNum +1;
%         plot(f2(1:15000,:),Pxx2(1:15000,:)) 
%         title('FFT of Input Stimulus')
%         xlabel('f (Hz)')
%         ylabel('|Pxx2(f)|')

        %Power Pectrum of Stimulus Input
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stimulus_simulink = [t_total' input_stimulus'];

    end
    
    %% Generate Clean Signals for Evaluation
    %Execute Simulink Model with No Noise
    
    set_output_noise_power_clean = 0;
    set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power_clean')
    %output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('Paralyzed_Model_Simulink',time);
    
    %Get Output Signals from Simulink
    %Muscle Force
    force_simulink_clean = out.Paralyzed_Model_Force;
    
    %Output Displacement
    output_displacement_simulink_clean = out.Paralyzed_Model_Displacement;
    t_simulink_clean = out.tout;

    %Zcur = [input_stimulus',output_displacement_simulink];
    Zcur_clean = [amplitude_modulation',output_displacement_simulink_clean];
    Zcur_clean = nldat(Zcur_clean,'domainIncr',0.001,'comment','Input Amplitude Modulation, Clean Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Clean Displacement (m)'});
    
    
    %% Add Noise to the Model
    for noise_level = 1:noise_level_iters
        %Execute Simulink Model with Set Output Noise

        set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power')
        output_noise_power = [output_noise_power set_output_noise_power];

        %Run Simulink;
        out = sim('Paralyzed_Model_Simulink',time);

        set_output_noise_power = noise_multiplier*set_output_noise_power;
        %set_output_noise = set_output_noise_power;
    
        %%
        %Get Output Signals from Simulink

        %Muscle Force
        force_simulink = out.Paralyzed_Model_Force;

        %FFT of Muscle Force
        L = length(force_simulink);

        Y = fft(force_simulink);
        P_force = abs(Y/L);
        Pxx_force = P_force(1:L/2+1);
        Pxx_force(2:end-1) = 2*Pxx_force(2:end-1);
        % Pxx_force = Pxx_force';

        f_force = (Fs*(0:(L/2))/L)';
    %     figure(figNum)
    %     figNum = figNum +1;
    %     plot(f_force(1:50,:),Pxx_force(1:50,:)) 
    %     title('FFT of Muscle Force')
    %     xlabel('f (Hz)')
    %     ylabel('|Pxx_force(f)|')

        %Power Pectrum of Force
        %[Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);

        %%
        %Output Displacement
        output_displacement_simulink = out.Paralyzed_Model_Displacement;
        t_simulink = out.tout;

        %Zcur = [input_stimulus',output_displacement_simulink];
        Zcur = [amplitude_modulation',output_displacement_simulink];

        Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

 
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
        %Model Identification (LNL Model or NLN Model)

        if LNL_model == true

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            LNL=lnlbl;
            set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',2,'nLags2',750,... 
            'hkTolerance', 0.1, 'nhkMaxIts', 10, 'nhkMaxInner', 5);

            LNL=nlident(LNL,Zcur);

            figure(signal)
            [R, V, yp] = nlid_resid(LNL,Zcur);

            accuracy_identification = [accuracy_identification V];

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

            figure(signal)
            [R, V, yp] = nlid_resid(Hammerstein,Zcur);

            accuracy_identification = [accuracy_identification V];

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
            
            %Identify the model with the unclean output
            Weiner=nlident(Weiner,Zcur);
            
            % Evaluate by comparing with the clean Output
            figure(signal)
            %[R, V, yp] = nlid_resid(Weiner,Zcur);
            [R, V, yp] = nlid_resid(Weiner,Zcur_clean);

            accuracy_identification = [accuracy_identification V];

            Weiner_all = [Weiner_all Weiner];
            Weiner = [];
            Zcur_all = [Zcur_all Zcur];

        elseif Linear_IRF_model == true
            % Identify a two-sided IRF Model

            set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

            IRF_model = irf(Zcur,'nLags',1200,'nSides',1);

            figure(signal)
            [R, V, yp] = nlid_resid(IRF_model,Zcur);

            accuracy_identification = [accuracy_identification V];

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
            accuracy_identification = [accuracy_identification V]; 

            R = output_displacement_simulink - pred;

            %Plots just the superimposed part with %VAF in the title
            figure(signal);
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
            Zcur_all = [Zcur_all Zcur];

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
    elseif Linear_IRF_model == true
        models_all = [models_all;IRF_model_all];
    else
        models_all = [models_all;NLN_all];
    end
    
    accuracy_identification = [];
    noise_snr_clean = [];
    noise_snr = [];
    output_noise_power = [];
    LNL_all = [];
    Hammerstein_all = [];
    Weiner_all = [];
    IRF_model_all = [];
    NLN_all = [];

end

%% Model Validation with Physiological Signals
%Set initial Parameters
set_output_noise_power_validation = 0;

%For Power Spectrums
Fs = 1000; 
Nfft = 1000;

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

    set_param('Paralyzed_Model_Simulink/Output Noise','Cov','set_output_noise_power_validation')
    %output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('Paralyzed_Model_Simulink',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %%
    %Get Output Signals from Simulink

    %Muscle Force
    force_simulink = out.Paralyzed_Model_Force;

    %Power Pectrum of Force
    [Pxx_force,f_force] = pwelch(force_simulink,gausswin(Nfft),Nfft/2,Nfft,Fs);


    %%
    %Output Displacement
    output_displacement_simulink = out.Paralyzed_Model_Displacement;
    t_simulink = out.tout;

    %Zcur = [input_stimulus',output_displacement_simulink];
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    %noise_snr = [noise_snr signal_to_noise];
    
    %%
    %Model Validation
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

                figure(signal)
                [R, V, yp] = nlid_resid(models_all(signal,noise_level),Zcur);

                accuracy_validation(trial,noise_level,signal) = V;
                
            end
            
        end

    else

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
        set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');
        
        for signal = 1:signals
            
            for noise_level = 1:noise_level_iters

                [y,fit,ic] = compare(Zcur_NLN,models_all(signal,noise_level));

                pred = y.y(:,1);
                V = vaf(output_displacement_simulink,pred);
                accuracy_validation(trial,noise_level,signal) = V; 

                R = output_displacement_simulink - pred;

                %Plots just the superimposed part with %VAF in the title
                figure(signal);
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
    
end

%% Plot Accuracy (Identification & Validation) vs Noise
figNum = 800;

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
%% Stimulus Response System (SRS)

%Identifies the SRS model with simulated data from the SRS simulation, and 
%plots the results. Can identify models based on a Sinusoidal (Simple) Input, PRBS input 
%or a "Physiological" input.

%The Sinusoidal (Simple) input is a 10 second signal that includes a series of five basic 
%movement pulses of random amplitude

%The PRBS input is a pseudorandom binary sequence, modified to have a
%variable amplitude that randomly changes every 10 seconds

%The "Physiological input was created to simulate rat movements in an
%experimental trial. In includes trains of square pulses with amplitude and
%frequency that randomly change every 10 seconds

%When running the script, you need to provide the following input:
% 1. Type of Model Structure? LNL/Hammerstein/Wiener/IRF
%       The Model structure used for the Identified SRS (Default is Wiener)
% 2. Compare two models (PRBS & Physiological)? Y/N
%       Whether or not you want to identify models from both PRBS and Physiological Inputs.
%       running this will compare the two identifoed models at the end
%if N,
% 3. Type of Input? Simple/PRBS/Physiological
%       Choose the type of input you want to use for Model Identification

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
model_type = str1;

prompt2 = 'Compare two models (PRBS & Physiological)? Y/N [Y]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'Y') & ~strcmp(str2,'N') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Y';
end

if strcmp(str2,'Y')
    str2 = true;
elseif strcmp(str2,'N')
    str2 = false;
end

if str2 == false
    prompt3 = 'Type of Input? Simple/PRBS/Physiological(Phys) [PRBS]: ';
    str3 = input(prompt3,'s');
    if ~strcmp(str3,'Simple') & ~strcmp(str3,'PRBS') & ~strcmp(str3,'Phys') & ~isempty(str3)
        disp('Invalid Input')
        return
    elseif isempty(str3)
        str3 = 'PRBS';
    end
    
    if strcmp(str3,'Simple')
        sinusoidal_stimulus = true;
        PRBS_stimulus = false;
    elseif strcmp(str3,'PRBS')
        sinusoidal_stimulus = false;
        PRBS_stimulus = true;
    else
        sinusoidal_stimulus = false;
        PRBS_stimulus = false;
    end 
    input_type = str3;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];
figNum = 1;

%For Power Spectrums or FFTs
Fs = 1000;
Nfft = 200000;

%Sinusoidal (Simple) Signal Parameters
desired_displacement_frequency = 0.5;
sinusoidal_stimulus_max_amplitude = 0.02;
sine_stimulus_time = 9.999;
sine_stimulus_frequency = 50;
uniphasic = true;

%PRBS Signal Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;      %PRBS can either be constant amplitude or variable amplitude
N = PRBS_stimulus_time/10;      %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                      %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 20;            %Amplitude of PRBS Signal (mm)

%Set Physiological Signal Parameters
physiological_stimulus_time = 480;
physiological_stimulus_max_amplitude = 0.02;    %"Physiological" Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz)
sig = 0.8;                                      %Std of Frequency Distribution (Hz)
W = 0.45;                                       %Width of signal pulse (seconds) 
nf = physiological_stimulus_time/10;            %Number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%Initialize models
SRS_models = [];

Zcur_all = [];
accuracy = [];

%Compare Models from both PRBS and Physiological Input
compare_two_models = str2;

if compare_two_models == true
    sinusoidal_stimulus = false;
    PRBS_stimulus = [true false];
end

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

%% Generate Desired Displacement & Amplitude Modulation Signal for Model Identification
%(Sinusoidal, PRBS, or "Physiological") 

for num_signals = 1:length(PRBS_stimulus)

    if sinusoidal_stimulus == true
        
        Nfft = 10000;
        t_total = 0:0.001:sine_stimulus_time;
        time = sine_stimulus_time;

        AmplitudesRandom = makedist('Uniform','lower',0,'upper',sinusoidal_stimulus_max_amplitude);
        desired_displacement_amplitude = [];

        for ii = 1:ceil(time)/2
            rand_amp = random(AmplitudesRandom,1,1);
            desired_displacement_amplitude = [desired_displacement_amplitude ones(1,2000)*rand_amp];
        end

        desired_displacement = max(-desired_displacement_amplitude.*square(2*pi*desired_displacement_frequency.*t_total),0);

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
        Pxx1 = Pxx1';
        f1 = (Fs*(0:(L/2))/L)';

        stim_amplitude = desired_displacement*170;
        
        if uniphasic == true
            input_stimulus = max(stim_amplitude.*sin(2*pi*sine_stimulus_frequency.*t_total),0);
        else
            input_stimulus = stim_amplitude.*sin(2*pi*sine_stimulus_frequency.*t_total);
        end
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    elseif PRBS_stimulus(num_signals) == true
        
        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;

        A = [0];                                    %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %Initial interval has max PRBS amplitude
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

        %Generate a nonperiodic PRBS of length 100 samples.
        u = idinput(time*1000+1,'prbs',Band,Range);

        %Create an iddata object from the generated signal. 
        %For this example, specify the sample time as 1 second.
        u = iddata([],u,0.001);

        U = (u.InputData)';
        desired_displacement = A.*U;

        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement);
        ax = gca;
        ax.FontSize = 15;
        xlabel('Time (s)','Fontsize',18)
        ylabel('Displacement (m)','Fontsize',18)
        title('PRBS Desired Displacement','Fontsize',24)
        grid on

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
        Pxx1 = Pxx1';
        f1 = (Fs*(0:(L/2))/L)';

        stim_frequency = 50;
        stim_amplitude = desired_displacement*170;
        
        %Theoretical Electrical Stimulus
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    else
        
        t_total = 0:0.001:physiological_stimulus_time;
        time = physiological_stimulus_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 2.21;
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
%         figure(figNum)
%         figNum = figNum+1;
%         histogram(amp_distribution,100)
%         title('Amplitude Distribution')
%         xlabel('Amplitude (mm)')

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

        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement)
        ax = gca;
        ax.FontSize = 15;
        xlabel('Time (s)','Fontsize',18)
        ylabel('Displacement (m)','Fontsize',18)
        title('"Physiological" Desired Displacement','Fontsize',24)
        grid on

        %FFT of Desired Displacement
        L = length(desired_displacement);

        Y = fft(desired_displacement);
        P1 = abs(Y/L);
        Pxx1 = P1(1:L/2+1);
        Pxx1(2:end-1) = 2*Pxx1(2:end-1);
        f1 = (Fs*(0:(L/2))/L)';

        desired_displacement = desired_displacement';

        stim_amplitude = desired_displacement*170;
        stim_frequency = 50;
        
        %Theoretical Electrical Stimulus
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

    end

    %% Execute SRS Simulation (Simulink Model)
    
    %Set Output noise (set at zero)
    set_param('SRS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('SRS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %% Get Output Signals from SRS Simulation (Simulink)

    %Muscle Force
    force_simulink = out.SRS_Simulation_Force;

    %FFT of Muscle Force
    L = length(force_simulink);

    Y = fft(force_simulink);
    P_force = abs(Y/L);
    Pxx_force = P_force(1:L/2+1);
    Pxx_force(2:end-1) = 2*Pxx_force(2:end-1);
    f_force = (Fs*(0:(L/2))/L)';

    %Input Stimulus
    input_stimulus = out.SRS_Simulation_Stimulus;
    
    %FFT of Input Stimulus
    L = length(input_stimulus);

    Y = fft(input_stimulus);
    P2 = abs(Y/L);
    Pxx2 = P2(1:L/2+1);
    Pxx2(2:end-1) = 2*Pxx2(2:end-1);
    f2 = (Fs*(0:(L/2))/L)';
    
    %Paralyzed Displacement Output
    output_displacement_simulink = out.SRS_Simulation_Displacement;
    t_simulink = out.tout;

    %% Input/Output for Model Identification
    
    Zcur = [amplitude_modulation',output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Input Amplitude Modulation, Output Displacement','chanNames', {'Amplitude Modulation (V)' 'Displacement (m)'});

    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(t_total,amplitude_modulation)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Time (s)','Fontsize',18)
    ylabel('Amplitude','Fontsize',18);
    title('(a) Amplitude Modulation, A(t)','Fontsize',24)
    grid on

    subplot(2,1,2)
    plot(t_total,output_displacement_simulink)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Time (s)','Fontsize',18)
    ylabel('Displacement (m)','Fontsize',18);
    title('(b) Output Paralyzed Displacement, Pos_P(t)','Fontsize',24)
    grid on
    
    %% Plot the signals of the SRS Simulation
    
    if sinusoidal_stimulus == true
        
        %Plot Desired Displacement, FFT of Desired Displacement, Input
        %Stimulus, and FFT of Input Stimulus
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,1)
        plot(t_total,desired_displacement);
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Desired Displacement (m)','Fontsize',15)
        title('(a) Analog Signal of Desired Displacement','Fontsize',15)
        grid on

        subplot(2,2,2)
        plot(f1(1:50,1),Pxx1(1:50,1));
        ax = gca;
        ax.FontSize = 12;
        title('(b) FFT of Desired Displacement','Fontsize',15);
        ylabel('Displacement (m)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;

        subplot(2,2,3)
        plot(t_total,input_stimulus)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Amplitude (V)','Fontsize',15);
        title('(c) Electrical Stimulus Input','Fontsize',15)
        grid on

        subplot(2,2,4)
        plot(f2(1:800,1),Pxx2(1:800,1));
        ax = gca;
        ax.FontSize = 12;
        title('(d) FFT of Stimulus Input','Fontsize',15);
        ylabel('Amplitude (V)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;
        
        %Plots Input Stimulus, FFT of Input Stimulus, Muscle Force, and FFT
        %of Muscle Force
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,1)
        plot(t_total,input_stimulus)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Amplitude (V)','Fontsize',15);
        title('(a) Electrical Stimulus Input','Fontsize',15)
        grid on

        subplot(2,2,2)
        plot(f2(1:800,1),Pxx2(1:800,1));
        ax = gca;
        ax.FontSize = 12;
        title('(b) FFT of Stimulus Input','Fontsize',15);
        ylabel('Amplitude (V)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;

        subplot(2,2,3)
        plot(t_total,force_simulink)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Force (N)','Fontsize',15);
        title('(c) Muscle Force','Fontsize',15)
        grid on

        subplot(2,2,4)
        plot(f_force(1:50,1),Pxx_force(1:50,1));
        ax = gca;
        ax.FontSize = 12;
        title('(d) FFT of Muscle Force','Fontsize',15);
        ylabel('Force (N)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;
        
        %Plot Amplitude Modulation, Electrical Stimulus, and Muscle Force
        figure(figNum)
        figNum = figNum+1;
        subplot(3,1,1)
        plot(t_total,amplitude_modulation);
        ax = gca;
        ax.FontSize = 13;
        ylabel('Ampltiude (V)','Fontsize',16)
        title('(a) Amplitude Modulation Input, A(t)','Fontsize',18)
        grid on
        
        subplot(3,1,2)
        plot(t_total,input_stimulus)
        ax = gca;
        ax.FontSize = 13;
        ylabel('Amplitude (V)','Fontsize',16);
        title('(b) Electrical Stimulus','Fontsize',18)
        grid on
        
        subplot(3,1,3)
        plot(t_total,force_simulink)
        ax = gca;
        ax.FontSize = 13;
        xlabel('Time (s)','Fontsize',16)
        ylabel('Force (N)','Fontsize',16);
        title('(c) Muscle Force','Fontsize',18)
        grid on

    else
             
        %Plot Desired Displacement, FFT of Desired Displacement, Electrical 
        %Stimulus, FFT of Electrical Stimulus, Muscle Force and 
        %Paralyzed Displacement Output
        figure(figNum)
        figNum = figNum+1;
        
        if PRBS_stimulus == true
            subplot(3,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 12;
            xlabel('Time (s)','Fontsize',12)
            ylabel('Displacement (m)','Fontsize',12)
            title('(a) PRBS Desired Displacement','Fontsize',14)
            grid on
        else
            subplot(3,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 12;
            xlabel('Time (s)','Fontsize',12)
            ylabel('Displacement (m)','Fontsize',12)
            title('(a) "Physiological" Desired Displacement','Fontsize',14)
            grid on
        end

        subplot(3,2,2)
        plot(f1(1:5000,1),Pxx1(1:5000,1));
        ax = gca;
        ax.FontSize = 12;
        title('(b) FFT of Desired Displacement','Fontsize',14);
        ylabel('Displacement (m)','Fontsize',12); 
        xlabel('Frequency (Hz)','Fontsize',12);
        grid on;

        subplot(3,2,3)
        plot(t_total,input_stimulus)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',12)
        ylabel('Amplitude (V)','Fontsize',12);
        title('(c) Electrical Stimulus Input','Fontsize',14)
        grid on

        subplot(3,2,4)
        plot(f2(1:12000,1),Pxx2(1:12000,1));
        ax = gca;
        ax.FontSize = 12;
        title('(d) FFT of Stimulus Input','Fontsize',14);
        ylabel('Amplitude (V)','Fontsize',12); 
        xlabel('Frequency (Hz)','Fontsize',12);
        grid on;

        subplot(3,2,5)
        plot(t_total,force_simulink)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',12)
        ylabel('Force (N)','Fontsize',12);
        title('(e) Muscle Force','Fontsize',14)
        grid on

        subplot(3,2,6)
        plot(t_total,output_displacement_simulink)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',12)
        ylabel('Displacement (m)','Fontsize',12);
        title('(f) Output Displacement','Fontsize',14)
        grid on
        
        %Plot Desired Displacement, Amplitude Modulation, Muscle Force, and
        %Paralyzed displacement Output
        figure(figNum)
        figNum = figNum+1;

        if PRBS_stimulus == true
            subplot(2,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Time (s)','Fontsize',16)
            ylabel('Displacement (m)','Fontsize',16)
            title('(a) PRBS Desired Displacement','Fontsize',16)
            grid on
        else
            subplot(2,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Time (s)','Fontsize',16)
            ylabel('Displacement (m)','Fontsize',16)
            title('(a) "Physiological" Desired Displacement','Fontsize',16)
            grid on
        end

        subplot(2,2,2)
        plot(t_total,amplitude_modulation)
        ax = gca;
        ax.FontSize = 16;
        xlabel('Time (s)','Fontsize',16)
        ylabel('Amplitude (V)','Fontsize',16);
        title('(b) Amplitude Modulation Input, A(t)','Fontsize',16)
        grid on

        subplot(2,2,3)
        plot(t_total,force_simulink)
        ax = gca;
        ax.FontSize = 16;
        xlabel('Time (s)','Fontsize',16)
        ylabel('Force (N)','Fontsize',16);
        title('(c) Muscle Force','Fontsize',16)
        grid on

        subplot(2,2,4)
        plot(t_total,output_displacement_simulink)
        ax = gca;
        ax.FontSize = 14;
        xlabel('Time (s)','Fontsize',16)
        ylabel('Displacement (m)','Fontsize',16);
        title('(d) Output Paralyzed Displacement, Pos_P(t)','Fontsize',16)
        grid on

    end

    %% Calculate Signal to Noise Ratio
    
    %Get noise from SRS Simulation
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    %% Model Identification (LNL, Hammerstein, Wiener, or Linear IRF)

    if LNL_model == true

        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        LNL=lnlbl;
        set(LNL,'idMethod','hk','nLags1',750,'polyOrderMax',5,'nLags2',750,... 
        'hkTolerance', 0.1, 'nhkMaxIts', 4, 'nhkMaxInner', 4);
        
        LNL=nlident(LNL,Zcur);

        figure(figNum)
        figNum = figNum+1;
        subplot(1,3,1)
        plot(LNL{1,1}); 
        ax = gca;
        ax.FontSize = 13;
        title('(a) Input Linear Element','Fontsize', 16)
        xlabel('Lags (s)','Fontsize',16)
        ylabel('X1','Fontsize',16)
        grid on

        subplot(1,3,2)
        plot(LNL{1,2})
        ax = gca;
        ax.FontSize = 13;
        title('(b) Static Nonlinearity','Fontsize', 16)
        xlabel('Input','Fontsize',16)
        ylabel('Output','Fontsize',16)
        grid on

        subplot(1,3,3)
        plot(LNL{1,3}); 
        ax = gca;
        ax.FontSize = 13;
        title('(c) Output Linear Element','Fontsize', 16)
        xlabel('Lags (s)','Fontsize',16)
        ylabel('X1','Fontsize',16)
        grid on

        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(LNL,Zcur);

        accuracy = [accuracy V];
        
        SRS_models = [SRS_models LNL];
        Zcur_all = [Zcur_all Zcur];
        
    elseif Hammerstein_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        Hammerstein=nlbl;
        set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
        I=Hammerstein{1,2};
        set(I,'nLags',2400,'nSides',1); 
        Hammerstein{1,2}=I;

        Hammerstein=nlident(Hammerstein,Zcur);
        
        figure(figNum)
        figNum = figNum+1;
        subplot(1,2,1)
        plot(Hammerstein{1,1});
        ax = gca;
        ax.FontSize = 15;
        title('(a) Static Nonlinearity','Fontsize', 22)
        xlabel('Amplitude Modulation Input, A(t)','Fontsize',18)
        ylabel('Transformed Amplitude Modulation','Fontsize',18)
        grid on

        subplot(1,2,2)
        plot(Hammerstein{1,2})
        ax = gca;
        ax.FontSize = 15;
        title('Linear Element','Fontsize', 22)
        ylabel('X1', 'Fontsize',18)
        xlabel('Lags (s)','Fontsize',18)
        grid on
        
        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(Hammerstein,Zcur);

        accuracy = [accuracy V];
        
        SRS_models = [SRS_models Hammerstein];
        Zcur_all = [Zcur_all Zcur];
        
    elseif Weiner_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        
        Weiner = lnbl; %Wiener
        set(Weiner,'idMethod','hk');
        I = Weiner{1,1};
        set(I,'nLags',850,'nSides',1); % Set Number of lags and Sides in IRF
        Weiner{1,1} = I;
        
        Weiner=nlident(Weiner,Zcur);
        
        figure(figNum)
        figNum = figNum+1;
        subplot(1,2,1)
        plot(Weiner{1,1});
        ax = gca;
        ax.FontSize = 15;
        title('(a) Linear Element','Fontsize', 22)
        ylabel('X1', 'Fontsize',18)
        xlabel('Lags (s)','Fontsize',18)
        grid on

        subplot(1,2,2)
        plot(Weiner{1,2})
        ax = gca;
        ax.FontSize = 15;
        title('(b) Static Nonlinearity','Fontsize', 22)
        xlabel('Transformed Displacement Input (m)','Fontsize',18)
        ylabel('Paralyzed Displacement Output, Pos_P(t) (m)','Fontsize',18)
        grid on
                
        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(Weiner,Zcur);

        accuracy = [accuracy V];
        
        SRS_models = [SRS_models Weiner];
        Zcur_all = [Zcur_all Zcur];
        
    elseif Linear_IRF_model == true
        
        %Identify a Linear IRF Model
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        IRF_model = irf(Zcur,'nLags',1200,'nSides',1);
        
        figure(figNum)
        figNum = figNum+1;
        plot(IRF_model)
        title('IRF Model','Fontsize', 20)
        
        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(IRF_model,Zcur);

        accuracy = [accuracy V];
        
        if compare_two_models == true && PRBS_stimulus(num_signals) == true
            SRS_models{1} = IRF_model;
        elseif compare_two_models == true && PRBS_stimulus(num_signals) == false
            SRS_models{2} = IRF_model;
        else
            SRS_models{1} = IRF_model;
        end
        
        Zcur_all = [Zcur_all Zcur];

    end
    
    %Plot the Results of the Model Identifcation
    %Plots just the superimposed Displacements with %VAF in the title
    figure(figNum);
    figNum = figNum+1;
    pred = double(yp);
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 18;
    hold off
    title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 26)
    xlabel('Time (s)', 'Fontsize', 22)
    ylabel('Paralyzed Displacement, Pos_P(t) (m)', 'Fontsize', 22)
    legend('Predicted', 'Observed', 'Fontsize', 20)

    %Plots Residuals, Residual Distribution and Residual Spectrum
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(R)
    ax = gca;
    ax.FontSize = 15;
    title('(a) Residuals from SRS Identification','Fontsize',18)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Displacement (m)','Fontsize',20)
    grid on

    subplot(2,2,3)
    histogram(double(R))
    ax = gca;
    ax.FontSize = 15;
    xlabel('Displacement (m)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('(b) Residual Distribution','Fontsize',18)
    grid on

    R_zero = R - mean(R);

    if sinusoidal_stimulus == true

        [PxxR,fR] = pwelch(double(R_zero),Nfft,[],Nfft,Fs);
        subplot(2,2,4)
        semilogy(fR(1:200),PxxR(1:200,:),'LineWidth',1);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on

    else

        [PxxR,fR] = pwelch(double(R_zero),Nfft,[],Nfft,Fs);
        subplot(2,2,4)
        semilogy(fR(1:1300),PxxR(1:1300,:),'LineWidth',1);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
        
    end
    
end

%% Plot the Models (Compares the models Identified from PRBS Input and Physiological Input)

if compare_two_models == true && LNL_model == true
    
    LNL1 = SRS_models(1);
    LNL2 = SRS_models(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,3,1)
    hold on
    plot(LNL1{1,1}) 
    plot(LNL2{1,1})
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(a) Input Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    grid on

    subplot(1,3,2)
    hold on
    plot(LNL1{1,2})
    plot(LNL2{1,2})
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(b) Static Nonlinearity','Fontsize', 16)
    xlabel('Input','Fontsize',16)
    ylabel('Output','Fontsize',16)
    grid on

    subplot(1,3,3)
    hold on
    plot(LNL1{1,3}); 
    plot(LNL2{1,3});
    hold off
    ax = gca;
    ax.FontSize = 13;
    title('(c) Output Linear Element','Fontsize', 16)
    xlabel('Lags (s)','Fontsize',16)
    ylabel('X1','Fontsize',16)
    legend('PRBS', 'Physiological','Fontsize',13)
    grid on
    
elseif compare_two_models == true && Hammerstein_model == true
    
    Hammerstein1 = SRS_models(1);
    Hammerstein2 = SRS_models(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(Hammerstein1{1,1})
    hold on
    plot(Hammerstein{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 24)
    xlabel('Input','Fontsize',18)
    ylabel('Output','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(Hammerstein1{1,2})
    hold on
    plot(Hammerstein2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 24)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize', 14)
    grid on
    hold off
    
elseif compare_two_models == true && Weiner_model == true
        
    Weiner1 = SRS_models(1);
    Weiner2 = SRS_models(2);
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(Weiner1{1,1})
    hold on
    plot(Weiner2{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Linear Elements','Fontsize', 22)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on
    hold off

    subplot(1,2,2)
    plot(Weiner1{1,2})
    hold on
    plot(Weiner2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Static Nonlinearities','Fontsize', 22)
    xlabel('Transformed Displacement Input (m)','Fontsize',18)
    ylabel('Paralyzed Displacement Output, Pos_P(t) (m)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize', 14)
    grid on
    hold off
    
elseif compare_two_models == true && Linear_IRF_model == true
    
    IRF1 = SRS_models{1};
    IRF2 = SRS_models{2};
    
    figure(figNum)
    figNum = figNum+1;
    hold on
    plot(IRF1) 
    plot(IRF2)
    hold off
    ax = gca;
    ax.FontSize = 14;
    title('IRF Models','Fontsize', 20)
    xlabel('Lags (s)','Fontsize',20)
    ylabel('X1','Fontsize',20)
    legend('PRBS', 'Physiological','Fontsize',18)
    grid on
    
end

tEnd = toc(tStart)/60

%% EMG Response System (ERS)

%Identifies the ERS model with simulated data from the ERS simulation, and 
%plots the results. Can identify models based on a Simple Input, PRBS input or a 
%"Physiological" input.

%The Simple input is a 10 second signal that includes a series of four basic 
%movement pulses of random amplitude

%The PRBS input is a pseudorandom binary sequence, modified to have a
%variable amplitude that randomly changes every 10 seconds

%The "Physiological input was created to simulate rat movements in an
%experimental trial. In includes trains of square pulses with amplitude and
%frequency that randomly change every 10 seconds

%When running the script, you need to provide the following input:
% 1. Compare two models (PRBS & Physiological)? Y/N (Default is Y)
%       Whether or not you want to identify models from both PRBS and Physiological Inputs.
%       running this will compare the two identifoed models at the end
%if N,
% 2. Type of Input? Simple/PRBS/Physiological (Default is PRBS)
%       Choose the type of input you want to use for Model Identification

clc
clear all

%% User Input Prompts

prompt1 = 'Compare two models (PRBS & Physiological)? Y/N [Y]: ';
str1 = input(prompt1,'s');
if ~strcmp(str1,'Y') & ~strcmp(str1,'N') & ~isempty(str1)
    disp('Invalid Input')
    return
elseif isempty(str1)
    str1 = 'Y';
end

if strcmp(str1,'Y')
    str1 = true;
elseif strcmp(str1,'N')
    str1 = false;
end

if str1 == false
    prompt2 = 'Type of Input? Simple/PRBS/Physiological(Phys) [PRBS]: ';
    str2 = input(prompt2,'s');
    if ~strcmp(str2,'Simple') & ~strcmp(str2,'PRBS') & ~strcmp(str2,'Phys') & ~isempty(str2)
        disp('Invalid Input')
        return
    elseif isempty(str2)
        str2 = 'PRBS';
    end
    
    if strcmp(str2,'Simple')
        simple_movement = true;
        PRBS_movement = false;
    elseif strcmp(str2,'PRBS')
        simple_movement = false;
        PRBS_movement = true;
    else
        simple_movement = false;
        PRBS_movement = false;
    end      
    
    input_type = str2;
end

tStart = tic;

%% Set initial Parameters

%Noise Parameters
set_output_noise_power = 0;
noise_snr = [];
output_noise_power = [];
figNum = 1;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%PRBS Signal Parameters
PRBS_movement_time = 180;
variable_amplitude = true;     %PRBS can either be constant amplitude or variable amplitude
N = PRBS_movement_time/10;     %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                     %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 10;           %PRBS Amplitude (mm)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;    %"Physiological" Amplitude (m)
fr = 0.1;                                       %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                                      %Std of Frequency Distribution (Hz)
W = 0.55;                                       %Width of signal pulse (seconds) 
nf = physiological_movement_time/10;            %Number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (seconds)
chance_of_zero = false;

%Initialize
accuracy = [];
NHK_all = [];
Zcur_all = [];
emg_all = [];

%Compare Models Identified from PRBS and Physiological Inputs
compare_two_models = str1;

if compare_two_models == true
    simple_movement = false;
    PRBS_movement = [true false];
end

%% Generate Desired Displacement Signal (Simple, PRBS, or "Physiological") for Model Identification

for signal = 1:length(PRBS_movement)

    if simple_movement == true
        time = 10;
        t_total = 0:0.001:time;
        t_total_with_delay = 0:0.001:time-0.5;
        
        w = 0.8*pi;
        phi = 0;
        
        A = ones(1,length(t_total_with_delay));
        
        AR = makedist('Uniform','lower',0,'upper',0.01);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);
        A(1:2000) = random(AmplitudesRandom,1,1);
        A(2001:4500) = random(AmplitudesRandom,1,1);
        A(4501:7000) = random(AmplitudesRandom,1,1);
        A(7001:9501) = random(AmplitudesRandom,1,1);
        
        delay = zeros(1,500);
        
        desired_displacement = [delay A.*square(w*t_total_with_delay+phi)];
        desired_displacement = max(desired_displacement,0);
        
        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement);
        ax = gca;
        ax.FontSize = 20;
        xlabel('Time (s)','Fontsize',26)
        ylabel('Desired Displacement (m)','Fontsize',26)
        title('Analog Signal of Desired Displacement','Fontsize',28)
        grid on
        
        %[Pxx1,f1] = pwelch(desired_displacement,Nfft,[],Nfft,Fs);
        Nfft = Nfft/10;
        [Pxx1a,f1a] = pwelch(desired_displacement(1,499:1750),Nfft,[],Nfft,Fs);
        [Pxx1b,f1b] = pwelch(desired_displacement(1,2999:4250),Nfft,[],Nfft,Fs);
        [Pxx1c,f1c] = pwelch(desired_displacement(1,5499:6750),Nfft,[],Nfft,Fs);
        [Pxx1d,f1d] = pwelch(desired_displacement(1,7999:9250),Nfft,[],Nfft,Fs);
    
    elseif PRBS_movement(signal) == true
        
        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;

        A = 0;                                      %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;             %First interval is at max PRBS_amplitude
                else
                    R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS_amplitude
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

        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement);
        ax = gca;
        ax.FontSize = 16;
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        title('PRBS Desired Displacement','Fontsize',24)
        grid on
        
        %Power Pectrum of Desired Displacement
        desired_displacement_zero = desired_displacement - mean(desired_displacement);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);
              
    else   
        % "Physiological" Movement
        t_total = 0:0.001:physiological_movement_time;
        time = physiological_movement_time;

        FR = makedist('Normal','mu',fr,'sigma',sig);
        FrequenciesRandom_max = 1.8;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);
%         figure(figNum)
%         figNum = figNum+1;
%         histogram(freq_distribution,100)
%         title('Frequency Distribution')
%         xlabel('Frequency (Hz)')

        AR = makedist('Uniform','lower',0,'upper',physiological_movement_max_amplitude);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);
%         figure(figNum)
%         figNum = figNum+1;
%         histogram(amp_distribution,100)
%         title('Amplitude Distribution')
%         xlabel('Amplitude (mm)')

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
        
        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement)
        ax = gca;
        ax.FontSize = 16;
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        title('"Physiological" Desired Displacement','Fontsize',24)
        grid on

        %Power Pectrum of Desired Displacement
        desired_displacement_zero = desired_displacement - mean(desired_displacement);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

        desired_displacement = desired_displacement';

    end

    %% Create Neural Input for ERS Simulation
    
    %Create Frequency and Amplitude Parameters of the Neural Input based on
    %Desired Displacement Amplitude
    Amplitude = desired_displacement*100;    %mV
    Frequency = desired_displacement*14000;  %Hz

    %Generate Neural Command Signal with Frequency and Amplitude Parameters
    neural = (max(Amplitude.*square(2*pi*Frequency.*t_total),0))';

    neural_simulink = [t_total' neural]; %Input for the Simulink Model

    %Power Pectrum of Neural Input
    neural_zero = neural - mean(neural);
    
    if simple_movement == true
        Nfft = 1000;
        [Pxx2a,f2a] = pwelch(neural_zero(499:1750,1),Nfft,[],Nfft,Fs);
        [Pxx2b,f2b] = pwelch(neural_zero(2999:4250,1),Nfft,[],Nfft,Fs);
        [Pxx2c,f2c] = pwelch(neural_zero(5499:6750,1),Nfft,[],Nfft,Fs);
        [Pxx2d,f2d] = pwelch(neural_zero(7999:9250,1),Nfft,[],Nfft,Fs);
    else
        [Pxx2,f2] = pwelch(neural_zero,Nfft,[],Nfft,Fs);
    end

    %% Execute ERS Simulation (Simulink Model)

    %Set Output Noise as zero
    set_param('ERS_simulation/Output Noise','Cov','set_output_noise_power')
    output_noise_power = [output_noise_power set_output_noise_power];

    %Run Simulink;
    out = sim('ERS_simulation',time);

    %set_output_noise_power = ii*noise_multiplier*set_output_noise_power;
    %set_output_noise = set_output_noise_power;

    %% Get Output Signals from ERS Simulation Model
    
    %EMG
    emg_simulink = out.ERS_Simulation_EMG;
    
    %PDF of EMG
    figure(figNum)
    figNum = figNum+1;
    emg_pdf = emg_simulink;
    emg_pdf(emg_pdf==0) = []; %Removing the zero values
    histogram(emg_pdf,'Normalization','pdf')
    ax = gca;
    ax.FontSize = 15;
    xlabel('Volts (V)','Fontsize',18)
    ylabel('Density','Fontsize',18)
    title('EMG Distribution','Fontsize',24)
    grid on
    
    %Muscle Force
    force_simulink = out.ERS_Simulation_Force;
    
    %Power Spectrum of Muscle Force
    force_simulink_zero = force_simulink - mean(force_simulink);

    if simple_movement == true
        Nfft = 1000;
        [Pxx_force1,f_force1] = pwelch(force_simulink_zero(499:1750,1),Nfft,[],Nfft,Fs);
        [Pxx_force2,f_force2] = pwelch(force_simulink_zero(2999:4250,1),Nfft,[],Nfft,Fs);
        [Pxx_force3,f_force3] = pwelch(force_simulink_zero(5499:6750,1),Nfft,[],Nfft,Fs);
        [Pxx_force4,f_force4] = pwelch(force_simulink_zero(7999:9250,1),Nfft,[],Nfft,Fs);
    else
        [Pxx_force,f_force] = pwelch(force_simulink_zero,Nfft,[],Nfft,Fs);
    end
    
    %Healthy Output Displacement and Time
    output_displacement_simulink = out.ERS_Simulation_Displacement;
    t_simulink = out.tout;
    
    %% Plot ERS Simulation Signals

    if simple_movement == true
        
        %Plot the Desired Displacement, Desired Displacement Spectrum,
        %Neural Command, and Neural Command Input for Simple Movement
        figure(figNum)
        figNum = figNum+1;
        
        subplot(2,2,1)
        plot(t_total,desired_displacement,'Linewidth',1.5);
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Displacement (m)','Fontsize',15)
        title('(a) Desired Displacement','Fontsize',14)
        grid on
        
        subplot(2,2,2)
        semilogy(f1a(1:11,1),Pxx1a(1:11,1),f1b(1:11,1),Pxx1b(1:11,1),f1c(1:11,1),...
            Pxx1c(1:11,1),f1d(1:11,1),Pxx1d(1:11,1),'Linewidth',1);
        ax = gca;
        ax.FontSize = 12;
        title('(b) Power Spectrum of Displacement','Fontsize',14);
        ylabel('PSD (log)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;
        
        subplot(2,2,3)
        plot(t_total,neural)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Amplitude (% of MUs)','Fontsize',15);
        title('(c) Neural Command Input','Fontsize',14)
        grid on
        
        subplot(2,2,4)
        hold on
        semilogy(f2a(1:151,1),Pxx2a(1:151,1),'Linewidth',1);
        semilogy(f2b(1:151,1),Pxx2b(1:151,1),'Linewidth',1);
        semilogy(f2c(1:151,1),Pxx2c(1:151,1),'Linewidth',1);
        semilogy(f2d(1:151,1),Pxx2d(1:151,1),'Linewidth',1);
        hold off
        ax = gca;
        ax.FontSize = 12;
        title('(d) Power Spectrum of Neural Input','Fontsize',14);
        ylabel('PSD','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        legend('1st Movement','2nd Movement','3rd Movement','4th Movement','Fontsize',12)
        grid on;
        
        % Plot the Neural Input, Neural Input Spectrum, Muscle Force, and
        % Muscle Force Spectrum for a Simple Movement
        figure(figNum)
        figNum = figNum+1;
        
        subplot(2,2,1)
        plot(t_total,neural)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Amplitude (% of MUs)','Fontsize',15);
        title('(a) Neural Input Command','Fontsize',14)
        grid on
        
        subplot(2,2,2)
        hold on
        semilogy(f2a(1:151,1),Pxx2a(1:151,1),'Linewidth',1);
        semilogy(f2b(1:151,1),Pxx2b(1:151,1),'Linewidth',1);
        semilogy(f2c(1:151,1),Pxx2c(1:151,1),'Linewidth',1);
        semilogy(f2d(1:151,1),Pxx2d(1:151,1),'Linewidth',1);
        hold off
        ax = gca;
        ax.FontSize = 12;
        title('(b) Power Spectrum of Neural Input','Fontsize',14);
        ylabel('PSD','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        grid on;
        
        subplot(2,2,3)
        plot(t_total,force_simulink,'Linewidth',1.5)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',15)
        ylabel('Force (N)','Fontsize',15);
        title('(c) Muscle Force','Fontsize',14)
        grid on
        
        subplot(2,2,4)
        semilogy(f_force1(1:11,1),Pxx_force1(1:11,1),f_force2(1:11,1),Pxx_force2(1:11,1),...
            f_force3(1:11,1),Pxx_force3(1:11,1),f_force4(1:11,1),Pxx_force4(1:11,1),'Linewidth',1);
        ax = gca;
        ax.FontSize = 12;
        title('(d) Power Spectrum of Muscle Force','Fontsize',14);
        ylabel('PSD (log)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        legend('1st Movement','2nd Movement','3rd Movement','4th Movement','Fontsize',12)
        grid on;
        
    else

        %Plot Desired Displacement, Neural Command Input and Muscle Force on one figure
        figure(figNum)
        figNum = figNum+1;

        if PRBS_movement(1,signal) == true

            subplot(3,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 12;
            xlabel('Time (s)','Fontsize',12)
            ylabel('Displacement (m)','Fontsize',12)
            title('(a) PRBS Desired Displacement','Fontsize',14)
            grid on

            subplot(3,2,2)
            semilogy(f1(1:301,1),Pxx1(1:301,1),'LineWidth',1.5);
            ax = gca;
            ax.FontSize = 12;
            title('(b) Power Spectrum of Desired Displacement','Fontsize',14);
            ylabel('PSD (log)','Fontsize',12); 
            xlabel('Frequency (Hz)','Fontsize',12);
            grid on;

        else

            subplot(3,2,1)
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 12;
            xlabel('Time (s)','Fontsize',12)
            ylabel('Displacement (m)','Fontsize',12)
            title('(a) "Physiological" Desired Displacement','Fontsize',14)
            grid on

            subplot(3,2,2)
            semilogy(f1(1:60,1),Pxx1(1:60,1),'LineWidth',1.5);
            ax = gca;
            ax.FontSize = 12;
            title('(b) Power Spectrum of Desired Displacement','Fontsize',14);
            ylabel('PSD (log)','Fontsize',12); 
            xlabel('Frequency (Hz)','Fontsize',12);
            grid on;

        end

        subplot(3,2,3)
        plot(t_total,neural)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',12)
        ylabel('Amplitude (% of MUs)','Fontsize',12);
        title('(c) Neural Command Input','Fontsize',14)
        grid on

        subplot(3,2,4)
        semilogy(f2(1:1501,1),Pxx2(1:1501,1));
        ax = gca;
        ax.FontSize = 12;
        title('(d) Power Spectrum of Neural Input','Fontsize',14);
        ylabel('PSD (log)','Fontsize',12); 
        xlabel('Frequency (Hz)','Fontsize',12);
        grid on;

        subplot(3,2,[5 6])
        plot(t_total,force_simulink)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (s)','Fontsize',12)
        ylabel('Force (N)','Fontsize',12);
        title('(e) Muscle Force','Fontsize',14)
        grid on
    end

    %% Input/Output for Model Identification
    
    %EMG input and Healthy Displacement Output
    Zcur = [emg_simulink,output_displacement_simulink];
    Zcur = nldat(Zcur,'domainIncr',0.001,'comment','Output EMG, Output Displacement','chanNames', {'EMG (V)' 'Displacement (m)'});
    
    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(t_total,emg_simulink)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Time (s)','Fontsize',18)
    ylabel('EMG (V)','Fontsize',18);
    title('(a) Output EMG, E(t)','Fontsize',24)
    grid on
    
    subplot(2,1,2)
    plot(t_total,output_displacement_simulink)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Time (s)','Fontsize',18)
    ylabel('Displacement (m)','Fontsize',18);
    title('(b) Output Healthy Displacement, Pos_H(t)','Fontsize',24)
    grid on
    
    %% Calculate Signal to Noise Ratio
    
    %Get noise from ERS simulation
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    %% ERS Model Identification
    
    %Hammerstein System Identification (HK Method)
    set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
    NHK=nlbl;
    set(NHK,'idMethod','hk','displayFlag',true,'threshNSE',.001);

    %Set number of lags in IRF
    I=NHK{1,2};
    set(I,'nLags',2400);
    NHK{1,2}=I;
    
    %Identify the Hammerstein Model and Normalize the Linear Element
    NHK=nlident(NHK,Zcur);
    NHK = normCoefLE(NHK);
    
    %Set the location of the red dashed lines in plots below
    upper_limit = 0.00067;
    lower_limit = -0.00067;

    %Plots NL Hammerstein Part and EMG Distribution
    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(NHK{1,1});
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xline(upper_limit,'--r','LineWidth',1.5);
    xline(lower_limit,'--r','LineWidth',1.5);
    grid on

    subplot(2,1,2)
    histogram(emg_pdf,'Normalization','pdf')
    ax = gca;
    ax.FontSize = 15;
    xline(upper_limit,'--r','LineWidth',1.5);
    xline(lower_limit,'--r','LineWidth',1.5);
    xlabel('EMG Amplitude (V)','Fontsize',18)
    ylabel('Density','Fontsize',18)
    title('(b) EMG Distribution','Fontsize',24)
    grid on
    
    %Plot NL and Linear of Hammerstein System with only High Quality Values
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK{1,1});
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xlim([lower_limit upper_limit])
    grid on

    subplot(1,2,2)
    plot(NHK{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Element','Fontsize', 24)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on

    %Determine the accuracy and residuals of the model identification
    %Plots the Observed, Predicted, Superimposed, and Residuals
    figure(figNum);
    figNum = figNum+1;
    [R, V, yp] = nlid_resid(NHK,Zcur);
    
    accuracy = [accuracy V];

    %Plots just the superimposed part with %VAF in the title
    figure(figNum);
    figNum = figNum+1;
    pred = double(yp);
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 18;
    hold off
    title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 28)
    xlabel('Time (s)', 'Fontsize', 24)
    ylabel('Healthy Displacement, Pos_H(t) (m)', 'Fontsize', 24)
    legend('Predicted', 'Observed', 'Fontsize', 20)
    grid on

    %Plots the Residuals, Residuals Distribution, and Residual Power Spectrum
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(R)
    ax = gca;
    ax.FontSize = 15;
    title('(a) Residuals of ERS Hammerstein Model','Fontsize',20)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Displacement (m)','Fontsize',20)
    grid on

    subplot(2,2,3)
    p = pdf(R);
    plot(p)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Displacement (m)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('(b) Residual Distribution','Fontsize',20)
    grid on
    
    R_zero = R - mean(R);
    R_double = double(R_zero);
     
    [PxxR,fR] = pwelch(R_double,[],[],[],Fs);
    subplot(2,2,4)
    semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
    ax = gca;
    ax.FontSize = 15;
    title('(c) Power Spectrum of Residuals','Fontsize',20);
    ylabel('PSD (log)','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on

    %Plots the Residuals, Residuals Distribution (With Normal Distribution), and Residuals Power
    %Spectrum
    if PRBS_movement == true
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,[1 2])
        plot(R)
        ax = gca;
        ax.FontSize = 15;
        title('(a) Residuals of ERS Hammerstein Model','Fontsize',20)
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        grid on

        subplot(2,2,3)
        p = pdf(R);
        plot(p)
        hold on
        pd = fitdist(R_double,'Normal');
        x_values = -0.002:0.00004:0.002;
        z = pdf(pd,x_values);
        plot(x_values,z,'r', 'Linewidth',2.5)
        hold off
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',20)
        legend('Observed','Theoretical','Fontsize',15)
        grid on

        subplot(2,2,4)
        semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Power Spectrum of Residuals','Fontsize',20);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
    end
    
    %Plots the Residuals, Residuals Distribution, and Autocorrelation
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(R)
    ax = gca;
    ax.FontSize = 15;
    title('(a) Residuals of ERS Hammerstein Model','Fontsize',20)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Displacement (m)','Fontsize',20)
    grid on

    subplot(2,2,3)
    p = pdf(R);
    plot(p)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Displacement (m)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('(b) Residual Distribution','Fontsize',20)
    grid on

    L = length(R_double)-1;
    maxlag = L*0.10;
    [res_corr,lags] = xcorr(R_double,maxlag,'normalized');
    lags = lags/Fs;
    subplot(2,2,4),plot(lags,res_corr,'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 15;
    title('(c) Auto-Correlation of Residuals','Fontsize',20)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Correlation','Fontsize',20)
    grid on

    %Plots the Superimposed, Residuals Distribution, and Residuals Power
    %Spectrum
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,2,[1 2])
%     plot(t_total,pred);
%     hold on
%     plot(t_total, output_displacement_simulink)
%     ax = gca;
%     ax.FontSize = 16;
%     hold off
%     title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 20)
%     xlabel('Time (s)', 'Fontsize', 20)
%     ylabel('Healthy Displacement, Pos_H(t) (m)', 'Fontsize', 20)
%     legend('Predicted', 'Observed', 'Fontsize', 18)
% 
%     subplot(2,2,3)
%     p = pdf(R);
%     plot(p)
%     ax = gca;
%     ax.FontSize = 16;
%     xlabel('Displacement (m)','Fontsize',20)
%     ylabel('Density','Fontsize',20)
%     title('(b) Residual Distribution','Fontsize',18)
%     grid on
%     
%     R_zero = R - mean(R);
%     S = spect(R_zero);
%     subplot(2,2,4)
%     plot(S(1:150,:));
%     ax = gca;
%     ax.FontSize = 16;
%     title('(c) Power Spectrum of Residuals','Fontsize',18);
%     ylabel('PSD','Fontsize',20); 
%     xlabel('Frequency (Hz)','Fontsize',20);
%     grid on
    
    NHK_all = [NHK_all NHK];
    Zcur_all = [Zcur_all Zcur];
    emg_all = [emg_all emg_simulink];
    
end

%% Compare the Models Identifified from Different Input Types

if compare_two_models == true
    NHK1 = NHK_all(1);
    NHK2 = NHK_all(2);
    
    emg_pdf1 = emg_all(:,1);
    emg_pdf2 = emg_all(:,2);
    emg_pdf1(emg_pdf1==0) = []; %Removing the zero values
    emg_pdf2(emg_pdf2==0) = []; %Removing the zero values
    
    %Sets the Upper Limit and Lower limit for "High Quality" values of the
    %model nonlinear element
    upper_limit1 = 0.000586;
    lower_limit1 = -0.000613;
    upper_limit2 = 0.000602;
    lower_limit2 = -0.00063;
    upper_limit3 = max(upper_limit1, upper_limit2);
    lower_limit3 = min(lower_limit1, lower_limit2);
    
    x_lim = 0.0008;
    
    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(NHK1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xline(upper_limit1,'--r','LineWidth',2);
    xline(lower_limit1,'--r','LineWidth',2);
    xlim([-x_lim x_lim])
    grid on
    
    subplot(2,1,2)
    histogram(emg_pdf1,'Normalization','pdf')
    ax = gca;
    ax.FontSize = 15;
    xline(upper_limit1,'--r','LineWidth',2);
    xline(lower_limit1,'--r','LineWidth',2);
    xlabel('EMG Amplitude (V)','Fontsize',18)
    ylabel('Density','Fontsize',18)
    title('(b) EMG Distribution','Fontsize',24)
    xlim([-x_lim x_lim])
    grid on
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xlim([lower_limit1 upper_limit1])
    grid on
    
    subplot(1,2,2)
    plot(NHK1{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Element','Fontsize', 24)
    xlabel('Lags (s)','Fontsize',18)
    ylabel('X1','Fontsize',18)
    grid on

    figure(figNum)
    figNum = figNum+1;
    subplot(2,1,1)
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xline(upper_limit2,'--r','LineWidth',2);
    xline(lower_limit2,'--r','LineWidth',2);
    xlim([-x_lim x_lim])
    grid on
    
    subplot(2,1,2)
    histogram(emg_pdf2,'Normalization','pdf')
    ax = gca;
    ax.FontSize = 15;
    xline(upper_limit2,'--r','LineWidth',2);
    xline(lower_limit2,'--r','LineWidth',2);
    xlabel('EMG Amplitude (V)','Fontsize',18)
    ylabel('Density','Fontsize',18)
    title('(b) EMG Distribution','Fontsize',24)
    xlim([-x_lim x_lim])
    grid on
    
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK2{1,1})
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xlim([lower_limit2 upper_limit2])
    grid on
    
    subplot(1,2,2)
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Element','Fontsize', 24)
    xlabel('Lags (s)','Fontsize',18)
    ylabel('X1','Fontsize',18)
    grid on

    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK1{1,1})
    hold on
    plot(NHK2{1,1})
    hold off
    ax = gca;
    ax.FontSize = 15;
    title('(a) Static Nonlinearities','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xlim([lower_limit3 upper_limit3])
    grid on

    subplot(1,2,2)
    plot(NHK1{1,2})
    hold on
    plot(NHK2{1,2})
    ax = gca;
    ax.FontSize = 15;
    title('(b) Linear Elements','Fontsize', 24)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    legend('PRBS', 'Physiological','Fontsize',16)
    grid on
    hold off
    
end

tEnd = toc(tStart)/60
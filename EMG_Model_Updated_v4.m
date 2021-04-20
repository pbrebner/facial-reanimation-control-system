%EMG Model v4 (Neural Command Input and Force/Movement Output)
%Simulink Method
clc
clear all

tStart = tic;

%%
%Set initial Parameters
%set_output_noise_power = 1e-10;   %Output Noise Power
set_output_noise_power = 0;
%noise_multiplier = 10;
accuracy = [];
noise_snr = [];
output_noise_power = [];
figNum = 10;

%For Power Spectrums
Fs = 1000; 
Nfft = 10000;

%Desired Displacement Signal Type (Simple, PRBS or Physiologically Based)
simple_movement = true;
PRBS_movement = true;
%physiological_movement = true;

%Signal Amplitude (PRBS Signal)
PRBS_movement_time = 180;
variable_amplitude = true;
N = PRBS_movement_time/10;     %Number of times the amplitude randomly changes (For Variable Amplitude Only)
M = 10000;                     %Number of each random value (For Variable Amplitude Only)
PRBS_amplitude = 10;           %Amplitude (For Constant Amplitude Only)

%Set Physiological Signal Parameters
physiological_movement_time = 180;
physiological_movement_max_amplitude = 0.01;
fr = 0.1;                 %Frequency distribution mean (Hz) (Max is 1.8 Hz)
sig = 0.6;                %Std of Frequency Distribution (Hz)
W = 0.55;                  
nf = 18;                  %number of random signal changes
t_interval = physiological_movement_time/nf;    %Length of random interval (s)
chance_of_zero = false;

NHK_all = [];
Zcur_all = [];
emg_all = [];

compare_two_models = false;

%%
if compare_two_models == true
    PRBS_movement = [true false];
end

for num_signals = 1:length(PRBS_movement)
    %Create Analog Signal (Desired Movement)
    if simple_movement == true
        time = 10;
        t_total = 0:0.001:time;
        t_total_with_delay = 0:0.001:time-0.5;
        
        w = 0.8*pi;
        phi = 0;
        
        %A = ones(1,length(t_total));
        A = ones(1,length(t_total_with_delay));
%         A(1:2000) = 0;
%         A(2001:20001) = 1;
%         A(20002:40001) = 3;
%         A(40002:70001) = 6;
        
        AR = makedist('Uniform','lower',0,'upper',0.01);  %Full Amplitude Range
        AmplitudesRandom = AR;
        amp_distribution = random(AmplitudesRandom,10000,1);
        A(1:2000) = random(AmplitudesRandom,1,1);
        A(2001:4500) = random(AmplitudesRandom,1,1);
        A(4501:7000) = random(AmplitudesRandom,1,1);
        %A(7001:10001) = random(AmplitudesRandom,1,1);
        A(7001:9501) = random(AmplitudesRandom,1,1);
        
        delay = zeros(1,500);
        
        %desired_displacement = A.*square(w*t_total+phi);
        desired_displacement = [delay A.*square(w*t_total_with_delay+phi)];
        desired_displacement = max(desired_displacement,0);
        
        figure(figNum)
        figNum = figNum+1;
        plot(t_total,desired_displacement);
        ax = gca;
        ax.FontSize = 26;
        xlabel('Time (s)','Fontsize',32)
        ylabel('Desired Displacement (m)','Fontsize',32)
        title('Analog Signal of Desired Displacement','Fontsize',36)
        grid on
        
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);
        Nfft = 1000;
        [Pxx1a,f1a] = pwelch(desired_displacement(1,499:1750),Nfft,[],Nfft,Fs);
        [Pxx1b,f1b] = pwelch(desired_displacement(1,2999:4250),Nfft,[],Nfft,Fs);
        [Pxx1c,f1c] = pwelch(desired_displacement(1,5499:6750),Nfft,[],Nfft,Fs);
        [Pxx1d,f1d] = pwelch(desired_displacement(1,7999:9250),Nfft,[],Nfft,Fs);
    
    
    elseif PRBS_movement(num_signals) == true
        %(PRBS with Random or Constant Amplitude)
        t_total = 0:0.001:PRBS_movement_time;
        time = PRBS_movement_time;

        A = 0;                    %Intialize amplitude
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
        %[Pxx1,f1] = pwelch(desired_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);
              
    else   

        %(Physiologically Based Movment - Square pulses with random amplitude and random frequency)
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

        %PWR = makedist('Normal','mu',0.001,'sigma',0.005);
        %PulseWidthRandom = truncate(PWR,0.001,0.010);   % 1ms - 10 ms
        %pw_distribution = random(PulseWidthRandom,10000,1);
        %figure(107)
        %histogram(pw_distribution,100)
        %title('Pulse Width Distribution')

        desired_displacement=[0];
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
            %A = 0.01;
            %W = random(PulseWidthRandom,1,1);

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
        %[Pxx1,f1] = pwelch(desired_displacement_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
        [Pxx1,f1] = pwelch(desired_displacement_zero,Nfft,[],Nfft,Fs);

        desired_displacement = desired_displacement';

    end

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
    neural_zero = neural - mean(neural);
    %[Pxx2,f2] = pwelch(neural_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);
    
    if simple_movement == true
        Nfft = 1000;
        [Pxx2a,f2a] = pwelch(neural_zero(499:1750,1),Nfft,[],Nfft,Fs);
        [Pxx2b,f2b] = pwelch(neural_zero(2999:4250,1),Nfft,[],Nfft,Fs);
        [Pxx2c,f2c] = pwelch(neural_zero(5499:6750,1),Nfft,[],Nfft,Fs);
        [Pxx2d,f2d] = pwelch(neural_zero(7999:9250,1),Nfft,[],Nfft,Fs);
    else
        [Pxx2,f2] = pwelch(neural_zero,Nfft,[],Nfft,Fs);
    end

    if simple_movement == true
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
%         hold on
%         semilogy(f1a(1:11,1),Pxx1a(1:11,1));
%         semilogy(f1b(1:11,1),Pxx1b(1:11,1));
%         semilogy(f1c(1:11,1),Pxx1c(1:11,1));
%         semilogy(f1d(1:11,1),Pxx1d(1:11,1));
        semilogy(f1a(1:11,1),Pxx1a(1:11,1),f1b(1:11,1),Pxx1b(1:11,1),f1c(1:11,1),Pxx1c(1:11,1),f1d(1:11,1),Pxx1d(1:11,1),'Linewidth',1);
%         hold off
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
%         semilogy(f2a(1:151,1),Pxx2a(1:151,1),f2b(1:151,1),Pxx2b(1:151,1),f2c(1:151,1),Pxx2c(1:151,1),f2d(1:151,1),Pxx2d(1:151,1));
        hold off
        ax = gca;
        ax.FontSize = 12;
        title('(d) Power Spectrum of Neural Input','Fontsize',14);
        ylabel('PSD','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        legend('1st Movement','2nd Movement','3rd Movement','4th Movement','Fontsize',12)
        grid on;
        
    end

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

    %Power Pectrum of EMG
    emg_simulink_zero = emg_simulink - mean(emg_simulink);
    [Pxx_emg,f_emg] = pwelch(emg_simulink_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);

    %Plot frequency spectrum
%     figure(figNum)
%     figNum = figNum+1;
%     semilogy(f_emg,Pxx_emg);
%     ax = gca;
%     ax.FontSize = 15;
%     title('Power Spectrum of EMG','Fontsize',24);
%     ylabel('PSD (log)','Fontsize',18); 
%     xlabel('Frequency (Hz)','Fontsize',18);
%     grid on;

    force_simulink = out.EMG_Model_Force;
    
%     figure(figNum)
%     figNum = figNum+1;
%     plot(t_total,force_simulink)
%     ax = gca;
%     ax.FontSize = 26;
%     xlabel('Time (s)','Fontsize',32)
%     ylabel('Force (N)','Fontsize',32);
%     title('Muscle Force','Fontsize',36)
%     grid on
    
    %Power Pectrum of FOrce
    force_simulink_zero = force_simulink - mean(force_simulink);
    %[Pxx_force,f_force] = pwelch(force_simulink_zero,gausswin(Nfft),Nfft/2,Nfft,Fs);

    if simple_movement == true
        Nfft = 1000;
        [Pxx_force1,f_force1] = pwelch(force_simulink_zero(499:1750,1),Nfft,[],Nfft,Fs);
        [Pxx_force2,f_force2] = pwelch(force_simulink_zero(2999:4250,1),Nfft,[],Nfft,Fs);
        [Pxx_force3,f_force3] = pwelch(force_simulink_zero(5499:6750,1),Nfft,[],Nfft,Fs);
        [Pxx_force4,f_force4] = pwelch(force_simulink_zero(7999:9250,1),Nfft,[],Nfft,Fs);
    else
        [Pxx_force,f_force] = pwelch(force_simulink_zero,Nfft,[],Nfft,Fs);
    end
    
    if simple_movement == true
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
%         semilogy(f2a(1:151,1),Pxx2a(1:151,1),f2b(1:151,1),Pxx2b(1:151,1),f2c(1:151,1),Pxx2c(1:151,1),f2d(1:151,1),Pxx2d(1:151,1));
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
        %hold on
%         semilogy(f_force1(1:11,1),Pxx_force1(1:11,1));
%         semilogy(f_force2(1:11,1),Pxx_force2(1:11,1));
%         semilogy(f_force3(1:11,1),Pxx_force3(1:11,1));
%         semilogy(f_force4(1:11,1),Pxx_force4(1:11,1));
        semilogy(f_force1(1:11,1),Pxx_force1(1:11,1),f_force1(1:11,1),Pxx_force2(1:11,1),f_force1(1:11,1),Pxx_force3(1:11,1),f_force1(1:11,1),Pxx_force4(1:11,1),'Linewidth',1);
        %hold off
        ax = gca;
        ax.FontSize = 12;
        title('(d) Power Spectrum of Muscle Force','Fontsize',14);
        ylabel('PSD (log)','Fontsize',15); 
        xlabel('Frequency (Hz)','Fontsize',15);
        legend('1st Movement','2nd Movement','3rd Movement','4th Movement','Fontsize',12)
        grid on;
        
    end

    %%
    %plot desired displacement, neural command and muscle force on one figure
    figure(figNum)
    figNum = figNum+1;

    if PRBS_movement(1,num_signals) == true
        
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
        %semilogy(f1(1:2000,1),Pxx1(1:2000,1),'LineWidth',1.5);
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
        %semilogy(f1(1:400,1),Pxx1(1:400,1),'LineWidth',1.5);
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
    %semilogy(f2(1:10000,1),Pxx2(1:10000,1),'LineWidth',1.5);
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


    %%
    output_displacement_simulink = out.EMG_Model_Displacement;
    t_simulink = out.tout;

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
    set(I,'nLags',2400); %(accuracy increases if nlags is increased)
    NHK{1,2}=I;

    NHK=nlident(NHK,Zcur);
    %NHK = normGainLE(NHK);
    NHK = normCoefLE(NHK);
    %NHK = normCoefNLE(NHK);
    
    upper_limit = 0.00067;
    lower_limit = -0.00067;

    %Plots Hammerstein System and EMG Distribution with red line limits 
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(2,2,1)
%     plot(NHK{1,1});
%     ax = gca;
%     ax.FontSize = 10;
%     title('(a) Static Nonlinearity','Fontsize', 14)
%     xlabel('Input','Fontsize',12)
%     ylabel('Output','Fontsize',12)
%     xline(upper_limit,'--r','LineWidth',1);
%     xline(lower_limit,'--r','LineWidth',1);
%     grid on

%     subplot(2,2,[2 4])
%     plot(NHK{1,2});
%     ax = gca;
%     ax.FontSize = 10;
%     title('(c) Linear Element','Fontsize', 14)
%     ylabel('X1', 'Fontsize',12)
%     xlabel('Time (s)','Fontsize',12)
%     grid on

%     subplot(2,2,3)
%     histogram(emg_pdf,'Normalization','pdf')
%     ax = gca;
%     ax.FontSize = 10;
%     xline(upper_limit,'--r','LineWidth',1);
%     xline(lower_limit,'--r','LineWidth',1);
%     title('(b) EMG Distribution','Fontsize',14)
%     xlabel('Volts (V)','Fontsize',12)
%     ylabel('Density','Fontsize',12)
%     grid on

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

    %Plots Linear Hammerstein Part
    %figure(figNum)
    %figNum = figNum+1;
    %plot(NHK{1,2});
    %grid on

    %Plots NL and Linear Hammerstein Parts with Red Line Limits on NL
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(1,2,1)
%     plot(NHK{1,1});
%     ax = gca;
%     ax.FontSize = 14;
%     title('(a) Static Nonlinearity','Fontsize', 24)
%     xlabel('Input','Fontsize',18)
%     ylabel('Output','Fontsize',18)
%     xline(upper_limit,'--r','LineWidth',1);
%     xline(lower_limit,'--r','LineWidth',1);
%     grid on
% 
%     subplot(1,2,2)
%     plot(NHK{1,2});
%     ax = gca;
%     ax.FontSize = 14;
%     title('(b) Linear Element','Fontsize', 24)
%     ylabel('X1', 'Fontsize',18)
%     xlabel('Lag (s)','Fontsize',18)
%     grid on
    
    %Plot NL and Linear of Hammerstein System with only High Quality Values
    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(NHK{1,1});
    ax = gca;
    ax.FontSize = 15;
    %yticks([0 0.4 0.8 1.2])
    title('(a) Static Nonlinearity','Fontsize', 24)
    xlabel('EMG Input, E(t) (V)','Fontsize',18)
    ylabel('Transformed EMG (V)','Fontsize',18)
    xlim([lower_limit upper_limit])
    %xline(upper_limit,'--r','LineWidth',1);
    %xline(lower_limit,'--r','LineWidth',1);
    grid on

    subplot(1,2,2)
    plot(NHK{1,2})
    ax = gca;
    ax.FontSize = 15;
    %yticks([0 0.1 0.2 0.3])
    title('(b) Linear Element','Fontsize', 24)
    ylabel('X1', 'Fontsize',18)
    xlabel('Lags (s)','Fontsize',18)
    grid on

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
    
    %Box Plot (Expanded Time Scale) of previous plot
%     figure(figNum);
%     figNum = figNum+1;
%     plot(t_total,pred);
%     hold on
%     plot(t_total, output_displacement_simulink)
%     ax = gca;
%     ax.FontSize = 14;
%     hold off
%     title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 24)
%     xlabel('Time (s)', 'Fontsize', 18)
%     ylabel('Displacement (m)', 'Fontsize', 18)
%     legend('Predicted', 'Observed')
%     %create a new pair of axes inside current figure
%     axes('position',[.65 .175 .25 .25])
%     box on % put box around new pair of axes
%     indexOfInterest = (t_total < 100) & (t_total > 60); 
%     plot(t_total(indexOfInterest),pred(indexOfInterest)) % plot on new axes
%     hold on
%     plot(t_total(indexOfInterest),output_displacement_simulink(indexOfInterest)) % plot on new axes
%     hold off
%     axis tight

    %Plots the Residuals, Residuals Distribution, and Residuals Power
    %Spectrum
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
%     S = spect(R_zero);
%     S_frequency = 0.0556:0.0556:0.0556*length(S);
%     subplot(2,2,4)
%     semilogy(S_frequency(:,1:150),S(1:150,:),'LineWidth',1.5);
%     ax = gca;
%     ax.FontSize = 15;
%     title('(c) Power Spectrum of Residuals','Fontsize',20);
%     ylabel('PSD (log)','Fontsize',20); 
%     xlabel('Frequency (Hz)','Fontsize',20);
%     grid on
     
    %[PxxR,fR] = pwelch(double(R_zero),gausswin(Nfft),Nfft/2,Nfft,Fs);
    [PxxR,fR] = pwelch(double(R_zero),[],[],[],Fs);
    subplot(2,2,4)
    %semilogy(fR(1:110,:),PxxR(1:110,:),'LineWidth',1.5);
    semilogy(fR(1:650,:),PxxR(1:650,:),'LineWidth',1.5);
    ax = gca;
    ax.FontSize = 15;
    title('(c) Power Spectrum of Residuals','Fontsize',20);
    ylabel('PSD (log)','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on

    %Plots the Residuals, Residuals Distribution (With Normal Distribution), and Residuals Power
    %Spectrum
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
    R_double = double(R_zero);
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
    %[res_corr,lags] = xcov(R,maxlag,'coeff');
    lags = lags/Fs;
    subplot(2,2,4),plot(lags,res_corr,'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 15;
    title('(c) Auto-Correlation of Residuals','Fontsize',20)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Correlation','Fontsize',20)
    grid on
    
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
    hold on
    R_double = double(R);
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

    L = length(R_double)-1;
    maxlag = L*0.10;
    [res_corr,lags] = xcorr(R_double,maxlag,'normalized');
    %[res_corr,lags] = xcov(R,maxlag,'coeff');
    lags = lags/Fs;
    subplot(2,2,4),plot(lags,res_corr,'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 15;
    title('(c) Auto-Correlation Residuals','Fontsize',20)
    xlabel('Time (s)','Fontsize',20)
    ylabel('Correlation','Fontsize',20)
    grid on
    
    %Plots the Superimposed, Residuals, Residuals Distribution, and Residuals Power
    %Spectrum
%     figure(figNum)
%     figNum = figNum+1;
%     subplot(3,1,1)
%     plot(t_total,pred);
%     hold on
%     plot(t_total, output_displacement_simulink)
%     ax = gca;
%     ax.FontSize = 15;
%     hold off
%     title(['(a) Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 22)
%     xlabel('Time (s)', 'Fontsize', 18)
%     ylabel('Displacement (m)', 'Fontsize', 18)
%     legend('Predicted', 'Observed', 'Fontsize', 16)
%     
%     subplot(4,1,2)
%     plot(R)
%     ax = gca;
%     ax.FontSize = 14;
%     title('Residuals of EMG/Whisk Hammerstein Model','Fontsize',24)
%     xlabel('Time (s)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',20)
%     grid on
% 
%     subplot(3,1,2)
%     p = pdf(R);
%     plot(p)
%     ax = gca;
%     ax.FontSize = 15;
%     xlabel('Displacement (m)','Fontsize',18)
%     ylabel('Density','Fontsize',18)
%     title('(b) Residual Distribution','Fontsize',22)
%     grid on
% 
%     S = spect(R);
%     subplot(3,1,3)
%     plot(S(1:50,:),'LineWidth',1.5);
%     ax = gca;
%     ax.FontSize = 15;
%     title('Power Spectrum of Residuals','Fontsize',22);
%     ylabel('PSD','Fontsize',18); 
%     xlabel('Frequency (Hz)','Fontsize',18);
%     grid on
    

    %Plots the Superimposed, Residuals Distribution, and Residuals Power
    %Spectrum
    figure(figNum)
    figNum = figNum+1;
    subplot(2,2,[1 2])
    plot(t_total,pred);
    hold on
    plot(t_total, output_displacement_simulink)
    ax = gca;
    ax.FontSize = 16;
    hold off
    title(['(a) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 20)
    xlabel('Time (s)', 'Fontsize', 20)
    ylabel('Healthy Displacement, Pos_H(t) (m)', 'Fontsize', 20)
    legend('Predicted', 'Observed', 'Fontsize', 18)
    
%     subplot(4,1,2)
%     plot(R)
%     ax = gca;
%     ax.FontSize = 14;
%     title('Residuals of EMG/Whisk Hammerstein Model','Fontsize',24)
%     xlabel('Time (s)','Fontsize',20)
%     ylabel('Displacement (m)','Fontsize',20)
%     grid on

    subplot(2,2,3)
    p = pdf(R);
    plot(p)
    ax = gca;
    ax.FontSize = 16;
    xlabel('Displacement (m)','Fontsize',20)
    ylabel('Density','Fontsize',20)
    title('(b) Residual Distribution','Fontsize',18)
    grid on
    
    R_zero = R - mean(R);
    S = spect(R_zero);
    subplot(2,2,4)
    plot(S(1:150,:));
    ax = gca;
    ax.FontSize = 16;
    title('(c) Power Spectrum of Residuals','Fontsize',18);
    ylabel('PSD','Fontsize',20); 
    xlabel('Frequency (Hz)','Fontsize',20);
    grid on
    
    NHK_all = [NHK_all NHK];
    Zcur_all = [Zcur_all Zcur];
    emg_all = [emg_all emg_simulink];
    
end

%%

if compare_two_models == true
    NHK1 = NHK_all(1);
    NHK2 = NHK_all(2);
    
    emg_pdf1 = emg_all(:,1);
    emg_pdf2 = emg_all(:,2);
    emg_pdf1(emg_pdf1==0) = []; %Removing the zero values
    emg_pdf2(emg_pdf2==0) = []; %Removing the zero values
    
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
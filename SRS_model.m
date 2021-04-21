%Paralyzed Model Updated (Stimulus Input and Force/Movement Output)
%Simulink Method

clc
clear all

tStart = tic;

%%
%Set initial Parameters
%set_output_noise_power = 1e-10;        %Output Noise Power
set_output_noise_power = 0;
%noise_multiplier = 10;
accuracy = [];
noise_snr = [];
output_noise_power = [];
figNum = 100;

%For Power Spectrums or FFT
Fs = 1000; 
Nfft = 10000;

simple_stimulus = false;
sinusoidal_stimulus = false;
PRBS_stimulus = true;
% physiological_stimulus = true;

%simple stimulus (Constant Amp square pulse with Modulated Frequency)
desired_displacement_frequency = 0.5;
simple_stimulus_time = 9.999;
simple_stimulus_amplitude = 6.5;
simple_stimulus_pulsewidth = 0.001;

%Sinusoidal (Simple) Stimulus (Constant Freq Sine with Modulated Amp)
desired_displacement_frequency = 0.5;
sinusoidal_stimulus_max_amplitude = 0.02;
sine_stimulus_time = 9.999;
sine_stimulus_frequency = 50;
uniphasic = true;

%PRBS Stimulus
PRBS_stimulus_time = 480;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm
% PRBS_stimulus_amplitude = 6.5;  %volts

%Set Physiological Signal Parameters
physiological_stimulus_time = 480;
physiological_stimulus_max_amplitude = 0.02;
fr = 0.1;                 %Frequency distribution mean (Hz)
sig = 0.8;                %Std of Frequency Distribution (Hz)
W = 0.45;                  
nf = physiological_stimulus_time/10;                  %number of random signal changes
t_interval = physiological_stimulus_time/nf;    %Length of random interval (s)
chance_of_zero = false;

LNL_all = [];
Hammerstein_all = [];
Weiner_all = [];
IRF_model_all = [];
NLN_all = [];
Zcur_all = [];

compare_two_models = true;

LNL_model = false;
Hammerstein_model = false;
Weiner_model = true;
Linear_IRF_model = true;

%%
if compare_two_models == true
    PRBS_stimulus = [true false];
end

for num_signals = 1:length(PRBS_stimulus)
    %Create Analog Signal (Desired Movement)
    if simple_stimulus == true

        t_total = 0:0.001:simple_stimulus_time;
        time = simple_stimulus_time;

        AmplitudesRandom = makedist('Uniform','lower',0,'upper',0.01);
        desired_displacement_amplitude = [];

        for ii = 1:ceil(time)/2
            rand_amp = random(AmplitudesRandom,1,1);
            desired_displacement_amplitude = [desired_displacement_amplitude ones(1,2000)*rand_amp];
        end
        %desired_displacement_amplitude = [ones(1,2000)*0.002 ones(1,2000)*0.004 ones(1,2000)*0.006 ones(1,2000)*0.008 ones(1,2000)*0.01];

        desired_displacement = max(-desired_displacement_amplitude.*square(2*pi*desired_displacement_frequency.*t_total),0);
%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,desired_displacement)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Displacement (m)','Fontsize',18);
%         title('Desired Displacement','Fontsize',24)
%         grid on

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
%         plot(f1(1:50,:),Pxx1(1:50,:)) 
%         title('FFT of Desired Displacement')
%         xlabel('f (Hz)')
%         ylabel('|Pxx1(f)|')

        %Power Spectrum of Desired Displacement
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stim_amplitude = zeros(1,length(desired_displacement));
        for i = 1:length(desired_displacement)
            if desired_displacement(i) == 0
                stim_amplitude(i) = 0;
            else
                stim_amplitude(i) = simple_stimulus_amplitude;
            end
        end

        stim_frequency = desired_displacement*5000;  %Hz

    %     input_stimulus = stim_amplitude.*square(2*pi*stim_frequency.*t_total);
    %     input_stimulus = max(input_stimulus,0);
    % 
    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(t_total,input_stimulus)
    %     ax = gca;
    %     ax.FontSize = 15;
    %     xlabel('Time (s)','Fontsize',18)
    %     ylabel('Volts (V)','Fontsize',18);
    %     title('Stimulus Input','Fontsize',24)
    %     grid on

        D = [];
        for j = 1:1001:length(desired_displacement)
            if desired_displacement(j) == 0
    %             D = [D zeros(1,1000)];
                D = D;
            else
                g = 1/stim_frequency(j);
                D = [D (j/1000+1/stim_frequency(j)):g:1+j/1000];
            end
        end

        data = (simple_stimulus_amplitude*pulstran(t_total,D,@rectpuls,simple_stimulus_pulsewidth))';
        input_stimulus = data';
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,input_stimulus)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Volts (V)','Fontsize',18);
%         title('Stimulus Input','Fontsize',24)
%         grid on

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
%         plot(f2(1:800,:),Pxx2(1:800,:)) 
%         title('FFT of Input Stimulus')
%         xlabel('f (Hz)')
%         ylabel('|Pxx2(f)|')

        %Power Spectrum of Electrical Stimulus
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

    %     g = 1/simple_stimulus_frequency;
    %     D = (1:g:simple_stimulus_time)';     % pulse delay times
    %     data = (simple_stimulus_amplitude*pulstran(t_total,D,@rectpuls,simple_stimulus_pulsewidth))';
    %     
    %     input_stimulus = [data; zeros(1000,1)]';
    %     t_total = 0:0.001:simple_stimulus_time+1;
    %     time = simple_stimulus_time+1;
    %     figure(figNum)
    %     figNum = figNum+1;
    %     plot(t_total,input_stimulus)
    %     
         stimulus_simulink = [t_total' input_stimulus'];

    elseif sinusoidal_stimulus == true

        t_total = 0:0.001:sine_stimulus_time;
        time = sine_stimulus_time;

        AmplitudesRandom = makedist('Uniform','lower',0,'upper',sinusoidal_stimulus_max_amplitude);
        desired_displacement_amplitude = [];

        for ii = 1:ceil(time)/2
            rand_amp = random(AmplitudesRandom,1,1);
            desired_displacement_amplitude = [desired_displacement_amplitude ones(1,2000)*rand_amp];
        end

        desired_displacement = max(-desired_displacement_amplitude.*square(2*pi*desired_displacement_frequency.*t_total),0);

%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,desired_displacement);
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
        Pxx1 = Pxx1';

        f1 = (Fs*(0:(L/2))/L)';
%         figure(figNum)
%         figNum = figNum+1;
%         plot(f1(1:50,:),Pxx1(1:50,:)) 
%         title('FFT of Desired Displacement')
%         xlabel('f (Hz)')
%         ylabel('|Pxx1(f)|')

        %Power Spectrum of Desired Displacement
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stim_amplitude = desired_displacement*170;
        if uniphasic == true
            input_stimulus = max(stim_amplitude.*sin(2*pi*sine_stimulus_frequency.*t_total),0);
        else
            input_stimulus = stim_amplitude.*sin(2*pi*sine_stimulus_frequency.*t_total);
        end
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];


%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,input_stimulus)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Volts (v)','Fontsize',18);
%         title('Electrical Stimulus Input','Fontsize',24)
%         grid on

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
%         plot(f2(1:800,:),Pxx2(1:800,:)) 
%         title('FFT of Input Stimulus')
%         xlabel('f (Hz)')
%         ylabel('|Pxx2(f)|')

        %Power Pectrum of Neural Input
        %[Pxx2,f2] = pwelch(input_stimulus,gausswin(Nfft),Nfft/2,Nfft,Fs);

        stimulus_simulink = [t_total' input_stimulus'];

    elseif PRBS_stimulus(num_signals) == true

        t_total = 0:0.001:PRBS_stimulus_time;
        time = PRBS_stimulus_time;

        A = [0];                    %Intialize amplitude
        if variable_amplitude == true   
            for k = 1:N
                if k == 1
                    R = PRBS_amplitude;
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
%         figure(figNum)
%         figNum = figNum+1;
%         plot(f1(1:5000,:),Pxx1(1:5000,:)) 
%         title('FFT of Desired Displacement')
%         xlabel('f (Hz)')
%         ylabel('|Pxx1(f)|')

        %Power Pectrum of Desired Displacement
        %[Pxx1,f1] = pwelch(desired_displacement,gausswin(Nfft),Nfft/2,Nfft,Fs);

        %stim_frequency = desired_displacement*5000;  %Hz
        stim_frequency = 50;
        %stim_amplitude = PRBS_stimulus_amplitude;
        stim_amplitude = desired_displacement*170;

        %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,input_stimulus)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Volts (v)','Fontsize',18);
%         title('Electrical Stimulus Input','Fontsize',24)
%         grid on

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
        FrequenciesRandom_max = 2.21;
        FrequenciesRandom = truncate(FR,0,FrequenciesRandom_max);
        freq_distribution = random(FrequenciesRandom,10000,1);
        figure(figNum)
        figNum = figNum+1;
        histogram(freq_distribution,100)
        title('Frequency Distribution')
        xlabel('Frequency (Hz)')

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

%         figure(figNum)
%         figNum = figNum+1;
%         subplot(2,1,1)
%         plot(t_total,desired_displacement)
%         ax = gca;
%         ax.FontSize = 14;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Desired Displacement (m)','Fontsize',18)
%         title('Analog Signal of Desired Movement','Fontsize',24)
%         grid on

        %Plot frequency spectrum
        % subplot(2,1,2)
        % semilogy(f1(1:6,1),Pxx1(1:6,1));
        % ax = gca;
        % ax.FontSize = 14;
        % title('Power Spectrum of Desired Displacement','Fontsize',24);
        % ylabel('PSD (log)','Fontsize',18); 
        % xlabel('Frequency (Hz)','Fontsize',18);
        % grid on;

        desired_displacement = desired_displacement';

        stim_amplitude = desired_displacement*170;  %mV
    %     stim_amplitude = physiological_stimulus_amplitude;
    %     stim_frequency = desired_displacement*5000;  %Hz
        stim_frequency = 50;

        %input_stimulus = max(stim_amplitude.*square(2*pi*stim_frequency.*t_total),0);
        input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);
        
        amplitude_modulation = stim_amplitude;
        amplitude_modulation_simulink = [t_total' amplitude_modulation'];

%         figure(figNum)
%         figNum = figNum+1;
%         plot(t_total,input_stimulus)
%         ax = gca;
%         ax.FontSize = 15;
%         xlabel('Time (s)','Fontsize',18)
%         ylabel('Volts (V)','Fontsize',18);
%         title('Stimulus Input','Fontsize',24)
%         grid on

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

%     figure(figNum)
%     figNum = figNum+1;
%     plot(t_total,force_simulink)
%     ax = gca;
%     ax.FontSize = 15;
%     xlabel('Time (s)','Fontsize',18)
%     ylabel('Force (N)','Fontsize',18);
%     title('Muscle Force','Fontsize',24)
%     grid on

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
    %Input Stimulus and Output Displacement
    input_stimulus = out.Paralyzed_Model_Stimulus;
    output_displacement_simulink = out.Paralyzed_Model_Displacement;
    t_simulink = out.tout;

    %Zcur = [input_stimulus',output_displacement_simulink];
    %Zcur = [input_stimulus,output_displacement_simulink];
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
    
    %%
    %Plot Desired Displacement, Electrical Stimulus, Muscle Force and Output
    %Displacement
    
    L = length(input_stimulus);

    Y = fft(input_stimulus);
    P2 = abs(Y/L);
    Pxx2 = P2(1:L/2+1);
    Pxx2(2:end-1) = 2*Pxx2(2:end-1);

    f2 = (Fs*(0:(L/2))/L)';

    if sinusoidal_stimulus == true || simple_stimulus == true

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
        
        figure(figNum)
        figNum = figNum+1;
        subplot(3,1,1)
        plot(t_total,amplitude_modulation);
        ax = gca;
        ax.FontSize = 13;
        %xlabel('Time (s)','Fontsize',15)
        ylabel('Ampltiude (V)','Fontsize',16)
        title('(a) Amplitude Modulation Input, A(t)','Fontsize',18)
        grid on
        
        subplot(3,1,2)
        plot(t_total,input_stimulus)
        ax = gca;
        ax.FontSize = 13;
        %xlabel('Time (s)','Fontsize',15)
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
        
        %Plot Desired Displacement, input stimulus, muscle force, and
        %output displacement
        figure(figNum)
        figNum = figNum+1;

        subplot(2,2,1)
        if PRBS_stimulus == true
            plot(t_total,desired_displacement);
            ax = gca;
            ax.FontSize = 16;
            xlabel('Time (s)','Fontsize',16)
            ylabel('Displacement (m)','Fontsize',16)
            title('(a) PRBS Desired Displacement','Fontsize',16)
            grid on
        else
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

    %%
    %Calculate Signal to Noise Ratio
    output_noise_simulink = out.Output_Noise;

    signal_to_noise = snr(output_displacement_simulink, output_noise_simulink);
    noise_snr = [noise_snr signal_to_noise];

    %%
    %Model Identification (LNL Model or NLN Model)

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
        title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
        xlabel('Time (s)', 'Fontsize', 24)
        ylabel('Displacement (m)', 'Fontsize', 24)
        legend('Predicted', 'Observed', 'Fontsize', 20)

        %Plots Residuals, Residual Distribution and Residual Spectrum
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,[1 2])
        plot(R)
        ax = gca;
        ax.FontSize = 15;
        title('(a) Residuals from LNL Identification','Fontsize',18)
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        grid on

        subplot(2,2,3)
    %     p = pdf(R);
    %     plot(p)
        histogram(double(R))
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',18)
        grid on

        S = spect(R);
        subplot(2,2,4)
        S_frequency = 0.0556:0.0556:0.0556*length(S);
        subplot(2,2,4)
        semilogy(S_frequency(:,1:180),S(1:180,:),'LineWidth',1.5);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on

        disp(LNL.idMethod)
        
        LNL_all = [LNL_all LNL];
        Zcur_all = [Zcur_all Zcur];
        
    elseif Hammerstein_model == true
        
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});

        Hammerstein=nlbl;
        set(Hammerstein,'idMethod','hk','displayFlag',true,'threshNSE',.001);
        I=Hammerstein{1,2};
        set(I,'nLags',2400,'nSides',1); %(accuracy increases if nlags is increased)
        Hammerstein{1,2}=I;

        Hammerstein=nlident(Hammerstein,Zcur);
        
        figure(figNum)
        figNum = figNum+1;
        subplot(1,2,1)
        plot(Hammerstein{1,1});
        ax = gca;
        ax.FontSize = 15;
        title('(a) Static Nonlinearity','Fontsize', 24)
        xlabel('Input','Fontsize',18)
        ylabel('Output','Fontsize',18)
        grid on

        subplot(1,2,2)
        plot(Hammerstein{1,2})
        ax = gca;
        ax.FontSize = 15;
        title('Linear Element','Fontsize', 24)
        ylabel('X1', 'Fontsize',18)
        xlabel('Lags (s)','Fontsize',18)
        grid on
        
        figure(figNum)
        figNum = figNum+1;
        [R, V, yp] = nlid_resid(Hammerstein,Zcur);

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
        title(['Superimposed, VAF = ' num2str(V) '%'], 'Fontsize', 28)
        xlabel('Time (s)', 'Fontsize', 24)
        ylabel('Displacement (m)', 'Fontsize', 24)
        legend('Predicted', 'Observed', 'Fontsize', 20)

        %Plots Residuals, Residual Distribution and Residual Spectrum
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,[1 2])
        plot(R)
        ax = gca;
        ax.FontSize = 15;
        title('(a) Residuals from Hammerstein Identification','Fontsize',18)
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        grid on

        subplot(2,2,3)
    %     p = pdf(R);
    %     plot(p)
        histogram(double(R))
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',18)
        grid on

        S = spect(R);
        subplot(2,2,4)
        S_frequency = 0.0556:0.0556:0.0556*length(S);
        subplot(2,2,4)
        semilogy(S_frequency(:,1:180),S(1:180,:),'LineWidth',1.5);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
    
        Hammerstein_all = [Hammerstein_all Hammerstein];
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
    %     p = pdf(R);
    %     plot(p)
        histogram(double(R))
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',18)
        grid on
        
        R_zero = R - mean(R);
        S = spect(R_zero);
        subplot(2,2,4)
        S_frequency = 0.0556:0.0556:0.0556*length(S);
%         subplot(2,2,4)
%         semilogy(S_frequency(:,1:180),S(1:180,:),'LineWidth',1.5);
%         ax = gca;
%         ax.FontSize = 15;
%         title('(c) Residual Power Spectrum','Fontsize',18);
%         ylabel('PSD (log)','Fontsize',20); 
%         xlabel('Frequency (Hz)','Fontsize',20);
%         grid on
        
        Nfft = 200000;
        [PxxR,fR] = pwelch(double(R_zero),Nfft,[],Nfft,Fs);
        subplot(2,2,4)
        semilogy(fR(1:1300),PxxR(1:1300,:),'LineWidth',1);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
        
        Weiner_all = [Weiner_all Weiner];
        Zcur_all = [Zcur_all Zcur];
        
    elseif Linear_IRF_model == true
        % Identify a two-sided IRF Model
        
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
        title(['Superimposed (IRF Model), VAF = ' num2str(V) '%'], 'Fontsize', 28)
        xlabel('Time (s)', 'Fontsize', 24)
        ylabel('Displacement (m)', 'Fontsize', 24)
        legend('Predicted', 'Observed', 'Fontsize', 20)

        %Plots Residuals, Residual Distribution and Residual Spectrum
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,[1 2])
        plot(R)
        ax = gca;
        ax.FontSize = 15;
        title('(a) Residuals from IRF Identification','Fontsize',18)
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        grid on

        subplot(2,2,3)
    %     p = pdf(R);
    %     plot(p)
        histogram(double(R))
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',18)
        grid on

        S = spect(R);
        subplot(2,2,4)
        S_frequency = 0.0556:0.0556:0.0556*length(S);
        subplot(2,2,4)
        semilogy(S_frequency(:,1:180),S(1:180,:),'LineWidth',1.5);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
        
        if compare_two_models == true && PRBS_stimulus(num_signals) == true
            IRF1 = IRF_model;
        elseif compare_two_models == true && PRBS_stimulus(num_signals) == false
            IRF2 = IRF_model;
        else
            IRF1 = IRF_model;
        end
        
        Zcur_all = [Zcur_all Zcur];

    else
        set(Zcur, 'chanNames', {'Predicted (m)' 'Displacement (m)'});
        Zcur_NLN = iddata(output_displacement_simulink,input_stimulus',Fs);
        %'domainIncr',0.001,'comment','Input Stimulus, Output Displacement','chanNames', {'Stimulus (V)' 'Displacement (m)'});

        set(Zcur_NLN,'InputName','Input Stimulus','OutputName','Output Displacement');

        NLN = nlhw(Zcur_NLN,[2 2 2]);
        plot(NLN);

        figure(figNum)
        figNum = figNum+1; 
        compare(Zcur_NLN,NLN,'r');
        grid on

        [y,fit,ic] = compare(Zcur_NLN,NLN);

        pred = y.y(:,1);
        V = vaf(output_displacement_simulink,pred);
        accuracy = [accuracy V]; 

        R = output_displacement_simulink - pred;

        %Plots just the superimposed part with %VAF in the title
        figure(figNum);
        figNum = figNum+1;
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

        %Plots Residuals, Residual Distribution and Residual Spectrum
        figure(figNum)
        figNum = figNum+1;
        subplot(2,2,[1 2])
        plot(R)
        ax = gca;
        ax.FontSize = 15;
        title('(a) Residuals from NLN Identification','Fontsize',18)
        xlabel('Time (s)','Fontsize',20)
        ylabel('Displacement (m)','Fontsize',20)
        grid on

        subplot(2,2,3)
    %     p = pdf(R);
    %     plot(p)
        histogram(R)
        ax = gca;
        ax.FontSize = 15;
        xlabel('Displacement (m)','Fontsize',20)
        ylabel('Density','Fontsize',20)
        title('(b) Residual Distribution','Fontsize',18)
        grid on

        S = spect(R);
        subplot(2,2,4)
        S_frequency = 0.0556:0.0556:0.0556*length(S);
        subplot(2,2,4)
        semilogy(S_frequency(:,1:180),S(1:180,:),'LineWidth',1.5);
        ax = gca;
        ax.FontSize = 15;
        title('(c) Residual Power Spectrum','Fontsize',18);
        ylabel('PSD (log)','Fontsize',20); 
        xlabel('Frequency (Hz)','Fontsize',20);
        grid on
        
        %NLN_all = [NLN_all NLN];
        %Zcur_all = [Zcur_all Zcur];

    end
    
    
end

%%
if compare_two_models == true && LNL_model == true
    
    LNL1 = LNL_all(1);
    LNL2 = LNL_all(2);
    
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
    
elseif compare_two_models == true && LNL_model == false && Hammerstein_model == true
    
    Hammerstein1 = Hammerstein_all(1);
    Hammerstein2 = Hammerstein_all(2);
    
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
    
elseif compare_two_models == true && LNL_model == false && Hammerstein_model == false && Weiner_model == true
        
    Weiner1 = Weiner_all(1);
    Weiner2 = Weiner_all(2);
    
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
    
elseif compare_two_models == true && LNL_model == false && Hammerstein_model == false && Weiner_model == false && Linear_IRF_model == true
    
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


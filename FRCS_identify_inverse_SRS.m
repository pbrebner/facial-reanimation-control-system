%% FRCS: Identify the Inverse SRS

%INSTRUCTIONS: Must Run FRCS_identify_models before running this script

%Takes the SRS model identified in FRCS_identify_models and estimates the
%inverse SRS

%Identifying the Inverse SRS Model:
%   1. Simulate the response of the SRS you identified to a PRBS input to generate a simulated output.
%   2. Then identify an inverse SRS model between the simulated output, as
%   input, and the simulated PRBS input, as output.

%When running the script, you need to provide the following input:
% 1. Which SRS Model do you want to Inverse? PRBS/Phys
%       Select which identified SRS model you want to Inverse. Default is
%       Physiological or if only one model was identified it'll skip this
%       question and default to the model identified
% 2. Inverse SRS Model Structure? LNL/Hammerstein/Wiener/IRF
%       Select the Inverse SRS Model Structure. Default is Hammerstein
% 3. Number of Validation Trials?
%       Number of Validation trials to test the Inverse SRS

%% User Input Prompts

if compare_two_models == true
    prompt1 = 'Which SRS Model do you want to Inverse? PRBS/Phys [Phys]: ';
    str1 = input(prompt1,'s');
    if ~strcmp(str1,'PRBS') & ~strcmp(str1,'Phys') & ~isempty(str1)
        disp('Invalid Input')
        return
    elseif isempty(str1)
        str1 = 'Phys';
    end
    
    if strcmp(str1,'PRBS')
        if strcmp(model_type,'IRF')
            SRS = SRS_all{1};
        else
            SRS = SRS_all(1);
        end
        disp(['Inversing ' model_type ' SRS Model Identified from PRBS Input'])
    elseif strcmp(str1,'Phys')
        if strcmp(model_type,'IRF')
            SRS = SRS_all{2};
        else
            SRS = SRS_all(2);
        end
        disp(['Inversing ' model_type ' SRS Model Identified from Physiological Input'])
    end
    SRS_inverse_input_type = str1;
else
    if strcmp(model_type,'IRF')
        SRS = SRS_all{1};
    else
        SRS = SRS_all(1);
    end
    disp(['Inversing ' model_type ' SRS Model Identified from ' input_type ' Input'])
    SRS_inverse_input_type = input_type;
end

prompt2 = 'Inverse SRS Model Structure? LNL/Hammerstein(Hamm)/Wiener/IRF [Hamm]: ';
str2 = input(prompt2,'s');
if ~strcmp(str2,'LNL') & ~strcmp(str2,'Hamm') & ~strcmp(str2,'Wiener') & ~strcmp(str2,'IRF') & ~isempty(str2)
    disp('Invalid Input')
    return
elseif isempty(str2)
    str2 = 'Hamm';
end
SRS_inverse_model_structure = str2;

prompt3 = 'Number of Validation Trials? 1-30 [1]: ';
str3 = input(prompt3);
if str3<1 | str3>30
    disp('Invalid Input')
    return
elseif isempty(str3)
    str3 = 1;
end

tStart = tic;

%% Set Initial Parameters

figNum = 100;
Fs = 1000;

%Simulated Input (PRBS) Parameters
PRBS_stimulus_time = 480;
variable_amplitude = true;      %Can either be constant or variable amplitude
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20;            %PRBS Signal Amplitude (mm)

%Number of Signals (Identification and Validation)
num_signals = str3+1;

%Inverse Model Structure
SRS_inverse_LNL = false;
SRS_inverse_Hamm = false;
SRS_inverse_Weiner = false;
SRS_inverse_IRF = false;

if strcmp(str2,'LNL')
    SRS_inverse_LNL = true;
elseif strcmp(str2,'Hamm')
    SRS_inverse_Hamm = true;
elseif strcmp(str2,'Wiener')
    SRS_inverse_Weiner = true;
elseif strcmp(str2,'IRF')
    SRS_inverse_IRF = true;
end

%% Simulated Input (PRBS Input)

simulated_input = [];

for signal = 1:num_signals

    t_total = 0:0.001:PRBS_stimulus_time;
    time = PRBS_stimulus_time;

    A = [0];                                    %Intialize amplitude
    if variable_amplitude == true   
        for k = 1:N
            if k == 1
                R = PRBS_amplitude;
            else
                R = rand(1,1)*PRBS_amplitude;   %Randomly generate a number between 0 and PRBS Amplitude
            end

            for j = 1:M
                A = [A R];
            end
        end
    else
        A = PRBS_amplitude;              %Else set as Constant Amplitude
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

    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    simulated_input(:,signal) = stim_amplitude;

end

%% Simulate the reponse of the SRS model to the Input Realizations

%Simulated Outputs of the SRS Model
simulated_output = [];
simulated_outputs = [];

%Runs through all input realizations
for signal = 1:num_signals

    simulated_output = nlsim(SRS,simulated_input(:,signal));
    set(simulated_output, 'domainIncr',0.001);
    
    simulated_output = max(simulated_output,0);
    simulated_outputs(:,signal) = simulated_output.dataSet;

end

%% Identify a model between a simulated output (as input) and a simulated input (as output)

%Input/Output for Inverse SRS Identification
Zcur_simulated = [simulated_outputs(1000:end,1) simulated_input(1000:end,1)];
Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

%plot the Input/Output
figure(figNum)
figNum = figNum+1;
subplot(2,1,1)
plot(t_total(1,1:end-999),simulated_input(1000:end,1))
ax = gca;
ax.FontSize = 15;
title('Simulated Amplitude Modulation, A(t)','Fontsize',20)
ylabel('Amplitude (V)','Fontsize',18)
grid on

subplot(2,1,2)
plot(t_total(1,1:end-999),simulated_outputs(1000:end,1))
ax = gca;
ax.FontSize = 15;
title('Simulated Paralyzed Displacement, Pos_P(t)','Fontsize',20)
ylabel('Displacement (m)','Fontsize',18)
xlabel('Time (s)','Fontsize',18)
grid on

%Outputs to other input realizations 
validation_output = [];

for signal = 2:num_signals

    Validation_Output = nldat(simulated_outputs(1000:end,signal));
    set(Validation_Output,'domainIncr',0.001);
    
    validation_output(:,signal-1) = Validation_Output.dataSet;

end

%% Identify the Inverse SRS

%Number of Lags of Inverse SRS IRF(s)
nLags = 400;

if SRS_inverse_LNL == true

    SRS_inverse = lnlbl;  %LNL Model
    set(SRS_inverse,'idMethod','hk','hkTolerance', 0.1,...
        'nhkMaxIts', 4, 'nhkMaxInner', 4);

    I1 = irf;
    set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
    P = polynom;
    I3 = irf;
    set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
    SRS_inverse.elements = {I1 P I3};
    
    SRS_inverse = nlident(SRS_inverse,Zcur_simulated);
    
elseif SRS_inverse_Hamm == true
    
    SRS_inverse = nlbl;  %Hammerstein
    set(SRS_inverse,'idMethod','hk','displayFlag',true,'threshNSE',.001);
    I2 = irf;
    set(I2,'nLags',nLags, 'nSides', 2,'domainIncr',0.001); % Set number of lags and Sides in IRF
    SRS_inverse{1,2} = I2;
    
    SRS_inverse = nlident(SRS_inverse,Zcur_simulated);
    
elseif SRS_inverse_Weiner == true

    SRS_inverse = lnbl; %Wiener
    set(SRS_inverse,'idMethod','hk');
    I1 = irf;
    set(I1,'nLags',nLags,'nSides',2,'domainIncr',0.001); % Set Number of lags and Sides in IRF
    SRS_inverse{1,1} = I1;

    SRS_inverse = nlident(SRS_inverse,Zcur_simulated);
    
elseif SRS_inverse_IRF == true

    %Two-Sided IRF
    SRS_inverse = irf(Zcur_simulated,'nLags',nLags,'nSides',2,'domainIncr',0.001);
    
end

figure(figNum)
figNum = figNum+1;
[R, V, yp] = nlid_resid(SRS_inverse,Zcur_simulated);

%% Plot the Input/Output, Inverse SRS, and Residuals in One Figure 
%(LNL, Hammerstein, Wiener, or IRF)

if SRS_inverse_LNL == true
    
    figure(figNum)
    figNum = figNum+1;

    Zcur_simulated_double = Zcur_simulated.dataSet;

    subplot(3,3,1)
    Zcur_simulated_double = Zcur_simulated.dataSet;
    Zcur_input = Zcur_simulated_double(:,1);
    plot(t_total(1,1000:end),Zcur_input)
    title('Simulated Paralyzed Displacement, Pos_P(t) (Used as Input)')
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Displacement (m)','Fontsize', 10)

    subplot(3,3,3)
    Zcur_output = Zcur_simulated_double(:,2);
    plot(t_total(1,1000:end),Zcur_output)
    title('Simulated Amplitude Modulation, A(t) (Used as Output)')
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)

    subplot(3,3,4)
    plot(SRS_inverse{1,1})

    subplot(3,3,5)
    plot(SRS_inverse{1,2})

    subplot(3,3,6)
    plot(SRS_inverse{1,3})

    subplot(3,3,7)
    pred = double(yp);
    hold on
    plot(t_total(1,1000:end), Zcur_output)
    plot(t_total(1,1000:end),pred);
    hold off
    V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
    title(['Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)', 'Fontsize', 10)
    legend('Observed', 'Predicted', 'Fontsize', 6)

    subplot(3,3,9)
    plot(R)
    title('Residuals','Fontsize', 10)
    xlabel('Time (s)','Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on
    
elseif SRS_inverse_Hamm == true
    
    figure(figNum)
    figNum = figNum+1;

    subplot(3,2,1)
    Zcur_simulated_double = Zcur_simulated.dataSet;
    Zcur_input = Zcur_simulated_double(:,1);
    plot(t_total(1,1:end-999),Zcur_input)
    title('(a) Simulated Paralyzed Displacement, Pos_P(t) (Used as Input)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Displacement (m)','Fontsize', 10)
    grid on

    subplot(3,2,2)
    Zcur_output = Zcur_simulated_double(:,2);
    plot(t_total(1,1:end-999),Zcur_output)
    title('(b) Simulated Amplitude Modulation, A(t) (Used as Output)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on

    subplot(3,2,3)
    plot(SRS_inverse{1,1})
    title('(c) Static Nonlinear Element','Fontsize', 10)
    xlabel('Paralyzed Displacement Input, Pos_P(t) (m)', 'Fontsize', 10)
    ylabel('Transformed Displacement Output (m)','Fontsize', 10)
    grid on

    subplot(3,2,4)
    plot(SRS_inverse{1,2})
    title('(d) Linear Element','Fontsize', 10)
    xlabel('Lags (s)', 'Fontsize', 10)
    ylabel('X1','Fontsize', 10)
    grid on

    subplot(3,2,5)
    pred = double(yp);
    hold on
    plot(t_total(1,1:end-2998),pred(1000:end-1000,:));
    plot(t_total(1,1:end-2998), Zcur_output(1000:end-1000,:))
    hold off
    V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
    title(['(e) Superimposed A(t), VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)', 'Fontsize', 10)
    legend('Predicted', 'Observed', 'Fontsize', 8)
    grid on

    subplot(3,2,6)
    plot(R(1000:end-1000,:))
    title('(f) Residuals','Fontsize', 10)
    xlabel('Time (s)','Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on

    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(SRS_inverse{1,1})
    ax = gca;
    ax.FontSize = 14;
    title('(a) Static Nonlinear Element','Fontsize', 18)
    xlabel('Paralyzed Displacement Input, Pos_P(t) (m)', 'Fontsize', 18)
    ylabel('Transformed Displacement Output (m)','Fontsize', 18)
    grid on

    subplot(1,2,2)
    plot(SRS_inverse{1,2})
    ax = gca;
    ax.FontSize = 14;
    title('(b) Linear Element','Fontsize', 18)
    xlabel('Lags (s)', 'Fontsize', 18)
    ylabel('X1','Fontsize', 18)
    grid on
    
elseif SRS_inverse_Weiner == true
    
    figure(figNum)
    figNum = figNum+1;

    subplot(3,2,1)
    Zcur_simulated_double = Zcur_simulated.dataSet;
    Zcur_input = Zcur_simulated_double(:,1);
    plot(t_total(1,1:end-999),Zcur_input)
    title('(a) Simulated Paralyzed Displacement, Pos_P(t) (Used as Input)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Displacement (m)','Fontsize', 10)
    grid on

    subplot(3,2,2)
    Zcur_output = Zcur_simulated_double(:,2);
    plot(t_total(1,1:end-999),Zcur_output)
    title('(b) Simulated Amplitude Modulation, A(t) (Used as Output)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on

    subplot(3,2,3)
    plot(SRS_inverse{1,1})
    title('(c) Linear Element','Fontsize', 10)
    xlabel('Lags (s)', 'Fontsize', 10)
    ylabel('X1','Fontsize', 10)
    grid on

    subplot(3,2,4)
    plot(SRS_inverse{1,2})
    title('(d) Static Nonlinear Element','Fontsize', 10)
    xlabel('Transformed Modulation Input', 'Fontsize', 10)
    ylabel('Amplitude Modulation Output (V)','Fontsize', 10)
    grid on

    subplot(3,2,5)
    pred = double(yp);
    hold on
    plot(t_total(1,1:end-2998),pred(1000:end-1000,:));
    plot(t_total(1,1:end-2998), Zcur_output(1000:end-1000,:))
    hold off
    V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
    title(['(e) Superimposed A(t), VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)', 'Fontsize', 10)
    legend('Predicted', 'Observed', 'Fontsize', 8)
    grid on

    subplot(3,2,6)
    plot(R(1000:end-1000,:))
    title('(f) Residuals','Fontsize', 10)
    xlabel('Time (s)','Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on

    figure(figNum)
    figNum = figNum+1;
    subplot(1,2,1)
    plot(SRS_inverse{1,1})
    ax = gca;
    ax.FontSize = 14;
    title('(a) Linear Element','Fontsize', 18)
    xlabel('Lags (s)', 'Fontsize', 18)
    ylabel('X1','Fontsize', 18)
    grid on

    subplot(1,2,2)
    plot(SRS_inverse{1,2})
    ax = gca;
    ax.FontSize = 14;
    title('(b) Static Nonlinear Element','Fontsize', 18)
    xlabel('Transformed Modulation Input', 'Fontsize', 18)
    ylabel('Amplitude Modulation Output (V)','Fontsize', 18)
    grid on

elseif SRS_inverse_IRF == true

    figure(figNum)
    figNum = figNum+1;

    subplot(3,2,1)
    Zcur_simulated_double = Zcur_simulated.dataSet;
    Zcur_input = Zcur_simulated_double(:,1);
    plot(t_total(1,1:end-999),Zcur_input)
    title('(a) Simulated Paralyzed Displacement, Pos_P(t) (Used as Input)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Displacement (m)','Fontsize', 10)
    grid on

    subplot(3,2,2)
    Zcur_output = Zcur_simulated_double(:,2);
    plot(t_total(1,1:end-999),Zcur_output)
    title('(b) Simulated Amplitude Modulation, A(t) (Used as Output)','Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on

    subplot(3,2,[3 4])
    plot(SRS_inverse)
    title('(c) Two-Sided Linear IRF','Fontsize', 10)
    xlabel('Lags (s)', 'Fontsize', 10)
    ylabel('X1','Fontsize', 10)
    grid on

    subplot(3,2,5)
    pred = double(yp);
    hold on
    plot(t_total(1,1:end-999),pred);
    plot(t_total(1,1:end-999), Zcur_output)
    hold off
    V = vaf(Zcur_output(15000:end-14999,:),pred(15000:end-14999,:));
    title(['(d) Superimposed, VAF = ' num2str(round(V,1)) '%'], 'Fontsize', 10)
    xlabel('Time (s)', 'Fontsize', 10)
    ylabel('Amplitude (V)', 'Fontsize', 10)
    legend('Predicted', 'Observed', 'Fontsize', 8)
    grid on

    subplot(3,2,6)
    plot(R)
    title('(e) Residuals','Fontsize', 10)
    xlabel('Time (s)','Fontsize', 10)
    ylabel('Amplitude (V)','Fontsize', 10)
    grid on
    
end

%% Analysis of Inverse SRS Identification Residuals

figure(figNum)
figNum = figNum+1;
sgtitle('Analysis of Inverse SRS Identification Residuals','Fontsize',14)

subplot(2,2,[1 2])
R_temp = double(R);
plot(R(1000:end-1000,:))
ax = gca;
ax.FontSize = 13;
title('(a) Residuals','Fontsize',18)
ylabel('Amplitude (V)','Fontsize',18)
xlabel('Time (s)','Fontsize',18)
grid on

subplot(2,2,3)
histogram(R_temp(1000:end-1000,:))
ax = gca;
ax.FontSize = 13;
title('(b) Residual Distribution','Fontsize',18)
xlabel('Amplitude (V)','Fontsize',18)
ylabel('Density','Fontsize',18)
grid on

S = spect(R);
subplot(2,2,4)
S_frequency = 0.0556:0.0556:0.0556*length(S);
subplot(2,2,4)
semilogy(S_frequency(:,1:500),S(1:500,:),'LineWidth',1.5);
ax = gca;
ax.FontSize = 13;
title('(c) Residual Power Spectrum','Fontsize',18);
ylabel('PSD (log)','Fontsize',18); 
xlabel('Frequency (Hz)','Fontsize',18);
grid on

%% Inverse SRS Validation

control_system_test_output_double = [];
inverse_SRS_validation_accuracy = [];

for signal = 1:num_signals-1
    
    control_system_test_output = nlsim(SRS_inverse,validation_output(:,signal));
    set(control_system_test_output,'domainIncr',0.001);
    control_system_test_output_double(:,signal) = control_system_test_output.dataSet;
    
    inverse_SRS_validation_accuracy(signal,:) = vaf(simulated_input(1999:end-1000,signal+1),control_system_test_output_double(1000:end-1000,signal));

end

% Mean and Standard Deviation of Inverse SRS Validation Accuracy
validation_accuracy_mean = mean(inverse_SRS_validation_accuracy);
validation_accuracy_std = std(inverse_SRS_validation_accuracy);

disp(['Validation Accuracy Mean: ' num2str(round(validation_accuracy_mean,1)) '%'])
disp(['Validation Accuracy Std: ' num2str(round(validation_accuracy_std,1)) '%'])

%Plot one validation example
figure(figNum)
figNum = figNum+1;
hold on
plot(t_total(1,1:end-2998),control_system_test_output_double(1000:end-1000,1))
plot(t_total(1,1:end-2998),simulated_input(1999:end-1000,2))
hold off
title(['Superimposed Amplitude Modulation Signals of Validation Trial, VAF = ' num2str(round(inverse_SRS_validation_accuracy(1,1),1)) '%'], 'Fontsize', 19)
xlabel('Time (s)', 'Fontsize', 20)
ylabel('Amplitude Modulation, A(t) (V)', 'Fontsize', 20)
legend('Predicted', 'Observed', 'Fontsize', 16)
grid on

%% Analysis of Inverse SRS Validation Residuals

R_validation = simulated_input(1999:end-1000,2) - control_system_test_output_double(1000:end-1000,1);

figure(figNum)
figNum = figNum+1;
sgtitle('Analysis of Inverse SRS Residuals from Validation Trial','Fontsize',14)

subplot(2,2,[1 2])
plot(t_total(1,1:end-2998),R_validation)
ax = gca;
ax.FontSize = 13;
title('(a) Residuals of Validation Trial','Fontsize',18)
ylabel('Amplitude (V)','Fontsize',18)
xlabel('Time (s)','Fontsize',18)
grid on

subplot(2,2,3)
histogram(R_validation(1000:end-1000,:))
ax = gca;
ax.FontSize = 13;
title('(b) Residual Distribution','Fontsize',18)
xlabel('Amplitude (V)','Fontsize',18)
ylabel('Density','Fontsize',18)
grid on

S = spect(R_validation);
subplot(2,2,4)
S_frequency = 0.0556:0.0556:0.0556*length(S);
subplot(2,2,4)
semilogy(S_frequency(:,1:500),S(1:500,:),'LineWidth',1.5);
ax = gca;
ax.FontSize = 13;
title('(c) Residual Power Spectrum','Fontsize',18);
ylabel('PSD (log)','Fontsize',18); 
xlabel('Frequency (Hz)','Fontsize',18);
grid on

tEnd = toc(tStart)/60
%% Control system: Identify the Inverse SRS

%Identify the Inverse SRS Model:
%   1. Simulate the response of the SRS you identify to a white input to generate a simulated output.
%   2. Then identify a model between the simulated output, as an input, and the simulated input, to estimate the inverse model.
%   3. By using the simulated signals you can avoid problems with noise.

tStart = tic;

%% Set Initial Parameters
time = 180;
set_output_noise_power = 0;
figNum = 300;
Fs = 1000;

%Select which Model to Use
if compare_two_models == true && LNL_model == true
    SRS_model = LNL2;
elseif compare_two_models == true && Hammerstein_model == true
    SRS_model = Hammerstein2;
elseif compare_two_models == true && Weiner_model == true
    SRS_model = Weiner2;
elseif compare_two_models == true && Linear_IRF_model == true
    SRS_model = IRF2;
elseif compare_two_models == false && LNL_model == true
    SRS_model = LNL;
elseif compare_two_models == false && Weiner_model == true
    SRS_model = Weiner;
else
    SRS_model = IRF_model;
end

%% Simulated Input (PRBS Input)

simulated_input = [];

num_signals = 2;

PRBS_stimulus_time = 480;
variable_amplitude = true;
N = PRBS_stimulus_time/10;
M = 10000;
PRBS_amplitude = 20; %mm

for signal = 1:num_signals

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

    U = (u.InputData)';
    desired_displacement = A.*U;

    stim_frequency = 50;
    stim_amplitude = desired_displacement*170;

    input_stimulus = max(stim_amplitude.*sin(2*pi*stim_frequency.*t_total),0);

    simulated_input(:,signal) = stim_amplitude;

end
    
%% Simulate the reponse of the SRS model to the input

simulated_output = [];
simulated_outputs = [];

for signal = 1:num_signals

    simulated_output = nlsim(SRS_model,simulated_input(:,signal));
    set(simulated_output, 'domainIncr',0.001)
    
    simulated_output = max(simulated_output,0);

    simulated_outputs(:,signal) = simulated_output.dataSet;

end

%%
% Identify a model between the simulated output (as an input) and the
% simulated input (as an output)

Zcur_simulated = [simulated_outputs(1000:end,1) simulated_input(1000:end,1)];
Zcur_simulated = nldat(Zcur_simulated,'domainIncr',0.001,'comment','Simulated Output (as Input), Simulated Input (as Output)','chanNames', {'Simulated Output (as Input)' 'Simulated Input (as Output)'});

figure(figNum)
figNum = figNum+1;
subplot(2,1,1)
plot(t_total(1,1:end-999),simulated_input(1000:end,1))
ax = gca;
ax.FontSize = 15;
title('Simulated Amplitude Modulation','Fontsize',20)
ylabel('Amplitude (V)','Fontsize',18)
grid on

subplot(2,1,2)
plot(t_total(1,1:end-999),simulated_outputs(1000:end,1))
ax = gca;
ax.FontSize = 15;
title('Simulated Displacement','Fontsize',20)
ylabel('Displacement (m)','Fontsize',18)
xlabel('Time (s)','Fontsize',18)
grid on

test_output = [];

for signal = 2:num_signals

    Test_Output = nldat(simulated_outputs(1000:end,signal));
    set(Test_Output,'domainIncr',0.001);
    
    test_output(:,signal-1) = Test_Output.dataSet;

end

%%

nLags = 400;
nLags_plot = [];

accuracy_identification = [];
accuracy_validation = [];

while nLags > 50

    SRS_inverse = lnlbl;  %LNL Model
    set(SRS_inverse,'idMethod','hk','hkTolerance', 0.5,...
        'nhkMaxIts', 3, 'nhkMaxInner', 3);
    
    I1 = irf;
    set(I1,'nLags',1,'nSides',2,'domainIncr',0.001);
    % SRS_inverse{1,1} = I1;
    P = polynom;
    I3 = irf;
    set(I3,'nLags',nLags,'nSides', 2,'domainIncr',0.001);
    % SRS_inverse{1,3} = I3;
    SRS_inverse.elements = {I1 P I3};


%     SRS_inverse = nlbl;  %Hammerstein
%     set(SRS_inverse,'idMethod','hk','displayFlag',true,'threshNSE',.001);
%     I2 = irf;
%     set(I2,'nLags',nLags, 'nSides', 2); % Set number of lags and Sides in IRF
%     SRS_inverse{1,2} = I2;

%     SRS_inverse = lnbl; %Wiener
%     set(SRS_inverse,'idMethod','hk');
%     I1 = irf;
%     set(I1,'nLags',nLags,'nSides',2); % Set Number of lags and Sides in IRF
%     SRS_inverse{1,1} = I1;

    SRS_inverse = nlident(SRS_inverse,Zcur_simulated);

    %Two-Sided IRF
%     SRS_inverse = irf(Zcur_simulated,'nLags',nLags,'nSides',2);

    figure(figNum)
    figNum = figNum+1;
    [R, V, yp] = nlid_resid(SRS_inverse,Zcur_simulated);
    
    accuracy_identification = [accuracy_identification V];

    %% Inverse SRS Validation

    control_system_test_output_double = [];
    inverse_SRS_validation_accuracy = [];

    for signal = 1:num_signals-1

        control_system_test_output = nlsim(SRS_inverse,test_output(:,signal));

        control_system_test_output_double(:,signal) = control_system_test_output.dataSet;

        inverse_SRS_validation_accuracy(signal,:) = vaf(simulated_input(1999:end-1000,signal+1),control_system_test_output_double(1000:end-1000,signal));

    end
    
    accuracy_validation = [accuracy_validation inverse_SRS_validation_accuracy(1,:)];

    nLags_plot = [nLags_plot nLags];
    nLags = nLags - 100;

end

%% Plot the Results
figure(figNum)
figNum = figNum+1;
hold on
plot(nLags_plot,accuracy_identification,'Linewidth',2)
plot(nLags_plot,accuracy_validation,'Linewidth',2)
hold off
title('Accuracy vs NLags','Fontsize',20)
xlabel('NLags', 'Fontsize',18)
ylabel('Accuracy (%VAF)','Fontsize',18)
legend('Identification','Validation','Fontsize',16)
grid on

tEnd = toc(tStart)/60
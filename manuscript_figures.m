% Journal Manuscript Figures

%% ERS Identification
clear all
load('ERS_ident')

figNum = 1;

NHK1 = NHK_all(1);
NHK2 = NHK_all(2);

emg_pdf1 = emg_all(:,1);
emg_pdf2 = emg_all(:,2);
emg_pdf1(emg_pdf1==0) = []; %Removing the zero values
emg_pdf2(emg_pdf2==0) = []; %Removing the zero values

x_lim = 0.0008;

%Sets the Upper Limit and Lower limit for "High Quality" values of the
%model nonlinear element
upper_limit1 = 0.000586;
lower_limit1 = -0.000613;
upper_limit2 = 0.000602;
lower_limit2 = -0.00063;
upper_limit3 = max(upper_limit1, upper_limit2);
lower_limit3 = min(lower_limit1, lower_limit2);

% figure(figNum)
% figNum = figNum+1;
% subplot(2,1,1)
% plot(NHK1{1,1})
% ax = gca;
% ax.FontSize = 15;
% title('(a) Static Nonlinearity','Fontsize', 24)
% xlabel('EMG Input, E(t) (V)','Fontsize',18)
% ylabel('Transformed EMG (V)','Fontsize',18)
% xline(upper_limit1,'--r','LineWidth',2);
% xline(lower_limit1,'--r','LineWidth',2);
% xlim([-x_lim x_lim])
% grid on
% 
% subplot(2,1,2)
% histogram(emg_pdf1,'Normalization','pdf')
% ax = gca;
% ax.FontSize = 15;
% xline(upper_limit1,'--r','LineWidth',2);
% xline(lower_limit1,'--r','LineWidth',2);
% xlabel('EMG Amplitude (V)','Fontsize',18)
% ylabel('Density','Fontsize',18)
% title('(b) EMG Distribution','Fontsize',24)
% xlim([-x_lim x_lim])
% grid on

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

%% SRS Identification
load('SRS_ident')


%% ERS/SRS Noise
clear all
load('ERS_noise')

figNum = 400;

% accuracy_identification_all = max(0, accuracy_identification_all);
% 
% figure(figNum)
% figNum = figNum+1;
% plot(noise_snr_all_clean(2,:),accuracy_identification_all(2,:),'LineWidth',3);
% hold on
% plot(noise_snr_all_clean(1,:),accuracy_identification_all(1,:),'LineWidth',3);
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

accuracy_validation = max(0, accuracy_validation);
accuracy_validation = accuracy_validation(:,[1:8 10:end],:);
accuracy_validation_mean = mean(accuracy_validation);
accuracy_validation_std = std(accuracy_validation);

noise_snr_all_clean = noise_snr_all_clean(1, [1:8 10:end]);

figure(figNum)
figNum = figNum+1;
hold on

%plot(noise_snr_all_clean(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
%plot(noise_snr_all_clean(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
%patch([noise_snr_all_clean(2,:) fliplr(noise_snr_all_clean(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(noise_snr_all_clean(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(1,:) fliplr(noise_snr_all_clean(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'b','FaceAlpha',0.2)

%plot(noise_snr_all_clean(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
x1 = plot(noise_snr_all_clean(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3, 'Color', [0 0.4470 0.7410]);

load('SRS_noise')

%accuracy_identification_all = max(0, accuracy_identification_all);

% figure(figNum)
% figNum = figNum+1;
% plot(noise_snr_all_clean(2,:),accuracy_identification_all(2,:),'LineWidth',3);
% hold on
% plot(noise_snr_all_clean(1,:),accuracy_identification_all(1,:),'LineWidth',3);
% hold off
% ax = gca;
% ax.FontSize = 15;
% title('Identification Accuracy vs Output Noise','Fontsize',24);
% legend('PRBS','Physiological','Location','southeast','FontSize',18)
% ylabel('Accuracy (% VAF)','Fontsize',18); 
% xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
% grid on

accuracy_validation = max(0, accuracy_validation);
accuracy_validation_mean = mean(accuracy_validation);
accuracy_validation_std = std(accuracy_validation);

%plot(noise_snr_all_clean(2,:),min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100),'LineStyle','none','LineWidth',2)
%plot(noise_snr_all_clean(2,:),max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0),'LineStyle','none','LineWidth',2)
%patch([noise_snr_all_clean(2,:) fliplr(noise_snr_all_clean(2,:))], [min(accuracy_validation_mean(:,:,2)+accuracy_validation_std(:,:,2),100) fliplr(max(accuracy_validation_mean(:,:,2)-accuracy_validation_std(:,:,2),0))], 'b','FaceAlpha',0.2)

plot(noise_snr_all_clean(1,:),min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100),'LineStyle','none','LineWidth',2)
plot(noise_snr_all_clean(1,:),max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0),'LineStyle','none','LineWidth',2)
patch([noise_snr_all_clean(1,:) fliplr(noise_snr_all_clean(1,:))], [min(accuracy_validation_mean(:,:,1)+accuracy_validation_std(:,:,1),100) fliplr(max(accuracy_validation_mean(:,:,1)-accuracy_validation_std(:,:,1),0))], 'y','FaceAlpha',0.2)

%plot(noise_snr_all_clean(2,:),accuracy_validation_mean(:,:,2),'LineWidth',3);
x2 = plot(noise_snr_all_clean(1,:),accuracy_validation_mean(:,:,1),'LineWidth',3,'Color',[0.8500 0.3250 0.0980]);

hold off
ax = gca;
ax.FontSize = 15;
title('Validation Accuracy vs Output Noise','Fontsize',24);
legend([x1 x2], 'ERS','SRS','Location','southeast','FontSize',18)
ylabel('Accuracy (% VAF)','Fontsize',18); 
xlabel('Signal to Noise Ratio (dB)','Fontsize',18);
grid on


%% ERS/SRS Record Length
load('ERS_record_length')

validation_accuracy = max(validation_accuracy, 0);
validation_accuracy_mean = mean(validation_accuracy);
validation_accuracy_var = var(validation_accuracy);
validation_accuracy_std = std(validation_accuracy);

figure(figNum)
figNum = figNum+1;
hold on
plot(PRBS_movement_time,min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
plot(PRBS_movement_time,max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
patch([PRBS_movement_time fliplr(PRBS_movement_time)], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'b','FaceAlpha',0.2)

x1 = plot(PRBS_movement_time,validation_accuracy_mean,'LineWidth',3, 'Color', [0 0.4470 0.7410]);

load('SRS_record_length')

validation_accuracy = max(validation_accuracy, 0);
validation_accuracy_mean = mean(validation_accuracy);
validation_accuracy_var = var(validation_accuracy);
validation_accuracy_std = std(validation_accuracy);

plot(PRBS_stimulus_time,min(validation_accuracy_mean+validation_accuracy_std,100),'LineStyle','none','LineWidth',2)
plot(PRBS_stimulus_time,max(validation_accuracy_mean-validation_accuracy_std,0),'LineStyle','none','LineWidth',2)
patch([PRBS_stimulus_time fliplr(PRBS_stimulus_time)], [min(validation_accuracy_mean+validation_accuracy_std,100) fliplr(max(validation_accuracy_mean-validation_accuracy_std,0))], 'y','FaceAlpha',0.2)

x2 = plot(PRBS_stimulus_time,validation_accuracy_mean,'LineWidth',3,'Color',[0.8500 0.3250 0.0980]);

hold off
ax = gca;
ax.FontSize = 14;
title('Identification Signal Record Length vs Validation Accuracy','Fontsize', 24) 
xlabel('Record Length (s)','Fontsize', 18)
ylabel('Accuracy (%VAF)','Fontsize', 18)
legend([x1 x2], 'ERS', 'SRS')
grid on

%% FRCS
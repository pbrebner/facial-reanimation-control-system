% Journal Manuscript Figures

%% ERS Identification


%% SRS Identification


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
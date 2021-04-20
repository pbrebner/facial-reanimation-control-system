% This function simulates glucose levels following a meal ingestion at time 30 min.
% This function allows delivering a meal bolus at time 30 min.
% This function allows the implementation of closed-loop controllers to deliver variable basal rates.

%% clear
clc
clear all;

%% duration of simulation and sampling frequency
time = 0:10:480;

%% initial conditions (do not modify)
Gs0 = 5.5;

k12 = 8.6*10^(-2); %(l/min)
kis = .02; %(l/min)
ke = .076; %(l/min)
ci = .10/60; %(U/min)
w = 45; %(kg)
Vi = 190; %(ml/kg)
V = 160; %(ml/kg)
st = 11.4*10^(-4); %(/min/mU/l)
sd = 3.36*10^(-4); %(/min/mU/l)
se = 117*10^(-4); %(/mU/l)

u_ss = 1.6/60; %(U/min)
Qis1 = u_ss/kis; %(U)
Qis2 = Qis1; %(U)
Qi0 = (Qis2*kis + ci)/ke; %(U)
Ip = Qi0/Vi/w*1e6; %(mU/l)

x10 = Ip; %(mU/l)
x20 = Ip; %(mU/l)
x30 = Ip; %(mU/l)

Q10 = Gs0*V; %(umol/kg)
Q20 = Q10*x10*st/(x20*sd+k12); %(umol/kg)
y0 = [Qis1 Qis2 Qi0 x10 x20 x30 Q10 Q20 Gs0];

%% Simulation of meal ingestion at t = 30 min
T_total = [];
Y_total = [];
u_total = [];
Glucose = Gs0; %this will be an array that includes the glucose level values 
               %Glucose(1) will be the glucose value at time 0, Glucose(2)
               %will be the glucose value at time 10, etc.

for i = 1:length(time)-1 %for loop from time 0 to time 480 (8 hours)
    %Meal ingestion at t = 30 min
    if time(i) == 30
        %This simulates a meal that is ingested at time 30 min. 
        %The default is a meal of size 60 grams of carbohydrate. You can modify the meal size as needed. 
        meal = 60; % modify this to determine the meal size
    else
        meal = 0;
    end

    %insulin bolus at t = 30 min
    if time(i) == 30
        %This simulates an insulin bolus that is ingested at time 30 min.
        %The default is 0 units. You can modify the bolus size as needed. 
        bolus = 6; % modify this to determine the bolus size
        y0(1) = y0(1) + bolus;
    end

    %Implement your controller here as a function of the glucose level signal (Glucose).
    %Your controller should be added to the steady state basal rate: u_ss.
    %Your controller should replace the 0 below.
    basal_insulin = u_ss + 0; 

    %Remember: insulin delivery can not be negative.
    if basal_insulin < 0
        error('Insulin delivery can not be negative')
    end
    
    %This solves for the dynamics of the glucose system
    [T,Y] = ode45(@(t,y) model(t, y, basal_insulin, meal), [time(i) time(i+1)], y0);
    G = Y(end,9);
    Glucose = [Glucose G]; %this is glucose level signal. Use this array in your controller.

    %initial condition for next cycle
    y0 = Y(end,:);
    
    %for plotting
    T_total = [T_total ; T];
    Y_total = [Y_total ; Y];
    u_total = [u_total ; basal_insulin*60]; %60 to change from /min to /hr
end

%plotting
subplot(2,1,1);
plot(T_total, Y_total(:,9))
hold on;
plot(time,4*ones(size(time)),'--r','LineWidth',0.01);
plot(time,10*ones(size(time)),'--r','LineWidth',0.01);
title('Closed-loop simulation. Meal ingested at time 30 min.');
xlabel('Time (min)');
ylabel('Glucose Levels (mmol/L)');
axis([0 time(end) 0 20]);

subplot(2,1,2);
stairs(time(1:end-1), u_total)
xlabel('Time (min)');
ylabel('Insulin Delivery (U/min)');
axis([0 time(end) 0 10]);


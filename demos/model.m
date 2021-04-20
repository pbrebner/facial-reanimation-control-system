function dydt = model(t, y, basal_insulin, meal)

%This function solves the dynamics of the glucose system. DO NOT MODIFY.

global meal_time meal_value meal_added

%model parameters
kis = .02; %(l/min)
ke = .076; %(l/min)
ci = .10/60; %(U/min)
w = 45; %(kg)
ka1 = 4*10^(-2); %(l/min)
ka2 = 5.9*10^(-2); %(l/min)
ka3 = 6.3*10^(-2); %(l/min)
st = 11.4*10^(-4); %(/min/mU/l)
sd = 3.36*10^(-4); %(/min/mU/l)
se = 117*10^(-4); %(/mU/l)
f01 = 7.3; %(umol/kg/min)
k12 = 8.6*10^(-2); %(l/min)
EGP0 = 26.3; %(umol/kg/min)
d = 9.6;%(min)
km = .025; %(/min)
pm = .8 ; %(unitess)
k = 0.0614;
Vi = 190; %(ml/kg)
V = 160; %(ml/kg)

Qis1 = y(1);
Qis2 = y(2);
Qi = y(3);
x1 = y(4);
x2 = y(5);
x3 = y(6);
Q1 = y(7);
Q2 = y(8);
Gs = y(9);

dydt = zeros(9,1);
if meal ~= 0 && isempty(meal_added)
    meal_value = meal;
    meal_time = t;
    meal_added = 1;
end

%% subcutaneous insulin absorption subsystem
dydt(1) = basal_insulin - Qis1 * kis;
dydt(2) = Qis1 * kis - Qis2 * kis;

%% plasma insulin kinetics subsystem
dydt(3) = Qis2 * kis - Qi * ke + ci;

Ip = Qi * 10^6 / (Vi * w); %(mU/l)

%% insulin action subsystem
dydt(4) = -ka1 * x1 + ka1 * Ip;
dydt(5) = -ka2 * x2 + ka2 * Ip;
dydt(6) = -ka3 * x3 + ka3 * Ip;

%% gut absorption subsystem
Um = 0;
if t > meal_time 
    Um1 = km^2 * (t-meal_time) * exp(-km*(t-meal_time))*(meal_value*5551)/w * pm;
    if t > meal_time + d
       Um2 = km^2 * (t-meal_time-d) * exp(-km*(t-meal_time-d))*(meal_value*5551)/w * (1-pm);
    else
       Um2 = 0;           
    end
    Um = (Um1 + Um2);
end

%% glucose kinetics subsystem
if x3*se < 1
   dydt(7) = -f01 * (Q1/160)/(1+Q1/160) - x1*st*Q1 + k12*Q2 + EGP0*(1-x3*se) + Um;
   dydt(8) = x1*st * Q1 - (k12+x2*sd) * Q2;
else
   dydt(7) = -f01 * (Q1/160)/(1+Q1/160) - x1*st*Q1 + k12*Q2 + Um;
   dydt(8) = x1*st * Q1 - (k12+x2*sd) * Q2;
end

Gt = Q1 / V; %(mmol/l)

%% glucose sensor
dydt(9) = k*Gt - k*Gs;


clc
clear all

global u y bolus t loc_num 

%input
u = [0.0 0.2 1.3	1.7	1.1	0.5	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.1	1.7	0.3	0.3	1.8	2.2	1.9	0.3	0.8	0.1	0.5	1.9	1.9	2.0	2.0	2.0	2.1	0.5	0.0	0.1	0.1	0.3	0.7	1.2	1.6	1.6	0.9	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.1	0.8	1.6	1.6	1.6	1.6	1.6	1.6	0.7	0.0	0.0	0.0	0.0];

%output
y = [27.31	nan	nan	26.39	nan	nan	26.99	nan	nan	19.35	nan	26.03	25.18	25.94	26.04	nan	18.84	nan	16.63	nan	13.31	nan	28.95	nan	32.75	nan	35.60	nan	35.07	nan	33.07	nan	31.31	nan	31.84	nan	31.71	nan	35.16	nan	32.61	nan	30.60	nan	31.16	nan	31.81	nan	30.98	nan	28.26	nan	25.31	nan	23.41	nan	24.61	nan	25.12	nan	23.36	nan	21.64	nan	14.84	nan	12.70	nan	16.17	nan	12.83	nan	13.70	nan	10.52	nan	8.89	nan	9.12	nan	11.03	nan	13.45	nan	19.18	nan	15.65	nan	15.79];
bolus = 7.1; %@19h20
t = 0:10:(length(u)-1)*10;

stairs(t,u);%plots a stair plot
%%
hold;

loc_num = ~isnan(y);
plot(t(loc_num),y(loc_num),'x')
ylabel('y');
xlabel('t');

%%
%least square estimate
x = lsqnonlin(@myfun,[1/50 10 15 1]); %you give it an array of least square errors (number array)
%x(1) = k, x(2) = C, x(3) = Ib; x(4) = x10;

%plotting the fit
global k i
k = x(1);
x10 = x(4);
y0 = [x10 x10]; %assume steady state
Y_total = [];
T_total = [];

for i=1:length(t)-1
    [T,Y] = ode45(@(t,y) odefcn(y), t(i:i+1), y0);
    y0 = Y(end,:);
    if i == 20
        y0(1) = y0(1)+bolus;
    end
    Y_total = [Y_total; Y];
    T_total = [T_total; T];
end
plot(T_total,Y_total(:,1))
plot(T_total,Y_total(:,2)*x(2)+x(3))

legend('input', 'data', 'Q_1', 'output')
%plot(t(loc_num),myfun(x))

function F = myfun(x)
global u y bolus t i k
k = x(1);
x10 = x(4); %x(3) + x(2)*u(1)/x(1);
y0 = [x10 x10]; %assume steady state
Y_total = y0;
F = y0(2)*x(2)+x(3) - y(1); %Cost Function (Estimated Ip minus Acutal Ip from Data)
T_total = [];
for i=1:length(t)-1
    [T,Y] = ode45(@(t,y) odefcn(y), t(i:i+1), y0);
    y0 = Y(end,:);
    if i == 20
        y0(1) = y0(1)+bolus;
    end

    Y_total = [Y_total; Y(end,:)];
    T_total = [T_total; T(end)];
    if ~isnan(y(i+1)) %Only for values of y that are not NAN
        F = [F Y_total(i+1,2)*x(2)+x(3)-y(i+1)]; %Cost Function
    end
end

end

function dydt = odefcn(y)
global u i k

dydt(1) = u(i)/60 - k*y(1);
dydt(2) = k*y(1) - k*y(2);

dydt = dydt';

end


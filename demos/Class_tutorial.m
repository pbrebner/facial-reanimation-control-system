% BMDE 502 - February 5th Class
clear all

y0 = [1 0];

u = [1, 2, 0, 1];
t = 0:10:30;
T_total = [];
Y_total = [];

options = odeset('MaxStep', 1e-20); %This step size takes the function a lomg time

for i = 1:length(u)-1
    
    [T, Y] = ode45(@(t,y) model(t, y, u(i)), [t(i) t(i+1)], y0);
    y0 = Y(end,:);
    T_total = [T_total; T];
    Y_total = [Y_total; Y];
    
end

%plotting
plot(T_total, Y_total)
grid

function dydt = model(t, y, u)

dydt = zeros(2,1);

dydt(1) = t - 1/70*y(1);
dydt(2) = -y(2)*1/70 + y(1)*1/70;

end
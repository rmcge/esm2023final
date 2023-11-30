% Set values
dx = 1;
x = 0:dx:50;
B = 50;
dt = 150;
t = 0:dt:5400;
g = 9.8;
S0 = 0.0002;
n = 0.01;
V = 4;
C = 1;

% Building
Q = zeros(length(t), length(x));

% Estimated Values
Q0 = 510;
amp = 1.5;
T = 10;

% Average depth of the Nile River is 8-11m
y = 9.5;

% Assuming the river is a rectangle for cross-sectional area
A = (B * y) * ones(size(x));

% Setting initial alpha and beta
alpha_initial = 5;
beta_initial = 15;

% Initialize Q at the upstream boundary
Q(:, 1) = Q0 * (1 + amp * sin((2 * pi * t) / T));

% Dam break perturbation
t_dam_break = 2000;  % Adjust the time when the dam breaks
for j = 1:length(t)-1
    for i = 2:length(x)-1
        R = A(i)/((2*y)+B) % Guestimate for perimeter
        Sf = sqrt((V * n) / R^(2/3));

        alpha = 2 * (Q(j, i) / A(i)) + (((g * A(i) / B) - (Q(j, i)^2 / A(i)^2)) / ((Q(j, i) / A(i)) * (5/3 - (4 * R / (3 * B)))));
        beta = g * A(i) * (Sf - S0);

        % Update equation for Q
        Q(j+1, i) = Q(j, i) + alpha * (dt / dx) * (Q(j, i) - Q(j, i-1)) - beta * dt;
    end

    % Inflow boundary condition
    Q(j+1, 1) = Q0 * (1 + amp * sin((2 * pi * t(j+1)) / T));

    % Reflective boundary condition at the downstream end
    Q(j+1, end) = Q(j+1, end-1);

    % Dam break perturbation
    if t(j) >= t_dam_break && t(j) < t_dam_break + dt
        Q(j+1, 1) = Q(j+1, 1) + 400 * exp(-(t(j) - t_dam_break)^2 / (2 * 100^2));
    end
end

% Plotting
figure;
plot(t, Q);
title('Channel Flow Over Time with Dam Break');
xlabel('Time (s)');
ylabel('Channel Flow (m^3/s)');

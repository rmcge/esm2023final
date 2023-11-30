% Set values
dx = 1;
x = 0:dx:50;
B = 1200;
dt = 150;
t = 0:dt:54000;
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
A = B * y;

% Setting initial alpha and beta
alpha = 5;
beta = 15;

% Initial condition for Q
Q(:, 1) = Q0 * (1 + amp * sin((2 * pi * t) / T));

for j = 1:length(t)-1
    for i = 1:length(x)-1
        % Parameters for the update equation
        R = A / 58; % Guestimate for perimeter
        Sf = sqrt((V * n) / R^(2/3));

        % Calculate alpha and beta
        alpha = 2 * (Q(j, i) / A) + (((g * A / B) - (Q(j, i)^2 / A^2)) / ((Q(j, i) / A) * (5/3 - (4 * R / 3 * B))));
        beta = g * A * (Sf - S0);

        % Update equation for Q
        Q(j+1, i+1) = Q(j, i) + alpha * Q(j, i) * (dt / dx) + beta;

        % Handle boundary conditions (you may need to adjust this based on your problem)
        Q(j+1, 1) = Q(j+1, 2);
    end
end

% Plotting
figure;
plot(t, Q(:, 1));
title('Elevation at x = 0');
xlabel('Time (s)');
ylabel('Elevation (m)');

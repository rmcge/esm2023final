%Given in Keskin paper:
L = 2000; %length of channel in meters
b = 5; %width of the bottom of the channel in meters
S0 = 0.0005; %slope of the bottom of the channel
n = 0.0138; %Manning's roughness coefficient for the bottom of the channel
g = 9.8; %acceleration due to gravity, m^2/s
Q0 = 3; %base flow of inflow hydrograph in m^3/s
Qp = 12; %peak flow of inflow hydrograph in m^3/s
tp = 11;
tb = 21;

%more constants that we need but have to make up
h = 1.5; %water depth

%set up time and space vectors
dt = 10; %seconds
t = 0:dt:3000; %time vector
dx = 100; %meters
x = 0:dx:L; %space vector

%pre-allocate Q matrix and A matrix
Q = zeros(length(x), length(t));
Q(:, 1) = Q0;
for z = 0:tp
    Q(1,z+1) = Q0 + ((Qp-Q0)/(tp-1))*z;
end
for z = 0:tp
    Q(1,z+10) = Qp - ((Qp-Q0)/(tb-tp))*z;
end

Q(1,tb:end) = Q0;
A = zeros(length(x), length(t));
A(:,1) = b*h;
R = zeros(length(x), length(t));
R(1,1) = A(1,1)/(2*h+b);

%set initial and boundary conditions
tp = 10;
tb = 20;

for i = 2:length(x)-1
    for j = 1:length(t)-1
        R(i,j) = A(i,j)/(2*h+b);

        alpha = 2*(Q(i,j)/A(i,j)) + (((g*A(i,j))/b)+(Q(i,j)/A(i,j))^2)/((Q(i,j)/A(i,j))*((5/3)-(4/3)*(R(i,j)/b)));
        beta = g*A(i,j)*(((Q(i,j)^2*n^2)/(A(i,j)*R(i,j)^(4/3)))-S0);

        Q(i, j+1) = Q(i,j) - (dt/dx) * alpha * (Q(i,j)-Q(i-1,j)) + beta*dt;
        A(i, j+1) = A(i,j) - (dt/dx) * (Q(i, j+1)-Q(i-1,j+1));
    end
end

figure;
hold on
plot(t,Q(1,:), 'b');
% plot(t,Q(7,:), 'g');
% plot(t,Q(13,:), 'm');
% plot(t,Q(21, :), 'c');
ylim([0, inf]);
xlabel('Time (s)');
ylabel('discharge (m^3/s')
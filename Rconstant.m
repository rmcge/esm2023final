% Set values
dx = 100;
x = 0:dx:2000; %channel length is 2000m
B = 5; %channel width is 5m
dt = 1;
t = 0:dt:50; %50 seconds
g = 9.8;
S0 = 0.0005; %slope of channel bed
n = 0.0138; %roughness coefficient 
V = 1.5; %not sure about this velocity completely made up number
h = 3.5; %completely made up for water level
C = 0.5;

% Building
Q = zeros(length(x), length(t));
A = ones(length(x), length(t));

%insert boundary conditions?
Q0 = 3;
Qp = 12;
Qb = 3;
tp = 11;
tb = 21;
Q(:,1) = 3;
A(:,1) = Q0/V;
% Q(1,1:tp) = Q0 + ((Qp-Q0)/(tp));
% Q(1,tp:tb) = Q0 - ((Qp-Q0)/(tb-tp));
% Q(1,tb:length(t)) = Q0;
 
R = 0.05


for i = 2:length(x)-1
    for j = 1:length(t)-1 
        %R = A(i,j)/((2*h)+B);
        Sf = sqrt((V * n) / R^(2/3));

        alpha = 2 * (Q(i,j) / A(i,j)) + (((g * A(i,j) / B) - (Q(i,j)^2 / A(i,j)^2)) / ((Q(i,j) / A(i,j)) * (5/3 - (4 * R / (3 * B)))));
        beta = g*A(i)*(Sf - S0);

        Q(i,j+1) = Q(i,j) - (dt/dx)*alpha*(Q(i,j)-Q(i-1,j)) + beta*dt;
        A(i,j+1) = A(i,j) - (dt/dx)*(Q(i,j+1)-Q(i-1,j+1)); %if you take out this line it changes a lot!!!!
    end
    Q(1,1:tp) = Q0 + ((Qp - Q0) / tp) * t(1:tp);
    Q(1,tp:tb) = Q0 - ((Qp - Q0)/(tb - tp)) * (t(tp:tb)-tb);
    Q(1,tb:end) = Q0;
end

figure;
plot(t,Q(1,:));
xlabel('Time (s)');
ylabel('Flow Rate (m/s)');


for y = 2:7
    figure;
    plot(t,Q(y,:));
    xlabel('Time (s)');
    ylabel('Flow Rate (m/s)');
end
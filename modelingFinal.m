%Set values
dx = 1
x = 0:dx:50;
B = 50
dt = 150
t = 0:dt:54000;
g = 9.8
S0 = .0002
n = .01
V = 4
C =1
%Building
% beta = zeros(length(t),length(x));
% alpha = zeros(length(t),length(x));
Q = zeros(length(t),length(x));

% Estimated Values
Q0 = 510
amp = 1.5
T = 10
%average depth of the nile river is 8-11m
y = 9.5
% assuming river is a rectangle for cross sectional area
A = B * y
%Setting initial 
alpha = 5
beta = 15
Q(:,1) = Q0*(1+ amp*sin((2*pi*t)/T))'


for j = 2:length(t)-1
    for i = 2:length(x)-1
    %R is guestimate lol for perimeter
    % R = A/58;
    % Sf = sqrt((V*n)/R^(2/3));
    alpha = 2*(Q(j,i)/A(j,i))+(((g*A(j,i)/B)-(Q(j,i)^2/A(i,j)^2))/((Q(j,i)/A(j,i))*((5/3)-(4*R/3*B))));
    beta = g*A*(Sf-S0);

    Q(j+1,i) = Q(j,i) + alpha*(dt/dx)*(Q(j,i)-Q(j,i-1)) - beta*dt;


    %Q(j,i+1) = Q(j,i) + alpha*Q(j+1,i)*(dt/dx) + beta;
    % alpha(i,j) = 2*Q(i,j)/A +((((g*A)/B)-Q(i,j).^2/A.^2))/((Q(i,j)/A)* (5/3-(4.*R)/3.*B));
    % beta(i,j) = g*A*(Sf-S0);

    % Q(i,j+1) = Q(i,j) - C*(Q(i,j)-Q(i,j-1));
    % Q(i,j+1) = Q(i,j)-alpha(i,j)*(150/1)*(Q(i,j)-Q(i-1,j))-beta(i,j)*150;
    end
    %Q(1,i+1) = Q(2,i+1);
end
% figure
plot(t,Q)
    % alphaM = (alpha(i+1,j)+alpha(i,j+1))/2
    % betaM =(beta(i+1,j)+beta(i,j+1))/2
    %
    % Q(i+1,j+1)= (Q(i+1,j)+ dt/dx*alphaM*Q(i,j+1)-betaM*dt)/(a+alphaM*dx/dt)

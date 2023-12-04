x = 0:1:50;
B = 50
t = 0:150:54000;
g = 9.8
S0 = .0002
n = .01
V = 4
C =1

beta = zeros(length(time),length(x));
alpha = zeros(length(time),length(x));
Q = zeros(length(time),length(x));
%random once again bc we cant find one
Q0 = 2
%random values for amplitude and period lol
amp = 3
T = 10
Q(:,1) = Q0*(1+ amp*sin((2*pi*t)/T))'


for i = 2:length(time)
    for j = 2:length(x)-1
   
    A =  B * x(j);
    %R is guestimate lol for perimeter
    R = A/58;
    Sf = sqrt((V*n)/R^(2/3));
    alpha(i,j) = 2*Q(i,j)/A +((((g*A)/B)-Q(i,j).^2/A.^2))/((Q(i,j)/A)* (5/3-(4.*R)/3.*B));
    beta(i,j) = g*A*(Sf-S0);

    % Q(i,j+1) = Q(i,j) - C*(Q(i,j)-Q(i,j-1));
    % Q(i,j+1) = Q(i,j)-alpha(i,j)*(150/1)*(Q(i,j)-Q(i-1,j))-beta(i,j)*150;
    end
    Q(1,j+1) = Q(2,j+1);
end

    % alphaM = (alpha(i+1,j)+alpha(i,j+1))/2
    % betaM =(beta(i+1,j)+beta(i,j+1))/2
    %
    % Q(i+1,j+1)= (Q(i+1,j)+ dt/dx*alphaM*Q(i,j+1)-betaM*dt)/(a+alphaM*dx/dt)

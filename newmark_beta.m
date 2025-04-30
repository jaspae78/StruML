% Parameters
M = [2 0; 0 1];           % Mass matrix
K = [400 -200; -200 200]; % Stiffness matrix
dt = 0.01; T = 2;         % Time step and total time
N = T/dt;

% Newmark parameters (unconditionally stable)
beta = 0.25;
gamma = 0.5;

% Initial conditions
u = zeros(2,N); v = zeros(2,N); a = zeros(2,N);
f = zeros(2,N); f(2,:) = 10*sin(2*pi*(0:dt:T-dt)); % Force on 2nd DOF

% Initial acceleration
a(:,1) = M \ (f(:,1) - K*u(:,1));

% Constants
a0 = 1/(beta*dt^2); a1 = gamma/(beta*dt); a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1; a4 = gamma/beta - 1; a5 = dt*(gamma/(2*beta) - 1);

% Effective stiffness
Keff = K + a0*M;

for n = 1:N-1
    f_eff = f(:,n+1) + M*(a0*u(:,n) + a2*v(:,n) + a3*a(:,n));
    u(:,n+1) = Keff \ f_eff;
    a(:,n+1) = a0*(u(:,n+1) - u(:,n)) - a2*v(:,n) - a3*a(:,n);
    v(:,n+1) = v(:,n) + dt*((1 - gamma)*a(:,n) + gamma*a(:,n+1));
end

% Plot displacement
t = 0:dt:T-dt;
plot(t, u(1,:), 'b', t, u(2,:), 'r--')
legend('DOF 1', 'DOF 2'); xlabel('Time [s]'); ylabel('Displacement');
title('Newmark-beta Method')

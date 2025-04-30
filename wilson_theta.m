% Wilson-theta Method for 2 DOF System

clc; clear;

% System parameters
M = [2 0; 0 1];                   % Mass matrix
K = [400 -200; -200 200];        % Stiffness matrix

% Time parameters
dt = 0.01; T = 2;                 % Time step and total time
N = floor(T/dt);                 % Number of steps
theta = 1.4;                     % θ parameter (≥ 1.37 for stability)
dt_theta = theta * dt;

% Initialize arrays
u = zeros(2,N);                  % Displacement
v = zeros(2,N);                  % Velocity
a = zeros(2,N);                  % Acceleration
f = zeros(2,N);                  % External force

% Apply sinusoidal force to 2nd DOF
t = 0:dt:T-dt;
f(2,:) = 10 * sin(2*pi*t);

% Initial acceleration
a(:,1) = M \ (f(:,1) - K*u(:,1));

% Precompute constants
a0 = 6 / (dt_theta^2);
a1 = 3 / dt_theta;
a2 = 2 * a0;

% Effective stiffness matrix
Keff = K + a0 * M;

for n = 1:N-1
    % Interpolate external force at t + θΔt
    if n+1 <= N
        f_theta = f(:,n) + theta * (f(:,n+1) - f(:,n));
    else
        f_theta = f(:,n);  % Avoid indexing past end
    end

    % Effective force
    R = f_theta + M * (a0 * u(:,n) + a1 * v(:,n) + 2 * a(:,n));

    % Solve for displacement at t + θΔt
    u_theta = Keff \ R;

    % Compute acceleration and velocity at t + Δt
    a_next = a0 * (u_theta - u(:,n)) - a1 * v(:,n) - 2 * a(:,n);
    v_next = v(:,n) + (dt / 2) * (a(:,n) + a_next);
    u_next = u(:,n) + dt * v(:,n) + (dt^2 / 6) * (a(:,n) + a_next + 4 * a(:,n));  % Midpoint rule

    % Store results
    u(:,n+1) = u_next;
    v(:,n+1) = v_next;
    a(:,n+1) = a_next;
end

% Plot displacement
figure;
plot(t, u(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, u(2,:), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Displacement (m)');
legend('DOF 1','DOF 2');
title('Wilson-\theta Method (θ = 1.4)');
grid on;

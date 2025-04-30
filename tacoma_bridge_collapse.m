% Tacoma Narrows Bridge - Torsional Flutter Simulation
clc; clear; close all;

%% Parameters
m = 1000;         % mass (kg)
I = 5000;         % moment of inertia (kg*m^2)
k_y = 20000;      % vertical stiffness (N/m)
k_theta = 10000;  % torsional stiffness (N*m/rad)
c_y = 200;        % vertical damping (N·s/m)
c_theta = 150;    % torsional damping (N·m·s/rad)

U = 40;           % wind speed (m/s)
rho = 1.225;      % air density (kg/m^3)
B = 10;           % bridge width (m)
CL_alpha = 5;     % lift slope (approximate, per rad)

% Time setup
dt = 0.01; T = 60; t = 0:dt:T;

% Preallocate
y = zeros(size(t));      % vertical motion
theta = zeros(size(t));  % torsional motion
v_y = 0; v_theta = 0;

% Initial conditions
y(1) = 0.01;             % small initial displacement
theta(1) = 0.01;         % small initial twist

%% Time Integration using Euler method
for i = 1:length(t)-1
    % Aerodynamic lift (force) and moment (simplified)
    L = 0.5 * rho * U^2 * B * CL_alpha * theta(i);          % Lift ~ θ
    M = 0.5 * rho * U^2 * B^2 * CL_alpha * theta(i);        % Moment ~ θ

    % Equations of motion
    a_y = ( -c_y * v_y - k_y * y(i) + L ) / m;
    a_theta = ( -c_theta * v_theta - k_theta * theta(i) + M ) / I;

    % Update velocities
    v_y = v_y + a_y * dt;
    v_theta = v_theta + a_theta * dt;

    % Update displacements
    y(i+1) = y(i) + v_y * dt;
    theta(i+1) = theta(i) + v_theta * dt;
end

%% Plot Results
figure;
subplot(2,1,1)
plot(t, y, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Vertical Displacement (m)');
title('Vertical Motion of Bridge Deck');
grid on;

subplot(2,1,2)
plot(t, theta, 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Torsional Displacement (rad)');
title('Torsional Flutter (Twisting Motion)');
grid on;

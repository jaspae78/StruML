
% Jeffcott Rotor Simulation with Orbit Plot and Campbell Diagram
% Jasmine Pae - Aerospace Application of Vibration

clc; clear; close all;

%% Parameters
m = 1;              % kg, rotor mass
k = 1000;           % N/m, shaft stiffness
c = 5;              % Ns/m, damping
e = 0.01;           % m, eccentricity (imbalance)
omega_range = linspace(5, 100, 300);  % rad/s

wn = sqrt(k/m);     % natural frequency
X_response = zeros(size(omega_range));

%% Steady-State Amplitude vs Speed (Campbell Data)
for i = 1:length(omega_range)
    omega = omega_range(i);
    num = m * e * omega^2;
    denom = sqrt((k - m * omega^2)^2 + (c * omega)^2);
    X_response(i) = num / denom;
end

%% Orbit Plot at Resonance
omega_res = wn;     % simulate at critical speed
t = linspace(0, 2*pi/omega_res*5, 1000);  % simulate 5 cycles
F = m * e * omega_res^2;
A = F / sqrt((k - m * omega_res^2)^2 + (c * omega_res)^2);
phi = atan2(c * omega_res, (k - m * omega_res^2));

x = A * cos(omega_res * t);
y = A * sin(omega_res * t + phi);  % 90° out of phase creates orbit

%% Plot Orbit
figure;
plot(x, y, 'b', 'LineWidth', 2);
axis equal;
xlabel('X displacement (m)');
ylabel('Y displacement (m)');
title('Orbit Plot at Critical Speed');
grid on;

%% Plot Campbell Diagram (Amplitude vs Speed)
figure;
plot(omega_range, X_response, 'r', 'LineWidth', 2);
xlabel('Rotor Speed ω (rad/s)');
ylabel('Amplitude X (m)');
title('Campbell Diagram - Steady-State Response');
grid on;

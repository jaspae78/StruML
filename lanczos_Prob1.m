
% --------- MAIN SCRIPT ---------
clc;
clear;

% This is a sample MATLAB code to solve an eigenvalue problem using the Lanczos method

% Define symmetric matrix A (4x4 stiffness-like)
A = [6 -2 0 0;
    -2 6 -2 0;
     0 -2 6 -2;
     0 0 -2 6];

% Initial guess vector (random and normalized)
q1 = rand(4, 1);
q1 = q1 / norm(q1);

% Number of Lanczos iterations (number of modes to extract)
m = 3;

% Call Lanczos Method
[T, Q] = lanczos_method(A, q1, m);

% Display results
fprintf('Tridiagonal matrix T:\n');
disp(T);

% Compute approximate eigenvalues of A from T
[evecs_T, evals_T] = eig(T);

fprintf('Approximate eigenvalues of A:\n');
disp(diag(evals_T));


% --------- FUNCTION DEFINITIONS BELOW ---------

% Lanczos method for symmetric matrix A
function [T, Q] = lanczos_method(A, q1, m)
    n = length(q1);
    Q = zeros(n, m);
    alpha = zeros(m, 1);
    beta = zeros(m-1, 1);

    q1 = q1 / norm(q1);  % Normalize the initial vector
    Q(:,1) = q1;

    for j = 1:m
        if j == 1
            w = A * Q(:,j);
        else
            w = A * Q(:,j) - beta(j-1) * Q(:,j-1);
        end

        alpha(j) = Q(:,j)' * w;
        w = w - alpha(j) * Q(:,j);

        if j < m
            beta(j) = norm(w);
            if beta(j) ~= 0
                Q(:,j+1) = w / beta(j);
            end
        end
    end

    % Construct tridiagonal matrix T
    T = diag(alpha) + diag(beta,1) + diag(beta,-1);
end

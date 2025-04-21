% --------- MAIN SCRIPT ---------
clc;
clear;

% Define symmetric matrix A (example)
A = [4 1 2; 1 3 0; 2 0 1];  % 3x3 symmetric matrix for test
tol = 1e-6;                 % Tolerance for convergence
max_iter = 100;             % Maximum number of iterations

% Call Jacobi Eigenvalue Solver
[D, V] = jacobi_eigen(A, tol, max_iter);

% Sort eigenvalues for clean output
D = sort(D);

% Display results
fprintf('Jacobi Method Results:\n');
fprintf('----------------------\n');
fprintf('Eigenvalues:\n');
disp(D);

fprintf('Eigenvectors (columns):\n');
disp(V);


% --------- FUNCTION DEFINITIONS BELOW ---------

% Jacobi Eigenvalue Solver
% Inputs:
%   A        - symmetric matrix
%   tol      - convergence tolerance
%   max_iter - maximum number of iterations
% Outputs:
%   D        - vector of eigenvalues
%   V        - matrix of eigenvectors (columns)
function [D, V] = jacobi_eigen(A, tol, max_iter)
    n = size(A, 1);
    V = eye(n);  % Initialize eigenvector matrix

    for iter = 1:max_iter
        [max_val, p, q] = max_offdiag(A);  % Find largest off-diagonal

        if abs(max_val) < tol
            break;
        end

        % Compute rotation angle
        if A(p,p) == A(q,q)
            theta = pi/4;
        else
            theta = 0.5 * atan(2*A(p,q)/(A(q,q) - A(p,p)));
        end

        % Construct rotation matrix
        c = cos(theta);
        s = sin(theta);
        J = eye(n);
        J([p, q], [p, q]) = [c, -s; s, c];

        % Apply rotation
        A = J' * A * J;
        V = V * J;  % Accumulate eigenvectors
    end

    D = diag(A);  % Extract eigenvalues from diagonal
end

% Helper function to find largest off-diagonal element
function [max_val, p, q] = max_offdiag(A)
    n = size(A,1);
    max_val = 0;
    p = 1;
    q = 2;

    for i = 1:n
        for j = i+1:n
            if abs(A(i,j)) > abs(max_val)
                max_val = A(i,j);
                p = i;
                q = j;
            end
        end
    end
end

%% Function for Legendre-Gauss abscissae and weights
function [x,w]=lgwt(N,a,b)
% lgwt computes the Legendre-Gauss nodes (abscissae) and weights for integration
% over the interval [a, b].

% Inputs:
% - N: Number of Legendre-Gauss points
% - a, b: Interval for the quadrature (integration over [a, b])

% Outputs:
% - x: The N Legendre-Gauss abscissae (nodes) over the interval [a, b]
% - w: The corresponding weights for Legendre-Gauss quadrature

N=N-1; % Adjust for zero-based indexing
N1=N+1; N2=N+2; % Constants for use in the algorithm
xu=linspace(-1,1,N1)'; % Initial equally spaced points between -1 and 1

% Step 1: Initial guess for the Legendre-Gauss nodes using a cosine approximation
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Step 2: Initialize the Legendre-Gauss Vandermonde Matrix (LGVM)
L=zeros(N1,N2); % Legendre polynomial values
Lp=zeros(N1,N2); % Derivatives of Legendre polynomials

% Step 3: Compute the zeros of the N+1 Legendre polynomial using
% the recursion relation and Newton-Raphson method
y0=2; % Initial value to ensure convergence

% Iterate until the new points are within a small tolerance (epsilon)
while max(abs(y-y0))>eps
    L(:,1)=1; % Initialize the first Legendre polynomial (P0 = 1)
    Lp(:,1)=0; % Derivative of P0
    
    L(:,2)=y; % Initialize the second Legendre polynomial (P1 = x)
    Lp(:,2)=1; % Derivative of P1
    
    % Recursively compute higher-order Legendre polynomials
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k; 
        % Standard recursion relation for Legendre polynomials
    end
    
    % Compute the derivative of the N+1th Legendre polynomial
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2); 
    % This computes the derivative based on the recursion relation
    
    % Newton-Raphson update for the roots of the Legendre polynomial
    y0=y; % Store the old values of y
    y=y0-L(:,N2)./Lp; % Update y using Newton's method
end

%% Step 4: Linear mapping from [-1, 1] to [a, b]
x=(a*(1-y)+b*(1+y))/2; % Map the computed nodes from the interval [-1,1] to [a,b]

%% Step 5: Compute the weights for the quadrature
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2; % Compute the weights for Legendre-Gauss quadrature

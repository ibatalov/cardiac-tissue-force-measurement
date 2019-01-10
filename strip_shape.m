%% returns coordinates of an elastic strip (or rod) of a lenght L, bent at both ends by the angle alpha in radians
% N - number of (X,Y) pairs to return
function [X,Y] = strip_shape(L, alpha, N)

% k = sqrt(F/EI), F - force, E - Young's modulus, I - inertia moment
k = 1/L*integral(@(x) (2.*(cos(x) - cos(alpha))).^(-0.5), -alpha, alpha);

X = zeros(N, 1);
Y = zeros(N, 1);

for n = 1 : N
    theta = -alpha*(2*(n - 1)/(N - 1) - 1);
    if(theta < alpha)
        X(n) = 1./k.*integral(@(x) cos(x)./(2.*(cos(x) - cos(alpha))).^0.5, theta, alpha);
    else
        X(n) = 0;
    end
    Y(n) = 1./k.*(2.*(cos(theta) - cos(alpha))).^0.5;
end
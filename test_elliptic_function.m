% y = 1/k*(2*(cos(Theta) - cos(Alpha))^0.5
% x = 1/k*(integral cos(theta')/(2(cos(theta') - cos(alpha))^0.5 d(Theta)from theta to alpha)
% k = (F/EI)^0.5
N = 200;

figure;
L_data = zeros(10,1);

hold on
n_curves = 10;
max_n = n_curves + 1;

all_x = zeros(n_curves, N);
all_y = zeros(n_curves, N);

for m = 1 : 10
    
    L = m/10;
    alpha = pi*2/4;
    k = 1/L*integral(@(x) (2.*(cos(x) - cos(alpha))).^(-0.5), -alpha, alpha);
    
    theta_space = linspace(-alpha, alpha, N);
    X = zeros(N, 1);
    Y = zeros(N, 1);
    
    for n = 1 : N
        theta = theta_space(n);
        if(theta < alpha)
            X(n) = 1./k.*integral(@(x) cos(x)./(2.*(cos(x) - cos(alpha))).^0.5, theta, alpha);
        else
            X(n) = 0;
        end
        Y(n) = 1./k.*(2.*(cos(theta) - cos(alpha))).^0.5;
    end
    
    all_x(m, :) = X;
    all_y(m, :) = Y;
    
    X = X - X(round(N/2));    
    plot(X, Y);
end
hold off
axis([-0.5 0.5 0 0.6]);  

L_coords = linspace(0.1, 1, 10);
ratios = (L_coords')./(all_x(:, N) - all_x(:, 1));


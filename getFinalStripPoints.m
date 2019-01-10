%% L - curved strip length
%  l - length of the straight segments of the strip
%  alpha = (0,pi) - angle of bending (measured at each end, total angle is thaice higher)
%  (r0, c0) - row and column of the beginning of the curved part of the strip
%  phi = [0,pi) - strip rotation angle
function [rows, cols] = getFinalStripPoints(N_all, L, l, alpha, r0, c0, phi)

N_end = round(N_all*l/(L + l)/2);
N_bent = N_all - N_end*2;

[x_bent, y_bent] = strip_shape(L, alpha, N_bent);

a1 = [-cos(alpha); -sin(alpha)];
a2 = [cos(alpha); -sin(alpha)];
left_end = zeros(N_end, 2);
right_end = zeros(N_end, 2);

for n = 1 : N_end
    left_end(n, :) = [x_bent(1); y_bent(1)] + a1*n*l/N_end;
    right_end(n, :) = [x_bent(N_bent); y_bent(N_bent)] + a2*n*l/N_end;
end

x_all = [left_end(:, 1); x_bent; right_end(:,1)];
y_all = [left_end(:, 2); y_bent; right_end(:,2)];

rotation_matrix = [cos(phi), -sin(phi); sin(phi), cos(phi)];

r = zeros(2, length(x_all));
r(1,:) = x_all;
r(2,:) = y_all;

r = rotation_matrix*r;
r(1, :) = r(1,:) + c0;
r(2, :) = r(2,:) - r0;

rows = - r(2, :);
cols =   r(1, :);

% remove elements <= 0
cols = cols(round(rows) > 0);
rows = rows(round(rows) > 0);

rows = rows(round(cols) > 0);
cols = cols(round(cols) > 0);

% figure('Position', [100, 100, 600, 600]);
% plot(cols, rows, 'Marker', '.', 'LineStyle', 'none');
% %axis([110, 160, -250, -200]);
% axis square;

% rows1 = round(rows);
% cols1 = round(cols);
% image = zeros(max(rows1), max(cols1));
% image(sub2ind(size(image), rows1, cols1)) = 1;
% imshow(image);

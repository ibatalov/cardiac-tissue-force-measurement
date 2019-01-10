function [slice_coords] = slice_3d_model(coords, slice_distance)
% Cnverts 3D triangle coords into a bunch of polygons for each z-slice
%   Detailed explanation goes here
%	coords - input triangle coordinates. Format:
%	x11, y11, z11;	- first point of triangle 1
%	x12, y12, z12;	- second point of triangle 1
%	x13, y13, z13;	- third point of triangle 1
%	x21, y22, z23;	- first point of triangle 2
%	etc.
%	slice_coords - cell(slice_count, 1)
%	Each cell is a cell array containing coordinates of all polygons in the slice. 
%	Each polygon is an array of the following structure:
%	x1, y1;
%	x2, y2;
%	etc.
%	Here is another representation of the data structure:
%	slice_coords{slice_number}{polygon_number}(point_number, coord_number(=1 for x, =2 for y))

%% This is a test array. Comment it out if want to test the function.
square_coords = [1,1,1; -1,1,1; 1,-1,1; -1,-1,1; 1,1,-1; -1,1,-1; 1,-1,-1; -1,-1,-1];
triangle_points = [1,2,4; 1,3,4; 5,6,8; 5,7,8; 1,5,7; 1,7,3; 3,7,8; 3,8,4; 4,8,6; 4,6,2; 1,5,6; 1,6,2];
triangles_coords = zeros(numel(triangle_points), 3);
for row = 1 : size(triangle_points,1)
	for col = 1 : size(triangle_points, 2)
		triangles_coords(col + (row-1)*3, :) = square_coords(triangle_points(row, col),:);
	end
end
angle = 0; % in radians
rotM = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
triangles_coords = (rotM * (triangles_coords.')).';
coords = triangles_coords;
slice_distance = 0.1;

%%
min_coords = min(coords, [], 1);
max_coords = max(coords, [], 1);
slice_count = ceil((max_coords(3) - min_coords(3))/slice_distance) - 1;
slice_z_positions = linspace(min_coords(3) + slice_distance/2, max_coords(3) - slice_distance/2, slice_count);
slice_z_positions = slice_z_positions(:);

slice_coords = cell(slice_count, 1);

for slice_num = 1 : slice_count
	curr_z = slice_z_positions(slice_num);
	slice_points = [];
	% find triangles intersecting the z-plane
	for triangle_num = 1 : (size(coords, 1)/3 - 1);
		triangle_coords = coords(3*(triangle_num - 1) + 1 : 3*triangle_num, :);
		z_deltas = triangle_coords(:, 3) - curr_z;
		z_deltas = z_deltas ./ abs(z_deltas);
		if min(z_deltas) < 0 && max(z_deltas) > 0
			% this triangle does cross the current z-plane
			% now finding the intersection of the triangle and the plane
			if sum(z_deltas) > 0
				% 2 points are above, 1 - below
				start_point_number = find(z_deltas < 0);
			else
				% 1 point is above, 2 - below
				start_point_number = find(z_deltas > 0);
			end
			finish_point_numbers = [1,2,3];
			finish_point_numbers(start_point_number) = []; % remove the starting point
			
			dirs = triangle_coords(finish_point_numbers, :) - triangle_coords(start_point_number, :);
			coeffs = (curr_z - triangle_coords(start_point_number, 3)) ./ dirs(:, 3);
			coeffs = coeffs(:);
			intersection_points = repmat(triangle_coords(start_point_number, :), 2, 1) + repmat(coeffs, 1, 3) .* dirs;
			slice_points = [slice_points; intersection_points];
		else
			if min(z_deltas) == 0 || max(z_deltas) == 0 && abs(sum(z_deltas)) < 0.5
				% at least 2 vertices of this triangle lay in the plane
				touching_point_numbers = find(z_deltas == 0);
				intersection_points = triangle_coords(touching_point_numbers, :);
				slice_points = [slice_points; intersection_points];
			end	
		end
	end
	slice_coords{slice_num} = slice_points;
end

%% visualize slices
slice_num = 10;
figure;
plot(slice_coords{slice_num}(:,1), slice_coords{slice_num}(:,2), 'o');


end


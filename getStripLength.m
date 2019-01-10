function strip_length = getStripLength(strip_image, attachment_points)

% keeping only the biggest object in the image
cleared_strip = strip_image;
CC = bwconncomp(strip_image);
if(CC.NumObjects > 1)
    max_size = numel(CC.PixelIdxList{1}(:));
    max_n = 1;
    for n = 2 : CC.NumObjects
        size_n = numel(CC.PixelIdxList{n}(:));
        if(size_n > max_size)
            max_size = size_n;
            cleared_strip(CC.PixelIdxList{max_n}(:)) = 0;
            max_n = n;
        else
            cleared_strip(CC.PixelIdxList{n}(:)) = 0;
        end
    end
end

strip_skel = bwmorph(cleared_strip, 'skel', Inf);
[branch_row, branch_col] = find(bwmorph(strip_skel, 'branchpoints'));

% in case the skeleton branches, remove all the side branches, keep
% only the longest branch
if(~isempty(branch_row))
    for point_num = 1 : length(branch_row)
        temp_image = strip_skel;
        min_row = max(branch_row(point_num) - 1, 1);
        max_row = min(branch_row(point_num) + 1, size(temp_image, 1));
        min_col = max(branch_col(point_num) - 1, 1);
        max_col = min(branch_col(point_num) + 1, size(temp_image, 2));
        
        temp_image(min_row : max_row, min_col : max_col) = 0;
        
        % go through all branches for the current branching point, keep
        % only 2 longest
        CC = bwconncomp(temp_image);
        if(CC.NumObjects > 2)
            size_info = zeros(CC.NumObjects, 2);
            for n = 1 : CC.NumObjects
                size_info(n, 2) = numel(CC.PixelIdxList{n}(:));
                size_info(n, 1) = n;
            end
            
            size_info = sortrows(size_info, -2); % sort rows of the matrix based on the 2nd column in the descending order
            for obj_n = 3 : CC.NumObjects
                strip_skel(CC.PixelIdxList{size_info(obj_n,1)}) = 0;
            end
            if(size_info(2,2)/size_info(1,2) < 0.05)
                strip_skel(CC.PixelIdxList{size_info(2,1)}) = 0;
            end
        end
    end
end

end_points = find(bwmorph(strip_skel, 'endpoints'));
while(length(end_points) > 2)
    strip_skel(end_points) = 0;
    end_points = find(bwmorph(strip_skel, 'endpoints'));
end

strip_skel(round(attachment_points{2}(1)), round(attachment_points{2}(2))) = 1;
temp_skel = strip_skel;
sorted_skel = [round(attachment_points{1}(1)), round(attachment_points{1}(2))];
while (sorted_skel(end, 1) ~= round(attachment_points{2}(1)) || sorted_skel(end, 2) ~= round(attachment_points{2}(2))) && nnz(temp_skel) > 0
    r = sorted_skel(end, 1);
    c = sorted_skel(end, 2);
    min_r = max(r - 1, 1);
    max_r = min(r + 1, size(temp_skel, 1));
    min_c = max(c - 1, 1);
    max_c = min(c + 1, size(temp_skel, 2));
    
    if(nnz(temp_skel(min_r : max_r, min_c : max_c)) > 0)
        [next_r, next_c] = find(temp_skel(min_r : max_r, min_c : max_c));
        next_r = next_r(1) + r - 2;
        next_c = next_c(1) + c - 2;
    else
        [row_set, col_set] = find(temp_skel);
        sq_dist_set = (row_set - r).^2 + (col_set - c).^2;
        m_for_sorting = zeros(length(sq_dist_set),3);
        m_for_sorting(:,1) = row_set(:);
        m_for_sorting(:,2) = col_set(:);
        m_for_sorting(:,3) = sq_dist_set;
        m_for_sorting = sortrows(m_for_sorting, 3); % ascending sorting by the 3rd column
        next_r = m_for_sorting(1,1);
        next_c = m_for_sorting(1,2);
    end
    temp_skel(next_r, next_c) = 0;
    sorted_skel = [sorted_skel; next_r, next_c];
end
clearvars temp_skel;

strip_rows = sorted_skel(:,1);
strip_cols = sorted_skel(:,2);
shift_row = [strip_rows(1); strip_rows(1:end-1)];
shift_col = [strip_cols(1); strip_cols(1:end-1)];
l_strip = ((strip_rows - shift_row).^2 + (strip_cols - shift_col).^2).^0.5;

l_cumul = l_strip;
for i = 1 : length(l_strip) - 1
    l_cumul(i+1:end) = l_cumul(i+1:end) + l_strip(i);
end

row_poly = polyfit(l_cumul, strip_rows, 7);
col_poly = polyfit(l_cumul, strip_cols, 7);
l_cumul = l_cumul(:).';
l_vector = [l_cumul.^7; l_cumul.^6; l_cumul.^5; l_cumul.^4; l_cumul.^3; l_cumul.^2; l_cumul; ones(1,length(l_cumul))];

% figure;
% plot(l_cumul, row_poly*l_vector)
% hold on
% plot(l_cumul, strip_rows)
% plot(l_cumul, col_poly*l_vector)
% plot(l_cumul, strip_cols)
% hold off

strip_length = 0;
fit_strip = zeros(size(strip_image));
for t = 0 : l_cumul(end)
    t_column = [t.^7; t.^6; t.^5; t.^4; t.^3; t.^2; t; 1];
    curr_row = row_poly*t_column;
    curr_col = col_poly*t_column;
    if(t > 0)
        strip_length = strip_length + ((curr_row - prev_row).^2 + (curr_col - prev_col).^2).^0.5;   
    end
    prev_row = curr_row;
    prev_col = curr_col;
    if(curr_row > 1 && curr_row < size(fit_strip, 1) && curr_col > 1 && curr_col < size(fit_strip, 2))
        fit_strip(round(curr_row), round(curr_col)) = 1;
    end
end

strip_overlay = zeros(size(strip_image, 1), size(strip_image, 2), 3);
strip_overlay(:,:,1) = strip_image;
strip_overlay(:,:,2) = fit_strip;

figure;
imshow(strip_overlay);

end
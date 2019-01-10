E_strip = 3.1*10^6; % strip elastic modulus, Pa
w_strip = 1.5; % strip width, mm
t_strip = 0.125; % strip thickness, mm

samplingFactor = 4; % integer > 0. Best if it's 2^n, n = N.
cutoff = 0.135; % fraction of pixels that should be above the threshold

mm_per_px = 0.00849*samplingFactor;
L_bent_min_mm = 12.6;
L_bent_max_mm = 16.6;
L_end_min_mm = 3;
L_end_max_mm = 5;

L_strip = 0;
%L_strip = (L_bent_min_mm + L_bent_max_mm)/2; % strip length between tissue attachment points
I_strip = (t_strip/1000)^3*(w_strip/1000)/12; % second moment of area of the strip's cross-section, m^4

%% opeining file
if(exist(path, 'dir'))
    [file, tempPath] = uigetfile({'*.avi; *.mov'},'Select file to open', path);
else
    [file, tempPath] = uigetfile({'*.avi; *.mov'});
end
if exist([tempPath file], 'file')
    path = tempPath;
    vidObj = VideoReader([path file]);
    frame_rate = vidObj.FrameRate;
    video = zeros(vidObj.Height/samplingFactor, vidObj.Width/samplingFactor, floor(vidObj.Duration*vidObj.FrameRate));
    size(video)
    frame = 1;
    
    while hasFrame(vidObj)
        frameImage = readFrame(vidObj);
        frameImage = rgb2gray(frameImage);
        video(:,:,frame) = frameImage(1:samplingFactor:end, 1:samplingFactor:end);
        frame = frame + 1;
    end
    
    video = video/max(video(:));
    %implay(video);
else
    disp('file not found');
end

%% thresholding the video based on the cutoff fraction
% randomly chose frame 52 as an example of a frame
image = video (:,:,52);
% figure;
% imshow(image);
bin_width = 0.01;

user_cutoff = cutoff;
while(~isempty(user_cutoff))
    cutoff = user_cutoff;
    temp_cutoff = 0;
    threshold = 1;
    above = 0;
    total = numel(image);
    
    while(temp_cutoff < cutoff)
        threshold = threshold - bin_width;
        above = nnz(image > threshold);
        temp_cutoff = above / total;
    end
    
    binary_video = video > threshold;
    implay(binary_video);
    user_cutoff = input(['current cutoff fraction: ' num2str(cutoff) '. Enter new cutoff or press enter if satisfied.']);
end

%% find optimal erode length
frame_num = 2;
frame = binary_video(:,:,frame_num);
frame = bwmorph(frame, 'dilate', ceil(size(video, 1)/200));
frame = bwmorph(frame, 'erode', ceil(size(video, 1)/200));

CC = bwconncomp(frame);
max_size = numel(CC.PixelIdxList{1}(:));
max_n = 1;
for n = 2 : CC.NumObjects
    size_n = numel(CC.PixelIdxList{n}(:));
    if(size_n > max_size)
        max_size = size_n;
        frame(CC.PixelIdxList{max_n}(:)) = 0;
        max_n = n;
    else
        frame(CC.PixelIdxList{n}(:)) = 0;
    end
end

CC = bwconncomp(~frame);
for n = 1 : CC.NumObjects
    size_n = numel(CC.PixelIdxList{n}(:));
    if(size_n < size(video, 1)*size(video, 2)/50)
        frame(CC.PixelIdxList{n}(:)) = 1;
    end
end

good_erode_sizes = [];
total_convex_size = nnz(bwconvhull(frame));
for erode_size = round(size(video,2)/70) : round(size(video,2)/10);
    disk_dilate = strel('disk', round(erode_size*1.5), 0);

    %tissue = imdilate(bwmorph(frame, 'erode', erode_size), disk);
    tissue = imdilate(custom_erode(frame, erode_size), disk_dilate);
    convex_tissue = bwconvhull(tissue);
    
    ratio1 = nnz(tissue)/nnz(convex_tissue);
    if(~isempty(good_erode_sizes))
        old_ratio2 = ratio2;
    end
    ratio2 = nnz(convex_tissue)/total_convex_size;
    if(ratio1 > 0.65 && ratio2 < 0.3 && ratio2 > 0.10)
        if(~isempty(good_erode_sizes))
            if(ratio2/old_ratio2 < 0.86)
                break;
            end
        end
        good_erode_sizes = [good_erode_sizes; erode_size];
        disp([num2str(erode_size) '; r1: ' num2str(ratio1) '; r2: ' num2str(ratio2) '; good']);
    else
        disp([num2str(erode_size) '; r1: ' num2str(ratio1) '; r2: ' num2str(ratio2) '; bad']);
        if(~isempty(good_erode_sizes))
            break;
        end
    end
end
erode_size = sum(good_erode_sizes)/numel(good_erode_sizes);
erode_size = round(erode_size);

erode_input = erode_size;
while(~isempty(erode_input) && isnumeric(erode_input))
    erode_size = erode_input;
    disk_dilate = strel('disk', round(erode_size*1.5), 0);
    [eroded_image, tissue_broken] = custom_erode(frame, erode_size);
    tissue = imdilate(eroded_image, disk_dilate);
    figure;
    imshow(tissue);
    erode_input = input(['erode length: ' num2str(erode_size) '. Press enter if satisfied or the new erode length if not: ']);
end

%use_convex_tissue = input('Do you want to use convex tissue? (1 = yes, 0 = no): ');
use_convex_tissue = tissue_broken;

%% main loop processing the video frame by frame
measured_data = zeros(size(video, 3), 8); %frame number, row1, col1, row2, col2, distance, alpha, force
tissue_video = zeros(size(binary_video));
%estimated_points_video = zeros(size(binary_video));
combined_video = zeros(size(video, 1), size(video, 2), 3, size(video, 3));
%fitted_strip_parameters = cell(size(binary_video,3));
%messed_video = zeros(size(binary_video));

%erode_size = round(size(video,2)/50);
user_erode_size = erode_size;
disk_dilate = strel('disk', round(erode_size*1.5), 0);
%N = 50;
t0 = clock;

for frame_num = 1 : size(video, 3)
    disp(frame_num);
    frame = binary_video(:,:,frame_num);
    %% connect objects that are close to each other
    frame = bwmorph(frame, 'dilate', ceil(size(video, 1)/200));
    frame = bwmorph(frame, 'erode', ceil(size(video, 1)/200));
    
    %% remove small shit. Keep only the biggest object - the tissue with the strip
    CC = bwconncomp(frame);
    if(CC.NumObjects > 0)
        max_size = numel(CC.PixelIdxList{1}(:));
        max_n = 1;
        for n = 2 : CC.NumObjects
            size_n = numel(CC.PixelIdxList{n}(:));
            if(size_n > max_size)
                max_size = size_n;
                frame(CC.PixelIdxList{max_n}(:)) = 0;
                max_n = n;
            else
                frame(CC.PixelIdxList{n}(:)) = 0;
            end
        end
        
        %% remove holes in the strip and the tissue
        max_hole_area = (2/mm_per_px)^2;
        CC = bwconncomp(~frame);
        for n = 1 : CC.NumObjects
            size_n = numel(CC.PixelIdxList{n}(:));
            if(size_n < max_hole_area)
                frame(CC.PixelIdxList{n}(:)) = 1;
            end
        end
        binary_video(:,:,frame_num) = frame;
        combined_video(:,:,1,frame_num) = frame;
        %        tissue = imdilate(bwmorph(frame, 'erode', erode_size), disk);
        [eroded_tissue, tissue_broke] = custom_erode(frame, erode_size);
        tissue = imdilate(eroded_tissue, disk_dilate);
        if(tissue_broke && use_convex_tissue == 0)
            use_convex_tissue = 1;
            disp('Tissue broke! Using convex tissue from now on.');
        end
        
        if(use_convex_tissue == 1)
            convex_tissue = bwconvhull(tissue);
        else
            convex_tissue = tissue;
        end
        
        strip = frame & ~convex_tissue;
        tissue_video(:,:,frame_num) = tissue;
        
        %% find strip within tissue
        anchors_image = strip & imdilate(convex_tissue, disk_dilate);
        center_list = [];
        CC = bwconncomp(anchors_image);
        for object = 1 : CC.NumObjects
            [rows, cols] = ind2sub(size(anchors_image),CC.PixelIdxList{object});
            num_pixels = numel(rows);
            center_row = sum(rows)/num_pixels;
            center_col = sum(cols)/num_pixels;
            center_list = [center_list; num_pixels, center_row, center_col];
        end
        center_list = sortrows(center_list, -1);
        center_list = center_list(1:4, :);
        
        point_nums = cell(2, 1);
        min_dist_1 = numel(anchors_image);
        min_dist_2 = min_dist_1;
        
        for i = 1 : 3
            for j = i + 1 : 4
                sq_dist = (center_list(i,2) - center_list(j,2))^2 + (center_list(i,3) - center_list(j,3))^2;
                if(sq_dist < max(min_dist_1, min_dist_2))
                    if(min_dist_1 < min_dist_2)
                        min_dist_2 = sq_dist;
                        point_nums{2} = [i, j];
                    else
                        min_dist_1 = sq_dist;
                        point_nums{1} = [i, j];
                    end
                end
            end
        end
        
        %% find tissue center and inertia moments
        CC = bwconncomp(convex_tissue);
        if(CC.NumObjects ~= 1)
            disp('ERROR! More than one object is detected for the tissue');
        end
        [raw_tensor, center] = getInertiaTensor(convex_tissue);
        [eigenVectors, diagMatrix] = eig(raw_tensor);
        %area = nnz(convex_tissue);
        
        if diagMatrix(1,1) > diagMatrix(2,2)
            maxMoment = diagMatrix(1,1);
            minMoment = diagMatrix(2,2);
            mainDirection = eigenVectors(:,2);
        else
            maxMoment = diagMatrix(2,2);
            minMoment = diagMatrix(1,1);
            mainDirection = eigenVectors(:,1);
        end
        
        if(mainDirection(1) < 0)
            mainDirection = - mainDirection;
        end
        
        % a - big semi-axis of the fitted ellipse
        % b - small semi-axis of the fitted ellipse
        % a = sqrt(area/pi*sqrt(maxMoment/minMoment));
        
        %% testing the image
        %         test_frame = zeros(size(frame, 1), size(frame, 2), 3);
        %         test_frame(:,:,1) = anchors_image;
        %         for center_num = 1 : 2
        %             test_frame(round(center_list(point_nums{1}(center_num), 2)), round(center_list(point_nums{1}(center_num), 3)), 2) = 1;
        %             test_frame(round(center_list(point_nums{2}(center_num), 2)), round(center_list(point_nums{2}(center_num), 3)), 3) = 1;
        %         end
        %         figure;
        %         imshow(test_frame);
        
        %% find intersection of the strip and the tissue main inertia axis
        final_points = cell(2, 1);
        for num = 1 : 2
            distance_points = cell(2,1);
            a0 = mainDirection;
            a1 = center_list(point_nums{num}(2), :) - center_list(point_nums{num}(1), :);
            a1 = a1(2:3);
            %a1 = a1./((a1(1)^2 + a1(2)^2)^0.5); % this doesn't change anything
            r0 = center';
            r1 = center_list(point_nums{num}(1),:);
            r1 = r1(2:3);
            coeff = (a1(2)*(r1(1) - r0(1)) - a1(1)*(r1(2) - r0(2)))/(a0(1)*a1(2) - a0(2)*a1(1));
            final_points{num} = r0 + a0*coeff;
            
            %estimated_points_video(round(final_points{num}(1)),round(final_points{num}(2)),frame_num) = 1;
            combined_video(round(final_points{num}(1)),round(final_points{num}(2)), 3,frame_num) = 1;
            combined_video(round(center(1)),round(center(2)), 3,frame_num) = 1;
            for k = 0 : 100
                combined_video(round(r0(1) + coeff*a0(1)*k/100),round(r0(2) + coeff*a0(2)*k/100), 2,frame_num) = 1;
                combined_video(round(r1(1) + a1(1)*k/100),round(r1(2) + a1(2)*k/100), 2,frame_num) = 1;
            end
            
        end
        
        measured_data(frame_num, 1) = frame_num/frame_rate; % time, seconds
        measured_data(frame_num, 2) = final_points{1}(1)*mm_per_px; % anchor point 1, row, mm
        measured_data(frame_num, 3) = final_points{1}(2)*mm_per_px; % anchor point 1, col, mm
        measured_data(frame_num, 4) = final_points{2}(1)*mm_per_px; % anchor point 2, row, mm
        measured_data(frame_num, 5) = final_points{2}(2)*mm_per_px; % anchor point 2, col, mm
        d = sqrt((final_points{1}(1) - final_points{2}(1))^2 + (final_points{1}(2) - final_points{2}(2))^2)*mm_per_px; % distance, mm
        measured_data(frame_num, 6) = d;
        
        if(L_strip <= 0)
            strip_length = getStripLength(strip, final_points)*mm_per_px
            accept_length = input('Press Enter if satisfied with this strip detection or any other key for manual measurement: ');
            if(isempty(accept_length))
                L_strip = strip_length;
            else
                f = figure;
                imshow(combined_video(:,:,:,frame_num));
                [x, y] = getline(f);
                close(f);
                L_strip = 0;
                for i = 2: length(x)
                    L_strip = L_strip + ((x(i) - x(i-1)).^2 + (y(i) - y(i-1)).^2).^0.5;
                end
                L_strip = L_strip*mm_per_px
            end
        end
        
        % finding the bending half-angle using strip length and d
        alpha = fminunc(@(x) abs(L_strip.*integral(@(t) cos(t)./(2.*(cos(t) - cos(x))).^0.5, -x, x) - d.*integral(@(t) (2.*(cos(t) - cos(x))).^(-0.5), -x, x)), pi/2, optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton'));
        %force = (E_strip*I_strip/((d/1000)^2))*integral(@(x) cos(x)./(2.*(cos(x) - cos(alpha))).^0.5, -alpha, alpha)^2
        force = (E_strip*I_strip/((L_strip/1000)^2))*integral(@(x) (2.*(cos(x) - cos(alpha))).^(-0.5), -alpha, alpha)^2;
        measured_data(frame_num, 7) = alpha;
        measured_data(frame_num, 8) = force;
        
        %% plot some graphs
%         l_data = linspace(10, 50, 40);
%         f_data_d1 = zeros(40, 1);
%         alpha_data_d1 = zeros(40, 1);
%         f_data_d2 = zeros(40, 1);
%         alpha_data_d2 = zeros(40, 1);
%         for idx = 1 : 40
%             d = 6;
%             alpha_data_d1(idx) = fminunc(@(x) abs(l_data(idx).*integral(@(t) cos(t)./(2.*(cos(t) - cos(x))).^0.5, -x, x) - d.*integral(@(t) (2.*(cos(t) - cos(x))).^(-0.5), -x, x)), pi/2, optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton'));
%             f_data_d1(idx) = (E_strip*I_strip/((l_data(idx)/1000)^2))*integral(@(x) (2.*(cos(x) - cos(alpha_data_d1(idx)))).^(-0.5), -alpha_data_d1(idx), alpha_data_d1(idx))^2;
%             d = 5.5;
%             alpha_data_d2(idx) = fminunc(@(x) abs(l_data(idx).*integral(@(t) cos(t)./(2.*(cos(t) - cos(x))).^0.5, -x, x) - d.*integral(@(t) (2.*(cos(t) - cos(x))).^(-0.5), -x, x)), pi/2, optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton'));
%             f_data_d2(idx) = (E_strip*I_strip/((l_data(idx)/1000)^2))*integral(@(x) (2.*(cos(x) - cos(alpha_data_d2(idx)))).^(-0.5), -alpha_data_d2(idx), alpha_data_d2(idx))^2;
%             
%         end
%         figure;
%         plot(l_data, f_data_d1);
%         title('f(L)');
%         
%         figure;
%         plot(l_data, alpha_data_d1*180/pi);
%         title('alpha(L)');
%         
%         figure;
%         plot(l_data, f_data_d2./f_data_d1);
%         title('F(d=5.5)/F(d=6)');
        
        %% fitting theoretical strip shape
        %         CC = bwconncomp(tissue);
        %         if(CC.NumObjects ~= 1)
        %             disp('ERROR! More than one object is detected for the tissue');
        %         end
        %             [raw_tensor, center] = getInertiaTensor(tissue);
        %             [eigenVectors, diagMatrix] = eig(raw_tensor);
        %             area = nnz(tissue);
        %
        %             % orientation - angle between the vector and the regular
        %             % x-axis
        %
        %             if diagMatrix(1,1) > diagMatrix(2,2)
        %                 maxMoment = diagMatrix(1,1);
        %                 minMoment = diagMatrix(2,2);
        %                 mainDirection = eigenVectors(:,2);
        %             else
        %                 maxMoment = diagMatrix(2,2);
        %                 minMoment = diagMatrix(1,1);
        %                 mainDirection = eigenVectors(:,1);
        %             end
        %
        %             if(mainDirection(1) < 0)
        %                 mainDirection = - mainDirection;
        %             end
        %
        %             % a - big semi-axis of the fitted ellipse
        %             % b - small semi-axis of the fitted ellipse
        %             a = sqrt(area/pi*sqrt(maxMoment/minMoment));
        %             % first iteration of the fitted curve will begin at point1 and
        %             % end at point2 and have 90 degree bending angle
        %             point1 = center' - mainDirection*a;
        %             point2 = center' + mainDirection*a;
        %
        %             strip_skel = frame & ~tissue;
        %             strip_skel = strip_skel & ~bwmorph(strip_skel, 'erode', 1);
        %             %strip_skel = bwboundaries(frame_strip & ~tissue);
        %             [strip_rows, strip_cols] = find(strip_skel);
        %
        %             if(optimal_conditions(1) == 0)
        %                 alpha = pi/2;
        %                 L = (sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2))*2.18844; % this coefficient corresponds to 90 deree angle
        %                 row0 = point1(1);
        %                 col0 = point1(2);
        %                 phi = atan(mainDirection(2)/mainDirection(1)) + pi/2;
        %                 while(phi < -pi/2)
        %                     phi = phi + pi;
        %                 end
        %                 while(phi > pi/2)
        %                     phi = phi - pi;
        %                 end
        %
        %                 initial_conditions(1) = L;
        %                 initial_conditions(2) = round(L/4);
        %                 initial_conditions(3) = alpha;
        %                 initial_conditions(4) = row0;
        %                 initial_conditions(5) = col0;
        %                 initial_conditions(6) = phi;
        %
        %                 lower_bound = [L_bent_min_mm/mm_per_px, L_end_min_mm/mm_per_px, alpha - pi/6, 0, 0, phi - pi/12];
        %                 upper_bound = [L_bent_max_mm/mm_per_px, L_bent_max_mm/mm_per_px, alpha + pi/18, size(frame_tissue,1), size(frame_tissue,2), phi + pi/12];
        
        %                 optimal_conditions = fmincon(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions, [],[],[],[],lower_bound, upper_bound, [], optimoptions('fmincon', 'Display', 'off'));
        %                 L = optimal_conditions(1);
        %                 l = optimal_conditions(2);
        %                 fitted_strip_parameters{frame_num} = optimal_conditions;
        %             else
        %                 initial_conditions = optimal_conditions(end - 3:end);
        
        %                 lower_bound = [L_bent_min_mm/mm_per_px, L_end_min_mm/mm_per_px, alpha - pi/6, 0, 0, phi - pi/12];
        %                 upper_bound = [L_bent_max_mm/mm_per_px, L_bent_max_mm/mm_per_px, alpha + pi/18, size(frame_tissue,1), size(frame_tissue,2), phi + pi/12];
        
        %                 optimal_conditions = fmincon(@(arg) getFittingError(strip_rows, strip_cols, N, L, l, arg(1), arg(2), arg(3), arg(4)), initial_conditions, [],[],[],[],lower_bound(3:end), upper_bound(3:end), [], optimoptions('fmincon', 'Display', 'off'));
        %
        %                 optimal_conditions = [L, l, optimal_conditions];
        %                 fitted_strip_parameters{frame_num} = optimal_conditions;
        %             end
        %             optimal_conditions = fminunc(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions, optimoptions('fminunc', 'MaxIter', 200, 'MaxFunEvals', 2000));
        %             optimal_conditions = fminsearch(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions);
        %             getFittingError(strip_rows, strip_cols, N, optimal_conditions(1), optimal_conditions(2), optimal_conditions(3), optimal_conditions(4), optimal_conditions(5), optimal_conditions(6))
    end
end
disp(['processing time: ', num2str(round(etime(clock,t0))), ' seconds']);

%combined_video(:,:,1,:) = binary_video;
implay(combined_video);
implay(tissue_video);

figure;
plot(measured_data(1:end-1,1), measured_data(1:end-1,8));
title('Force');
figure;
plot(measured_data(1:end-1,1), measured_data(1:end-1,6));
title('Tissue Length');
figure;
plot(measured_data(1:end-1,1), measured_data(1:end-1,7)*180/pi);
title('Strip Angle');

%% Find strip force and angle from distance data
d = 5;
alpha = fminunc(@(x) abs(L_strip.*integral(@(t) cos(t)./(2.*(cos(t) - cos(x))).^0.5, -x, x) - d.*integral(@(t) (2.*(cos(t) - cos(x))).^(-0.5), -x, x)), pi/2, optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton'));
        force = (E_strip*I_strip/((L_strip/1000)^2))*integral(@(x) (2.*(cos(x) - cos(alpha))).^(-0.5), -alpha, alpha)^2;
        measured_data(frame_num, 7) = alpha;
        measured_data(frame_num, 8) = force;


%% save video file
% [video_name, video_path] = uiputfile({'*.avi'},'Save file as: ', path);
% v = VideoWriter([video_path, video_name]);
% open(v);
% writeVideo(v, combined_video);
% close(v);

%% plot theoretical force vs tissue length
% test_angles = linspace(pi/200, pi*3/4, 150);
% test_tissue_length = zeros(150, 1);
% test_forces = zeros(150,1);
% n = 1;
% for alpha = test_angles
%     M_alpha = integral(@(x) (2.*(cos(x) - cos(alpha))).^(-0.5), -alpha, alpha);
%     N_alpha = integral(@(x) cos(x)./(2.*(cos(x) - cos(alpha))).^0.5, -alpha, alpha);
%     test_tissue_length(n) = L_avg*N_alpha/M_alpha; % distance, mm
%     test_forces(n) = (E_strip*I_strip/((L_avg/1000)^2))*M_alpha^2;
%     n = n+1;
% end
% figure;
% plot(test_tissue_length, test_forces*10^6);

%%
% messed_video = binary_video & ~tissue_video;
% %messed_video = tissue_video;
% %implay(binary_video);
%
% %implay(messed_video, 30);
%
% testFrame = 52;
% frame_strip = messed_video(:,:,testFrame);
% frame_tissue = tissue_video(:,:,testFrame);
%
% [raw_tensor, center] = getInertiaTensor(frame_tissue);
% [eigenVectors, diagMatrix] = eig(raw_tensor);
% area = nnz(frame_tissue);
%
% % orientation - angle between the vector and the regular
% % x-axis
%
% if diagMatrix(1,1) > diagMatrix(2,2)
%     maxMoment = diagMatrix(1,1);
%     minMoment = diagMatrix(2,2);
%     mainDirection = eigenVectors(:,2);
% else
%     maxMoment = diagMatrix(2,2);
%     minMoment = diagMatrix(1,1);
%     mainDirection = eigenVectors(:,1);
% end
%
% if(mainDirection(1) < 0)
%     mainDirection = - mainDirection;
% end
%
% % a - big semi-axis of the fitted ellipse
% % b - small semi-axis of the fitted ellipse
% a = sqrt(area/pi*sqrt(maxMoment/minMoment));
% b = sqrt(area/pi*sqrt(minMoment/maxMoment));
%
% % first iteration of the fitted curve will begin at point1 and
% % end at point2 and have 90 degree bending angle
% point1 = center' - mainDirection*a;
% point2 = center' + mainDirection*a;
%
% strip_skel = frame_strip & ~tissue;
% strip_skel = strip_skel & ~bwmorph(strip_skel, 'erode', 1);
% %strip_skel = bwboundaries(frame_strip & ~tissue);
% [strip_rows, strip_cols] = find(strip_skel);
%
% alpha = pi/2;
% L = (sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2))*2.18844; % this coefficient corresponds to 90 deree angle
% row0 = point1(1);
% col0 = point1(2);
% phi = atan(mainDirection(2)/mainDirection(1)) + pi/2;
% while(phi < -pi/2)
%     phi = phi + pi;
% end
% while(phi > pi/2)
%     phi = phi - pi;
% end
%
% initial_conditions(1) = L;
% initial_conditions(2) = round(L/4);
% initial_conditions(3) = alpha;
% initial_conditions(4) = row0;
% initial_conditions(5) = col0;
% initial_conditions(6) = phi;
%
% N = 50;
%
%
% [rows, cols] =  getFinalStripPoints(N, initial_conditions(1), initial_conditions(2), initial_conditions(3), initial_conditions(4), initial_conditions(5), initial_conditions(6));
%
% cols = cols(round(rows) < size(frame_strip, 1));
% rows = rows(round(rows) < size(frame_strip, 1));
%
% rows = rows(round(cols) < size(frame_strip, 2));
% cols = cols(round(cols) < size(frame_strip, 2));
% figure;
% init_strip = zeros(size(frame_strip));
% init_strip(sub2ind(size(frame_strip), round(rows), round(cols))) = 1;
% full_image = zeros(size(frame_strip, 1), size(frame_strip, 2), 3);
% full_image(:,:,1) = frame_strip | frame_tissue;
% full_image(:,:,2) = init_strip;
% imshow(full_image);
% hold on;
% quiver(center(2), center(1), mainDirection(2)*a, mainDirection(1)*a, '-b', 'linewidth',2);
% ellipse(a,b,atan(mainDirection(1)/mainDirection(2)), center(2), center(1), 'b');
% hold off;
%
% %constraint_matrix = [1,0,0,0,0,0; -1,0,0,0,0,0; 0,1,0,0,0,0; 0,-1,0,0,0,0; 0,0,1,0,0,0; 0,0,-1,0,0,0; 0,0,0,0,0,1; 0,0,0,0,0,-1];
% %constraint_vector = [L + 100; -(L+100)/2; L/2; - L/10; alpha+10; -alpha+30; phi+15; -phi+15];
% %lower_bound = [(L+100)/2, L/10, alpha - pi/6, 0, 0, phi - pi/12];
% %upper_bound = [L+100, L/2, alpha + pi/18, size(frame_tissue,1), size(frame_tissue,2), phi + pi/12];
%
% lower_bound = [L_bent_min_mm/mm_per_px, L_end_min_mm/mm_per_px, alpha - pi/6, 0, 0, phi - pi/12];
% upper_bound = [L_bent_max_mm/mm_per_px, L_bent_max_mm/mm_per_px, alpha + pi/18, size(frame_tissue,1), size(frame_tissue,2), phi + pi/12];
%
% %optimal_conditions = fminunc(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions, optimoptions('fminunc', 'MaxIter', 200, 'MaxFunEvals', 2000));
% %optimal_conditions = fminsearch(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions);
% t0 = clock;
% optimal_conditions = fmincon(@(arg) getFittingError(strip_rows, strip_cols, N, arg(1), arg(2), arg(3), arg(4), arg(5), arg(6)), initial_conditions, [],[],[],[],lower_bound, upper_bound);
% round(etime(clock,t0) * 1000)

%%
% strip_fitting_video = zeros(size(video));
%
% for frame_num = 1 : size(video, 3)
%     %getFittingError(strip_rows, strip_cols, N, optimal_conditions(frame_num, 1), optimal_conditions(frame_num, 2), optimal_conditions(frame_num, 3), optimal_conditions(frame_num, 4), optimal_conditions(frame_num, 5), optimal_conditions(frame_num, 6))
%     if(fitted_strip_parameters{frame_num}(1) > 0)
%         [rows, cols] =  getFinalStripPoints(N, fitted_strip_parameters{frame_num}(1), fitted_strip_parameters{frame_num}(2), fitted_strip_parameters{frame_num}(3), fitted_strip_parameters{frame_num}(4), fitted_strip_parameters{frame_num}(5), fitted_strip_parameters{frame_num}(6));
%
%         cols = cols(round(rows) < size(binary_video, 1));
%         rows = rows(round(rows) < size(binary_video, 1));
%
%         rows = rows(round(cols) < size(binary_video, 2));
%         cols = cols(round(cols) < size(binary_video, 2));
%
%         frame = zeros(size(video, 1), size(video, 2));
%         frame(sub2ind(size(binary_video(:,:,1)), round(rows), round(cols))) = 1;
%         strip_fitting_video(:,:,frame_num) = frame;
%     end
% end
% implay(strip_fitting_video)

%%
%     figure;
%     init_strip = zeros(size(frame_strip));
%     init_strip(sub2ind(size(frame_strip), round(rows), round(cols))) = 1;
%     full_image = zeros(size(frame_strip, 1), size(frame_strip, 2), 3);
%     %full_image(:,:,1) = frame_strip | frame_tissue;
%     full_image(:,:,2) = init_strip;
%     full_image(:,:,3) = strip_skel;
%     imshow(full_image);
%     hold on;
%     quiver(center(2), center(1), mainDirection(2)*a, mainDirection(1)*a, '-b', 'linewidth',2);
%     ellipse(a,b,atan(mainDirection(1)/mainDirection(2)), center(2), center(1), 'b');
%     hold off;

% [wavelet_result, S] = wavedec2(image,1,'db1');
% A = vec2mat(wavelet_result(1,1:numel(image)/4), size(image, 1)/2).';
% B = vec2mat(wavelet_result(1,numel(image)*1/4 + 1 :numel(image)*2/4), size(image, 1)/2).';
% C = vec2mat(wavelet_result(1,numel(image)*2/4 + 1 :numel(image)*3/4), size(image, 1)/2).';
% D = vec2mat(wavelet_result(1,numel(image)*3/4 + 1 :numel(image)*4/4), size(image, 1)/2).';
% 
% A = A/max(A(:));
% B = B/max(B(:));
% C = C/max(C(:));
% D = D/max(D(:));
% 
% figure; imshow([A, B; C, D]);
% 
% level = 2;
% [wavelet_result, S] = wavedec2(image,level,'db1');
% shift = numel(image)/((level*2)^2);
% for i = 0 : level - 1
%     l = level - 1 - i;
%     A = vec2mat(wavelet_result(1,1 + l*shift:(l+1)*shift), floor(size(image, 1)/(l+1)/2)).';
%     B = vec2mat(wavelet_result(1,1 + (l+1)*shift:(l+2)*shift), floor(size(image, 1)/(l+1)/2)).';
%     C = vec2mat(wavelet_result(1,1 + (l+2)*shift:(l+3)*shift), floor(size(image, 1)/(l+1)/2)).';
%     D = vec2mat(wavelet_result(1,1 + (l+3)*shift:(l+4)*shift), floor(size(image, 1)/(l+1)/2)).';
%     figure; imshow([A, B; C, D]);
%     title(['Level ', num2str(l + 1)]);
% end
% 
% 
% 
% A = A/max(A(:));
% B = B/max(B(:));
% C = C/max(C(:));
% D = D/max(D(:));
% 
% figure; imshow([A, B; C, D]);
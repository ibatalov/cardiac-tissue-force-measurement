function [output, tissue_broken] = custom_erode(image, erode_size)

output = image;

%output = bwmorph(output, 'erode', round(erode_size/2));
disk = strel('disk', round(erode_size/2));
output = imerode(output, disk);
tissue_broken = false; % becomes true when tissue break detected

for i = round(erode_size/2) + 1 : erode_size
    output = bwmorph(output, 'erode', 1);
    CC = bwconncomp(output);
    
    if(CC.NumObjects > 1)
        sorted_sizes = zeros(CC.NumObjects, 2);
        for n = 1 : CC.NumObjects
            sorted_sizes(n, :) = [n, numel(CC.PixelIdxList{n})];
        end
        sorted_sizes = sortrows(sorted_sizes, -2);
                
        for n = 3 : CC.NumObjects
            output(CC.PixelIdxList{sorted_sizes(n,1)}) = 0;
        end
        
        if(~tissue_broken)
            [rows1, cols1] = ind2sub(size(output), CC.PixelIdxList{sorted_sizes(1,1)});
            center1 = [sum(rows1(:)), sum(cols1(:))]./sorted_sizes(1,2);
            
            [rows2, cols2] = ind2sub(size(output), CC.PixelIdxList{sorted_sizes(2,1)});
            center2 = [sum(rows2(:)), sum(cols2(:))]./sorted_sizes(2,2);
            
            dir1 = center1 - center2;
            dir1 = dir1/norm(dir1);
            
            [raw_tensor, ~] = getInertiaTensor(output);
            [eigenVectors, diagMatrix] = eig(raw_tensor);
            eig_ratio = sqrt(max(diagMatrix(1,1), diagMatrix(2,2)) / min(diagMatrix(1,1), diagMatrix(2,2)));
            
            if diagMatrix(1,1) > diagMatrix(2,2)
                dir2 = eigenVectors(:,2);
            else
                dir2 = eigenVectors(:,1);
            end
            %figure;
            %imshow(output);
            
            cos_angle = (dir1(:).') * dir2(:); % - dot product = A'*B
            if(abs(cos_angle) > cos(pi/18) && eig_ratio > 5)
                tissue_broken = true;
                %title('Tissue broken');
            else
                output(CC.PixelIdxList{sorted_sizes(2,1)}) = 0;
                %title('Tissue not broken');
            end
        end
    end
end

% figure;
% imshow(output);
% title(['final image: ' num2str(tissue_broken)]);

end
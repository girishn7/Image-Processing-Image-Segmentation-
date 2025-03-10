A = imread('C:\Users\giris\Downloads\shapes4.pnm');
[equalised, hist] = histo(A);
figure;
subplot(2,2,1),bar(hist)
title('Histogram of Image Intensities');
xlabel('Intensity');
ylabel('Frequency');
filter = [1/9, 2/9, 3/9, 2/9, 1/9];
smooth_hist = conv(hist, filter, 'same');
[peaks, locs] = find_custom_peaks(smooth_hist);
[~, idx] = sort(peaks, 'descend');
if length(idx) < 2
    error('Could not find two peaks in the histogram.');
end
peak1 = locs(idx(1));
peak2 = locs(idx(2));
valley = peak1 + ceil((peak2 - peak1)/2);
threshold = valley / 255;
fprintf('Threshold value: %.3f\n', threshold);
binary_A = imbinarize(A, threshold);
subplot(2,2,2)
imshow(binary_A);
title('Binary Image');
label_matrix = zeros(size(binary_A));
label_counter = 0;
for i = 1:size(binary_A, 1)
    for j = 1:size(binary_A, 2)
        if binary_A(i,j) == 1 && label_matrix(i,j) == 0
            label_counter = label_counter + 1;
            stack = [i, j];
            while ~isempty(stack)
                current_pixel = stack(1,:);
                stack(1,:) = [];
                neighbors = [                    current_pixel(1)-1, current_pixel(2);                    current_pixel(1), current_pixel(2)-1;                    current_pixel(1), current_pixel(2)+1;                    current_pixel(1)+1, current_pixel(2)                ];
                for k = 1:size(neighbors,1)
                    if neighbors(k,1) > 0 && neighbors(k,1) <= size(binary_A,1) && neighbors(k,2) > 0 && neighbors(k,2) <= size(binary_A,2)
                        if binary_A(neighbors(k,1),neighbors(k,2)) == 1 && label_matrix(neighbors(k,1),neighbors(k,2)) == 0
                            label_matrix(neighbors(k,1),neighbors(k,2)) = label_counter;
                            stack = [stack; neighbors(k,:)];
                        end
                    end
                end
            end
        end
    end
end
color_map = rand(label_counter,3);
gray_image = label_matrix./max(label_matrix(:));
subplot(2,2,3)
imshow(gray_image);
title('Connected Components');
colorbar;
pseudo_color_image = label2rgb(label_matrix, color_map, 'k');
subplot(2,2,4)
imshow(pseudo_color_image);
title('Pseudo-Colored Connected Components');
props = regionprops(label_matrix, 'Area', 'Centroid', 'Perimeter');
perimeter_vec = cat(1, props.Perimeter);
area_vec = cat(1, props.Area);
compactness = (perimeter_vec .^ 2) ./ (4 * pi * area_vec);
elongation = zeros(length(props), 1);
for i = 1:length(props)
    mu20 = central_moment(binary_A, props(i).Centroid, 2, 0);
    mu02 = central_moment(binary_A, props(i).Centroid, 0, 2);
    mu11 = central_moment(binary_A, props(i).Centroid, 1, 1);
    elongation(i) = (mu20 + mu02 + sqrt((mu20 - mu02)^2 + 4*mu11^2)) / (mu20 + mu02 - sqrt((mu20 - mu02)^2 + 4*mu11^2));
end
for i = 1:length(props)
    fprintf('Component %d:\n', i);
    fprintf('  Area: %d pixels\n', props(i).Area);
    fprintf('  Centroid: (%.2f, %.2f)\n', props(i).Centroid);
    fprintf('  Second order moments: mu20=%.2f, mu02=%.2f, mu11=%.2f\n', mu20, mu02, mu11);
    fprintf('  Perimeter: %.2f pixels\n', props(i).Perimeter);
    fprintf('  Compactness: %.2f\n', compactness(i));
    fprintf('  Elongation: %.2f\n', elongation(i));
end


function [peaks, locs] =  find_custom_peaks(x)
peaks = [];
locs = [];
dx = diff(x);
sign_changes = sign(dx(1:end-1)) ~= sign(dx(2:end));
locs = find(sign_changes);
peaks = x(locs);
end
function m = central_moment(I, c, p, q)
I = double(I);
[x, y] = meshgrid(1:size(I,2), 1:size(I,1));

mu_pq = sum(sum((x-c(1)).^p .* (y-c(2)).^q .* I)) / sum(sum(I));
mu_00 = sum(sum(I));
m = mu_pq / mu_00^(1+(p+q)/2);
end
function [equalized, hist] = histo(A)

hist = zeros(256, 1, 'int32');
for i=1:size(A,1)
    for j=1:size(A,2)
        hist(A(i,j)+1) = hist(A(i,j)+1)+1;
    end
end
cdf = zeros(256, 1);
cdf(1) = hist(1);
for i=2:256
    cdf(i) = cdf(i-1) + hist(i);
end
cdf = cdf / numel(A);
equalized = zeros(size(A), 'uint8');
for i=1:size(A,1)
    for j=1:size(A,2)
        equalized(i,j) = round(cdf(A(i,j)+1)*255);
    end
end
subplot(1,2,1),imshow(A),title('Original image')
subplot(1,2,2),imshow(uint8(equalized)),title('Equalized image after histogram equalization')
hist = int32(hist);

end
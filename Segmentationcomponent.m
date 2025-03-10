A = imread('C:\Users\giris\Downloads\coins.pnm');
equalized = histeq(A);
hist = imhist(equalized);
filter = [1/9, 2/9, 3/9, 2/9, 1/9];
smooth_hist = zeros(size(hist));
for i = 1:length(filter)
    smooth_hist = smooth_hist + filter(i) * circshift(hist, [-(i-1) 0]);
end
[peaks, locs] = findpeaks(smooth_hist);
[~, idx] = sort(peaks, 'descend');
if length(idx) < 2
    error('Could not find two peaks in the histogram.');
end
peak1 = locs(idx(1));
peak2 = locs(idx(2));
valley = peak1 + ceil((peak2 - peak1)/2);
threshold = valley / 255;
binary_A = imbinarize(A, threshold);
label_matrix = bwlabel(binary_A);
color_map = rand(max(label_matrix(:)),3);
gray_image = label_matrix./max(label_matrix(:));
pseudo_color_image = label2rgb(label_matrix, color_map, 'k');
props = regionprops(label_matrix, 'Area', 'Centroid', 'Perimeter');
perimeter_vec = cat(1, props.Perimeter);
area_vec = cat(1, props.Area);
compactness = (perimeter_vec.^ 2)./ (4 * pi * area_vec);
elongation = zeros(length(props), 1);
for i = 1:length(props)
    x = props(i).Centroid(1);
    y = props(i).Centroid(2);
    [rows, cols] = find(label_matrix == i);
    mu20 = sum((cols - x).^2) / props(i).Area;
    mu02 = sum((rows - y).^2) / props(i).Area;
    mu11 = sum((cols - x).*(rows - y)) / props(i).Area;
    elongation(i) = (mu20 + mu02 + sqrt((mu20 - mu02)^2 + 4*mu11^2)) / (mu20 + mu02 - sqrt((mu20 - mu02)^2 + 4*mu11^2));
end
figure;
subplot(2,2,1), bar(hist), title('Histogram of Image Intensities'),xlabel('Intensity'), ylabel('Frequency');
subplot(2,2,2), imshow(binary_A), title('Binary Image');
subplot(2,2,3), imshow(gray_image), title('Connected Components'), colorbar;
subplot(2,2,4), imshow(pseudo_color_image), title('Colored Connected Components');
for i = 1:length(props)
    fprintf('Component %d:\n', i);
    fprintf('  Area: %d pixels\n', props(i).Area);
    fprintf('  Centroid: (%.2f, %.2f)\n', props(i).Centroid);
    fprintf('  Second order moments: mu20=%.2f, mu02=%.2f, mu11=%.2f\n', mu20, mu02, mu11);
    fprintf('  Perimeter: %.2f pixels\n', props(i).Perimeter);
    fprintf('  Compactness: %.2f\n', compactness(i));
    fprintf('  Elongation: %.2f\n', elongation(i));
end

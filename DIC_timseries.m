%% Code to read timeseries of DIC-pseudo Phase images and find lengths and widths
% Created by Tanvi Kale, IISER Pune 
%%------ USER SET PARAMETERS-----
ScalingFactor = 0.045; %scaling factor px to mu
T = 0.70; %binarizing threshold
%% load in the test data 
imgPath = '/Users/tanvikale/Desktop/CA/Bac_cell_size/MG1655_ceph_181023/ROI_2.tif';
numFrame = imfinfo(imgPath); 
disp(['Number of frames found ', num2str(length(numFrame))]); 

outfile = '/Users/tanvikale/Desktop/CA/Bac_cell_size/MG1655_ceph_181023/ROI_2overlay.tif';
%% Show segmentation output of each frame
all_data = table(); 
for f = 1:100  %length(numFrame)
%
inputImage = imread(imgPath, f); 
adjustI = imadjust(inputImage); %imadjust(resultingImage);
%figure(2), imshow(adjustI, []);

% using thresholding and otsu filter on image

% Calculate the histogram
hist = imhist(adjustI);

% Example: Adjust the histogram by stretching the contrast
low_in = 0.55; % Adjust these values as needed
high_in = 0.8;
adjustedHist = hist;
adjustedHist(1:round(low_in*256)) = 0;
adjustedHist(round(high_in*256):256) = 0;

% Calculate the Otsu threshold
threshold = otsuthresh(adjustedHist);

% Binarize the image using the Otsu threshold
binaryImage = imbinarize(adjustI, threshold*0.95);

% Display the binary image
%imshow(binaryImage);
%title('Binary Image (Otsu Threshold)');

% Inverting the image to get cells as forground

%invertedImage = ~binaryImage;
%figure(5),imshow(invertedImage);
%topObjects = bwareafilt(invertedImage, 15);
%figure(6),imshow(topObjects);


% filtering out background 

% Label connected components
labeledImage = bwlabel(binaryImage);%(clearObjects);

% Measure object properties
stats = regionprops(labeledImage, 'Area', 'solidity');

% Define your size threshold
sizeThreshold = 150; % Adjust this value as needed

% Create a binary mask of objects above the size threshold
objectsAboveThreshold = [stats.Area] > sizeThreshold;

% Apply the mask to keep only the objects above the threshold
filteredBinaryImage = ismember(labeledImage, find(objectsAboveThreshold));

% Display the filtered binary image
%figure(4),imshow(filteredBinaryImage);
%title('Objects Above Size Threshold');

% Label connected components
labeledImage = bwlabel(filteredBinaryImage);%(clearObjects);

% Measure object properties
stats = regionprops(labeledImage, 'Area', 'solidity');

% Define your size threshold
sizeThreshold = 5000; % Adjust this value as needed

% Create a binary mask of objects above the size threshold
objectsAboveThreshold = sizeThreshold > [stats.Area];

% Apply the mask to keep only the objects above the threshold
filteredBinaryImage1 = ismember(labeledImage, find(objectsAboveThreshold));

% filtering objects on the basis of solidity

% Label connected components (objects)
labeledImage1 = bwlabel(filteredBinaryImage1);

% Measure region properties, including solidity
stats1 = regionprops(labeledImage1, 'Solidity');

% Set a threshold for solidity
solidityThreshold = 0.5; % Adjust this threshold as needed

% Create a binary mask to filter objects based on solidity
filteredBinaryImage2 = ismember(labeledImage1, find([stats1.Solidity] >= solidityThreshold));

% Display the filtered binary image
%imshow(filteredBinaryImage1);
%title('Filtered Binary Image based on Solidity');

% filling holes in connected objects

BW = imfill(filteredBinaryImage2, 'holes');
%figure(5),imshow(BW);

% detecting edges of the binary objects


% Apply Canny edge detection
edgeImage = edge(BW, 'Canny');

% Display the edge-detected image
%figure(2),imshow(edgeImage);
%title('Binary Image with Outlines (Edges)');

binoverlay = imoverlay(adjustI,edgeImage,'red');
%figure(3),imshow(binoverlay);

% refine contours
act_cont  = activecontour(adjustI, BW, 30,'edge', 'SmoothFactor', 1.0, 'ContractionBias', -0.1); 


label_im = bwlabel(act_cont); 
%bwconn_ = bwconncomp(act_cont);
%region_conn = regionprops(bwconn_, 'PixelIdxList', 'PixelList'); 
%C = labeloverlay(imadjust(adjustI), label_im, 'Colormap', 'winter');

% filtering out by size (to remove single pixels detected by active contours)

% Label connected components
labeledImage = bwlabel(act_cont);%(clearObjects);

% Measure object properties
stats = regionprops(labeledImage, 'Area', 'solidity');

% Define your size threshold
sizeThreshold = 200; % Adjust this value as needed

% Create a binary mask of objects above the size threshold
objectsAboveThreshold = [stats.Area] > sizeThreshold;

% Apply the mask to keep only the objects above the threshold
finalfilteredBinaryImage = ismember(labeledImage, find(objectsAboveThreshold));
label_im = bwlabel(finalfilteredBinaryImage); 
bwconn_ = bwconncomp(finalfilteredBinaryImage);
region_conn = regionprops(bwconn_, 'PixelIdxList', 'PixelList'); 
C = labeloverlay(imadjust(adjustI), label_im, 'Colormap', 'winter');

    figure(2), imshow(C, 'Border', 'tight'), hold on 
    for l = 1:length(region_conn)
        pixel_id = region_conn(l).PixelList;
        conn_id = region_conn(l).PixelIdxList; 
        %C = insertText(C, pixel_id(end,:), num2str(label_im(conn_id(end)))); 
        text(pixel_id(end,1), pixel_id(end, 2), num2str(label_im(conn_id(end))),...
        'FontSize', 16, 'Color', [1,1,1]);
    end

conn_comps = bwconncomp(finalfilteredBinaryImage, 8); 
region_connected = regionprops(conn_comps,'Area', 'PixelList', 'PixelIdxList', ...
    'Centroid');
minimum_area =  90; 
remove_objs = find((cat(1, region_connected.Area) < minimum_area)); 
region_connected(remove_objs) = [];
final_contours = logical(zeros(size(finalfilteredBinaryImage))); 
for obj  = 1:length(region_connected)
    indices = region_connected(obj).PixelIdxList;
    final_contours(indices) = 1; 
end 
clean_overlay = show_overlay(adjustI, finalfilteredBinaryImage); 
%clean_overlay = imoverlay(inputImage, bwperim(final_contours), 'red'); 

%%fig = figure(1); imshow(clean_overlay, 'Border','tight');commented on
%%13/10/2023

% Enlarge figure to full screen.
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%g = gcf;
%g.WindowState = 'maximized'
%figure(2), imshow(C, 'Border', 'tight'), hold on 
%for l = 1:length(region_conn)
 %   pixel_id = region_conn(l).PixelList;
 %   conn_id = region_conn(l).PixelIdxList; 
        %C = insertText(C, pixel_id(end,:), num2str(label_im(conn_id(end)))); 
 %  text(pixel_id(end,1), pixel_id(end, 2), num2str(label_im(conn_id(end))),...
  %      'FontSize', 16, 'Color', [1,1,1]);
%end
%figure(2), imshow(C, 'Border', 'tight');


%cont_overlay = imoverlay(imadjust(adjustI), bwperim(act_cont), [1 0 0]); 
%figure(1), imshow(cont_overlay);
h=getframe(gcf); % converting frame to matrix
imwrite(h.cdata, outfile, 'tif', 'Compression', 'none',...
        'WriteMode', 'append');%, '

return_table = getLengthsAndWidths2(finalfilteredBinaryImage, inputImage); 
return_table.Frame = repmat(f, height(return_table),1);
%{repeatedArray = repelem(f, length(region_conn)); repeatedArray = reshape(repeatedArray, [length(region_conn), 1]);return_table.FrameNo = cat(1, repeatedArray);
%}
return_table.Centroid = cat(1, region_connected.Centroid);
newColumn = (1:length(region_conn))';
return_table.Cellno = cat(1, newColumn);
all_data= [all_data; return_table]; 


end 

%% saving the movie

outFolder = '/Users/tanvikale/Desktop/CA/Bac_cell_size/MG1655_ceph_181023/'; 
filename  =  '/Users/tanvikale/Desktop/CA/Bac_cell_size/MG1655_ceph_181023/ROI_1_overlay.tif';
all_data.Area = ones(height(all_data), 1); 
[tracks, obno, outmat] = DICOT_tracking(outFolder, all_data, 10, 1, 'um', 's'); 
finfo = imfinfo(filename); 
makemovie(outFolder, 1:length(finfo), filename, outmat)

%% Required functions 
function return_table = getLengthsAndWidths2(final_contours2, I)
return_table = array2table(zeros(0, 3)); 
return_table.Properties.VariableNames = {'ObjCoordinates', 'Length', 'Width'}; 

conn_comps = bwconncomp(final_contours2, 4); 
region_connected = regionprops(conn_comps, 'PixelList', 'PixelIdxList');
%L = bwlabel(final_contours2); 
%figure(1), subplot(2,2,4), imshow(crop_image); hold on 
%length_array = []; 
%width_array = []; 
for im = 1:length(region_connected)
    %%
    im_coordinates = region_connected(im).PixelList; 
    %im_indexes = region_connected(im).PixelIdxList; 
    reducedImage = isolateObject(im_coordinates,I);
    %%
    [~, lengths, widths] = get_mid_array(reducedImage);  
    %plot(mid_array(:,2)+ min(im_coordinates(:,1)) - 1, ...
        %mid_array(:,1) + min(im_coordinates(:,2)) - 1 , 'r-', 'MarkerSize', 1);
    %  Compute lengths 
    curr_data = {region_connected(im).PixelIdxList,lengths, widths}; 
    return_table = [return_table; curr_data]; 
    %length_array = [length_array; lengths]; 
    %width_array = [width_array; widths*2]; 
    if im == 82
        disp([lengths, widths]);
    end
end
%% Required functions 
function [mid_points, length_total, width_value] = get_mid_array(reducedImage)
%%
get_angle = regionprops(reducedImage, 'Orientation'); 
rotatetheImage = imrotate(reducedImage, -(get_angle.Orientation-90), 'nearest', 'loose'); 
rotatetheImage(rotatetheImage >= 1) = 1; 
size_image_f = size(rotatetheImage); 

[row, col] =find(rotatetheImage); 
[~, major_axis] = max([length(unique(row)), length(unique(col))]); 
length_dimension = major_axis; % If 1 then go through rows else columns 
mid_array = []; 
for n = 1:max(size_image_f)
  get_all_entries = find(rotatetheImage(n,:)); 
   if ~isempty(get_all_entries)
       mid_point =  get_all_entries(1) + (get_all_entries(end) - get_all_entries(1))/2;
       if isequal(major_axis, 1)
        mid_array = [mid_array;n, mid_point]; 
       else 
           mid_array = [mid_array; mid_point, n];
       end 
   else 
       continue
   end 
end 
%% 
output = fit(mid_array(:,1), mid_array(:,2), 'smoothingspline',...
    'SmoothingParam',0.2); 
mid_array(:,2) = smooth(mid_array(:,2), 0.3); 
mid_array(:,2)  = output(mid_array(:,1)); 
mid_points = mid_array; 
length_segment = hypot(diff(mid_array(:,1)), diff(mid_array(:,2))); 
length_total =  sum(length_segment) * 0.045;  
% get perimeter contours

[p_row, p_col] = find(bwperim(rotatetheImage)); 
D = pdist2(mid_array, [p_row, p_col], 'euclidean'); 

[min_distance, ~] =  min(D, [], 2); 
width_value = max(min_distance)

% for All cells other than A22 treated cells use median values
%median_value = median(min_distance) ; 
%width_value = mean(min_distance(min_distance >= median_value)); 

width_value = width_value*2* 0.045;
 %figure(1), histogram(min_distance,6), hold on 
 %figure(1), line([median_value, median_value], ylim, 'LineWidth', 2.0, 'Color', 'k');
 %figure(1), line([width_value,width_value], ylim, 'LineWidth', 2.0, 'Color', 'b')
 %hold on 
 
%{
 figure(4), imshow(rotatetheImage, [],'InitialMagnification', 5000); hold on 
plot(mid_array(:,2), mid_array(:,1)  , 'r-', ...
    'MarkerSize', 1, 'LineWidth', 2.5);
 for w = 1:length(mid_array)
     first_point = mid_array(w,:); 
     second_point = [p_row(indices(w)), p_col(indices(w))];     
     plot([mid_array(w,2) p_col(indices(w))], ...
         [mid_array(w,1) p_row(indices(w))], 'b--', 'LineWidth', 2.5); 
 end 
%}

%%%figure(4), subplot(2,2,2), imshow(rotatetheImage)
end 

function reducedImage = isolateObject(pixel_coordinates, I) 
%test_im  = I; 
rescale_dims = [pixel_coordinates(:,1) - min(pixel_coordinates(:,1)) + 1, ...
    pixel_coordinates(:,2) - min(pixel_coordinates(:,2)) + 1]; 
null_image = zeros(max(rescale_dims(:,2)), max(rescale_dims(:,1))); 
burn_indx = sub2ind(size(null_image), rescale_dims(:,2), rescale_dims(:,1)); 
null_image(burn_indx) = 1;
reducedImage = null_image; 

%width = max(pixel_coordinates(:,1)) - min(pixel_coordinates(:,1));
%height = max(pixel_coordinates(:,2)) - min(pixel_coordinates(:,2));
%reducedImage = contour_original; 
% figure(2) ,imshow(null_image); 
end 

end

%%
function im =  show_overlay(im1, im2)
im2 = bwperim(im2);
im2 = bwlabel(im2); 
im = labeloverlay(im1, im2, "Colormap",'autumn'); 
end 


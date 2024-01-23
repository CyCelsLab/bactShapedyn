% Cephalexin + A22 bulge data: Finding point of inflection: 
% close all; 
clear; clc;
% Importing full data: 
length_data = readmatrix("Ryth_CephA22_length.csv");
width_data = readmatrix("Ryth_CephA22_width.csv");
tend = length(width_data);
time_data = (0:tend-1)';
figure()
for i = 1:10
    subplot(2,5,i)
    cellwidth = width_data(:,i);
    inflectionList(i) = findInflection(cellwidth);
    hold on 
    plot(time_data, cellwidth)
    xline(time_data(inflectionList(i)))
    title(sprintf('Cell %i',i))
    ylim([0.5, 2])
    xlim([0,100])
    hold off
end
sgtitle('Finding inflection point:')

% Some inflections are incorrectly determined, Manually fixing them: 

% Cell 6 and 10 don't look right. As a result, we manually change them.
% Since inflectionList has the index values of the inflection point, we
% need to find the index values of the inflection points of 6 and 10 manually 

% inflection point of 6: t = 69 min
% inflection point of 10: t = 49 min

% inflection Point index for 6: 
inflectionIndex6 = find(time_data == 69);

% inflection point index for 10: 
inflectionIndex10 = find(time_data==49);

% replace the values in inflectionList to make it more accurate: 
inflectionList(6) = inflectionIndex6;
inflectionList(10) = inflectionIndex10;

clear inflectionIndex10;
clear inflectionIndex6

figure()
for i = 1:10
    subplot(2,5,i)
    hold on 
    cellwidth = width_data(:,i);
    plot(time_data, cellwidth)
    xline(time_data(inflectionList(i)))
    title(sprintf('Cell %i',i))
    ylim([0.5, 2])
    xlim([0,100])
    hold off
end
sgtitle('Finding inflection point (MANUALLY FIXED):')

%% Overlapping the width: Keeping inflection at t = 0

allignedMatrix = NaN(10, 2*tend+1); % dimensions of a NaN matrix

% t = t-inflection; 
figure()
for i = 1:10
    cellwidth = width_data(:,i);
    inflectionTime = time_data(inflectionList(i));
    hold on 
    plot(time_data-inflectionTime, cellwidth)
    xline(0)
    for t = 1:tend
        allignedMatrix(i,tend+t-inflectionList(i))= cellwidth(t);
    end
end
writematrix(allignedMatrix, 'pre-post_allignedMatrix_width.csv')

%% Checking if AllingedMatrix looks correct
% figure()
% for i = 1:10
%     j = 1:2*tend+1;
%     hold on 
%     plot(j - 90, allignedMatrix(i, :))
% end

%% Mean values after allignement 
for j = 1:length(allignedMatrix)
    meanWidth(j) = nanmean(allignedMatrix(:, j));
    stdWidth(j) = nanstd(allignedMatrix(:,j));
end

figure()
j = 1:2*tend+1;
shiftedTime = j-90;
hold on 
errorbar(shiftedTime, meanWidth, stdWidth, 'linewidth',2, 'Color','#c9c9c9', 'CapSize',0);
for i = 1:10
    j = 1:2*tend+1;
    plot(shiftedTime, allignedMatrix(i, :))
end
xline(0, 'LineStyle','--')
plot(shiftedTime, meanWidth, 'k', 'LineWidth',2)
% At meanWidth(90), the bulge starts. 

%% For lengths: 

figure()
for i = 1:10
    subplot(2,5,i)
    hold on 
    celllength = length_data(:,i);
    plot(time_data, celllength)
    xline(time_data(inflectionList(i)))
    title(sprintf('Cell %i',i))
%     ylim([0.5, 2])
    xlim([0,100])
    ylabel('Length')
    hold off
end
sgtitle('Finding inflection point (MANUALLY FIXED):')


allignedMatrix = NaN(10, 2*tend+1); % dimensions of a NaN matrix
figure()
for i = 1:10
    celllength = length_data(:,i);
    inflectionTime = time_data(inflectionList(i));
    hold on 
    plot(time_data-inflectionTime, celllength)
    xline(0)
    for t = 1:tend
        allignedMatrix(i,tend+t-inflectionList(i))= celllength(t);
    end
end

writematrix(allignedMatrix, 'Pre-post_alligned_matrix_length.csv')

for j = 1:length(allignedMatrix)
    meanLength(j) = nanmean(allignedMatrix(:, j));
    stdLength(j) = nanstd(allignedMatrix(:,j));
end

figure()
j = 1:2*tend+1;
shiftedTime = j-90;
hold on 
errorbar(shiftedTime, meanLength, stdLength, 'linewidth',2, 'Color','#c9c9c9', 'CapSize',0);
for i = 1:10
    j = 1:2*tend+1;
    plot(shiftedTime, allignedMatrix(i, :))
end
xline(0, 'LineStyle','--')
plot(shiftedTime, meanLength, 'k', 'LineWidth',2)
% At meanWidth(90), the bulge starts. 

%% Saving pre-bulge and post-bulge timelines separately: 

pre_bulge(:,1) = meanLength(1:89);
post_bulge(:,1) = meanLength(90:end);

pre_bulge(:,2) = meanWidth(1:89);
post_bulge(:,2) = meanWidth(90:end);

writematrix(pre_bulge,'pre_bulge.csv')
writematrix(post_bulge, 'post_bulge.csv')

%% findInflection function definition: 

% Values for dt, Dt and DT are changed until the algorithm finds the inflection
% point correctly. dt = 10 works well within the data we have. Others don't
% matter
function inflectionLoc = findInflection(CellWidthArray)

    dt = 10;
    Dt = 1;
    DT = 1;
    
    smoothedWidth = movmean(CellWidthArray,dt);
    grad1 = gradient(smoothedWidth,Dt);
    
    
    % Taking the second gradient of the width data which will give information
    % on the speed of grad changing. We want the location where grad2 is max. 
    grad2 = gradient(grad1, DT);
    [~, globalPeakLoc] = max(grad2);
    
    % Visualization
%         figure(1)
%         hold on 
%         plot(1:tend, cellwidth)
%         plot(1:tend, smoothedWidth, 'LineWidth',2)
%         plot(1:tend, gradient(grad,0.05))
%         plot(1:tend,grad)
%         xline(globalPeakLoc)
%         hold off

    inflectionLoc = globalPeakLoc;

end

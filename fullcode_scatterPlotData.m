% This code contains the full pipeline for particle swarm optimisation. for
% all data points instead of the mean. 
% Some functions listed below will need to have to be in the same directory
% as this code to run. 

% Functions in the same directory: 
%   - Dinner_ODEs.m
clc; clear;

%% Data import for SJ wildtype
% Input data format: 
% [data_time], [avg_length], [avg_radius] 
% All above should be row vectors of the same length
tic;

%% Data import for perturbations
figure()

for i = 1:10

% For a22:
% filename = 'MG1655_a22_5cells_dynamics';
% data = readmatrix(strcat(filename,'.csv'));
% data_time = data(1:85,1);
% avg_length = data(1:85,2*i);
% avg_radius = 0.5*data(1:85,2*i+1);

% For Ceph: 
% filename = 'MG1655_ceph_lw_11cell_dynamics';
% data = readmatrix(strcat(filename,'.csv'));
% data_time = data(:,1);
% avg_length = data(:,2*i)+0.3;
% avg_radius = (0.5*data(:,2*i+1)+0.3);



% For Wildtype:
% filename = 'wildtype_lw';
% data = readmatrix(strcat(filename,'.csv'));
% data_time = data(:,1);
% avg_length = data(:,2*i);
% avg_radius = 0.5*data(:,2*i+1);

% For pre bulge: 
% filename = 'pre_bulge';
% filename1 = 'preBulge_length';
% data_length = readmatrix(strcat(filename1,'.csv'));
% filename2 = 'preBulge_width';
% data_width = readmatrix(strcat(filename1,'.csv'));
% avg_length = data_length(:,i);
% avg_radius = 0.5*data_width(:,i);


% For post bulge: 
filename = 'post_bulge';
filename1 = 'postBulge_length';
data_length = readmatrix(strcat(filename1,'.csv'));
filename2 = 'postBulge_width';
data_width = readmatrix(strcat(filename1,'.csv'));
avg_length = data_length(:,i);
avg_radius = 0.5*data_width(:,i);


% Clear up the NaN values; 
NanIndex= find(~isnan(avg_length)); % only once, since they are the same for both width and length

% Remove the Nan values; 
avg_length = avg_length(NanIndex(1):NanIndex(end));
avg_radius = avg_radius(NanIndex(1):NanIndex(end));

data_time = 0:length(avg_radius)-1;

% Clear variables that arent required
clear NanIndex;


tspan = [0,length(avg_length)+1]; % The time for which the simulation calculates radius and length

% Renaming lists to make it easier to type and be consistent with
% previously written code. 
avg_L = avg_length; 
avg_r = avg_radius;

% params = [gamma_n; k_n; eta_L; eta_r; P; R0];
opt_params = [];

% X0 - initial r, initial L
X0 = [avg_r(1), avg_L(2)];

% %%%%%%%%%%%
% ONLY FOR PRE_BULGE, CHANGING INTIAL CONDITIONS MANUALLY: 
% X0 = [avg_r(1), 1.7];
% %%%%%%%%%%%%%


% params = [gamma_n; k_n; eta_L; eta_r; P; R0];
% lb = [0.1, 4.5, 50, 50, 0.1, 0.1];
% ub = [0.5, 7, 10000, 10000, 0.7, 0.4];

% Iteration 1
% lb = [0.2, 4.5, 1500, 4000, 0.5, 0.3619809063];
% ub = [0.4, 13, 1500, 4000, 0.5, 0.3619809063];

% 
% % Iteration 2
% lb = [0.0517, 4.5, 50, 50, 0.6860, 0.3572];
% ub = [0.0517, 13, 10000, 10000, 0.6860, 0.3572];

% For Perturbations: 
lb = [0.1, 2.5, 500, 500, 0.2534928964, 0.1128442326];
ub = [0.7, 7, 10000, 10000, 0.2534928964, 0.1128442326];

options = optimoptions('particleswarm', 'UseParallel',true);
[opt_params,fval] = particleswarm(@(params) ObjectiveFunction(params, X0, data_time, avg_length, avg_radius, tspan), 6, lb, ub, options);

%%%% ONLY FOR POST_BULGE%%%%%%%%%%%%%%
% opt_params(4) = opt_params(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


writematrix(opt_params, sprintf('NonEqualVariable_etas-constant_P&R0_Scatter_Opt_Params_%s_%i.csv', filename, i))

%% Particle Swarm output figure: 

hold on 
X0 = [avg_radius(1), avg_length(1)];

% %%%%%%%%%%%
% ONLY FOR PRE_BULGE, CHANGING INTIAL CONDITIONS MANUALLY: 
% X0 = [avg_r(1), 1.7];
% %%%%%%%%%%%%%


[t,y] = ode15s(@(t,y) Dinner_ODEs(t, y, opt_params),tspan, X0);
subplot(1,2,1)
hold on
set(gca, 'Fontsize',15)
plot(t, 2*y(:,1), 'linewidth',2,'Color','b', 'DisplayName','Width simulation')
scatter(data_time, 2*avg_r,'MarkerEdgeColor', 'k','Marker', 'x', 'DisplayName', 'Data') % Scatter plot of time vs width = 2*radius
xlabel('Time (min)')
ylabel('Width (\mum)')
legend('Location','northwest')
ylim([0, 3])
hold off 

subplot(1,2,2)
hold on
set(gca, 'Fontsize',15)
plot(t, y(:,2), 'linewidth',2,'Color','r', 'DisplayName','Length simulation')
scatter(data_time, avg_L,'MarkerEdgeColor', 'k', 'Marker', 'x', 'DisplayName', 'Data') % Scatter plot of time vs width = 2*radius
xlabel('Time (min)')
ylabel('Length (\mum)')
legend('Location','northwest')
converted_params = opt_params; 
converted_params(1) = opt_params(1)*opt_params(5)*opt_params(6);
converted_params(2) = opt_params(2)*opt_params(5)*opt_params(6)^3;
text(25,3, sprintf('\\gamma = %0.3f, k = %0.3f, \\eta_L = %0.3f,\n \\eta_r = %0.3f, P = %0.3f, R_0 = %0.3f', converted_params), 'FontSize',15)


sgtitle(strrep(filename, '_', ' '), fontsize= 20)
end
disp('-----')

for k = 1:5
    beep()
    pause(.25)
end

toc;
%% Output parameters (Opt_params) in csv
% writematrix(opt_params, sprintf('8Jan_Opt_Params_%s_var_eta.csv', filename))

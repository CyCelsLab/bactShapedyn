% Defining the objective function for data fitting: 
function obj = ObjectiveFunction(params, X0, data_time, avg_length, avg_radius, tspan)
    % Solving the ODEs 
        [t,y] = ode15s(@(t,y) Dinner_ODEs(t, y, params),tspan, X0) ; 
        
        % Interpolation of the solution of the ODEs so that the L2 norm can be
        % found between the ODE solution and all the data points at that time
        % point. 
        
        % using vq_length = interp1(t,y(:,1),query_time)
         
        % Defining the objective function:
        obj_width = norm(interp1(t,y(:,1),data_time) - avg_radius);
        obj_length = norm(interp1(t,y(:,2),data_time) - avg_length);
        obj = obj_width*obj_length;


        if isnan(obj)
            fprintf('Objective is NaN, setting value to 1000000\n')
            obj = 1000000;
        end
       


end
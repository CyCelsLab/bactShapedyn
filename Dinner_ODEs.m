function dXdt = Dinner_ODEs(t, X, params)

    % Initializing variables 
    r = X(1); 
    L = X(2); 

    % Initializing parameters
%     PI = 0.3; % Turgor pressure in E.coli (nN um^-2)
%     PI_a = 0.4; % Growth pressure in E.coli (nN um^-2)
%     R0 = 0.38; % Preferred radius of cross-section (um)
    h = 3*10^(-3); % Cell wall thickness (um)

    gamma_n = params(1);
    k_n = params(2);
    eta_L = params(3);
    eta_r = params(4); 
    P = params(5);
    R0 = params(6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONLY FOR POST_BULGE
%     eta_r = eta_L;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Defining other useful values based on the parameters
    mu_L = 1/(2*pi*h*eta_L); % page 4 
    mu_r = 1/(2*pi*h*eta_r); % page 4
%     P = PI + PI_a; % PI_a is the active growth pressure
    gamma = gamma_n*(P*R0);
    k = k_n*P*R0^3;


    % Defining Energy density and its derivatuve with respect to r
    U = -P*pi.*(r.^2) + 2*pi*gamma.*r + k*pi.*r.*(1./r - 1/R0).^2;
    dUdr = -2*P*pi.*r + 2*pi*gamma - k*pi*(1./(r.^2)-1/(R0^2));
    kappa = -mu_L*U/r;
    
    % Defining the ODE equations: eqn 7a and 7b from Banerjee et al 2016.
    dLdt = L*kappa;
    drdt = -mu_r*r*dUdr;


    dXdt = [drdt;dLdt];

end

function  [a, b, c, d, e, t1, t2, t3, t4, t5, purity, recovery, productivity, energy_requirments, TCR] = PSACyclesample(L, P_0, P_l, t_ads, v_0, alpha, beta, N, type, SimParam, IsothermPar)
    ndot_0      =   P_0*v_0/(8.314*313.15)  ;   % Inlet molar flux [mol/s/m^2]
    t_pres      = 20                    ;   % Maximum/time of pressurization step [s]
    t_CnCdepres = 30                    ;   % Maximum/time of depressurization step [s]
    t_CoCdepres = 70                    ;   % Maximum/time of depressurization step [s]
    t_LR        = t_ads                 ;   % Time of light reflux step [s]
    t_HR        = t_LR                  ;   % Time of heavy reflux step [s]
    tau         = 0.5                   ;   % Parameter used for determining speed of pressure change
    P_inlet     = 1.02                  ;   % Pressure of feed gas at the inlet of the adsorption step

    % Flue gas parameters and constants
    R          = 8.314                  ;   % Universal gas constant [J/mol/K : Pa*m^3/mol/K]
    T_0        = 313.15                 ;   % Feed temperature of flue gas [K]
    y_0        = 0.15                   ;   % Inlet gas CO2 mole fraction[-]
    Ctot_0     = P_0/R/T_0              ;   % Inlet total concentration [mol/m^3]
    %v_0        = ndot_0/Ctot_0          ;   % Inlet velocity and scaling parameter [m/s]
    mu         = 1.72e-5                ;   % Viscosity of gas [Pa*s]
    epsilon    = 0.37                   ;   % Void fraction
    D_m        = 1.2995e-5              ;   % Molecular diffusivity [m^2/s]
    K_z        = 0.09                   ;   % Thermal conduction in gas phase [W/m/k]
    C_pg       = 30.7                   ;   % Specific heat of gas [J/mol/k]
    C_pa       = 30.7                   ;   % Specific heat of adsorbed phase [J/mol/k]
    MW_CO2     = 0.04402                ;   % Molecular weight of CO2 [kg/mol]
    MW_N2      = 0.02802                ;   % Molecular weight of N2 [kg/mol]
    %feed_gas  = 'Constant Pressure'    ;   % Whether flue gas during the feed step has a constant pressure or velocity
    feed_gas   = 'Constant Velocity'    ;   % Whether flue gas during the feed step has a constant pressure or velocity

    % Adsorbent parameters
    ro_s        = SimParam(type,1)      ;   % Density of the adsorbent [kg/m^3]
    r_p         = 1e-3                  ;   % Radius of the pellets [m]
    C_ps        = 1070                  ;   % Specific heat capacity of the adsorbent [J/kg/K]
    q_s         = IsothermPar(type,7)   ;   % Molar loading scaling factor [mol/kg]
    q_s0        = q_s*ro_s              ;   % Molar loading scaling factor [mol/m^3]
    k_CO2_LDF   = 0.1631                ;   % Mass transfer coefficient for CO2 [1/s]
    k_N2_LDF    = 0.2044                ;   % Mass transfer coefficient for N2 [1/s]

    % Isotherm parameters
    q_s_b      = [IsothermPar(type,1), IsothermPar(type,7)]   ;   % Saturation loading on site b [mol/kg]
    q_s_d      = [IsothermPar(type,2),  IsothermPar(type,8)]   ;   % Saturation loading on site d [mol/kg]
    b1          = [IsothermPar(type,3),  IsothermPar(type,9)]   ;   % Pre-exponential factor for site b [Pa-1]
    d1          = [IsothermPar(type,4),  IsothermPar(type,10)]  ;   % Pre-exponential factor for site d [Pa-1]
    deltaU_b   = [IsothermPar(type,5),  IsothermPar(type,11)]  ;   % Heat of adsorption for site b [J/mol]
    deltaU_d   = [IsothermPar(type,6),  IsothermPar(type,12)]  ;   % Heat of adsorption for site d [J/mol]
    
    deltaU = [SimParam(type,2), SimParam(type,3)];

    %deltaU     = [-36000, -15800]     ;

    %Distribute the values to the necessary variables
    Params     = zeros(39, 1) ;
    Params(1)  = N			  ;
    Params(2)  = deltaU(1)    ;
    Params(3)  = deltaU(2)    ;
    Params(4)  = ro_s		  ;
    Params(5)  = T_0		  ;
    Params(6)  = epsilon	  ;
    Params(7)  = r_p		  ;
    Params(8)  = mu			  ;
    Params(9)  = R			  ;
    Params(10) = v_0		  ;
    Params(11) = q_s0		  ;
    Params(12) = C_pg		  ;
    Params(13) = C_pa		  ;
    Params(14) = C_ps		  ;
    Params(15) = D_m		  ;
    Params(16) = K_z		  ;
    Params(17) = P_0		  ;
    Params(18) = L			  ;
    Params(19) = MW_CO2		  ;
    Params(20) = MW_N2		  ;
    Params(21) = k_CO2_LDF	  ;
    Params(22) = k_N2_LDF	  ;
    Params(23) = y_0		  ;
    Params(24) = tau		  ;
    Params(25) = P_l		  ;
    Params(26) = P_inlet	  ;
    Params(27) = 1			  ;   % Place for y at outlet of Adsorption = y at inlet of Light Reflux: y_LP
                                  % y_LR = 1 - No initial guess for inlet CO2 mole fraction in Light Reflux step
    Params(28) = 1			  ;   % Place for T at outlet of Adsorption = T at inlet of Light Reflux: T_LP
                                  % T_LR = 1 - No initial guess for inlet temperature in Light Reflux step
    Params(29) = 1			  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of Light Reflux
                                  % ndot_LR = 1 - No initial guess for inlet ndotin Light Reflux step
    Params(30) = alpha    	  ;
    Params(31) = beta         ;
    
    %Params(32) = P_I          ;
    Params(33) = y_0          ;   % Place for y at outlet of CnC depressurization = y at inlet of of Heavy Reflux: y_HP
                                  % y_HR = y_0 - Initial guess for inlet CO2 mole fraction in Heavy Reflux step
    Params(34) = T_0          ;   % Place for T at outlet of CnC depressurization = T at inlet of of Heavy Reflux: T_HP
                                  % T_HR = T_0 - Initial guess for inlet temperature in Heavy Reflux step
    Params(35) = ndot_0*beta  ; % 300/30   % Place for ndot at outlet of CnC depressurization = ndot at inlet of Heavy Reflux
                                  % ndot_HR = ndot_0*beta - Initial guess for inlet ndotin Heavy Reflux step
    Params(36) = 0.01    	  ;   % Place for y at outlet of Adsorption = y at inlet of CnC pressurization: y_LP
                                  % y_LR = 0.01 - Initial guess for inlet CO2 mole fraction in CnC pressurization step
    Params(37) = T_0    	  ;   % Place for T at outlet of Adsorption = T at inlet of CnC pressurization: T_LP
                                  % T_LR = T_0 - Initial guess for inlet temperature in CnC pressurization step
    Params(38) = ndot_0  	  ;   % Place for ndot at outlet of Adsorption = ndot at inlet of CnC pressurization 
                                  % NOTE: not used, seems not necessary. ndot_LR = ndot_0 - Initial guess for inlet ndot
                                  % in CnC pressurization step
    y_HR            = Params(33)      ;
    T_HR            = Params(34)      ;
    ndot_HR         = Params(35)	  ;                            
      
                                  
    if strcmpi(feed_gas, 'Constant Pressure') == 1
        Params(end) = 1 ;
    elseif strcmpi(feed_gas, 'Constant Velocity') == 1
        Params(end) = 0 ;
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end

    Times          = [ t_pres; t_ads; t_CnCdepres; t_LR; t_CoCdepres; t_HR ] ;

    IsothermParams = [q_s_b, q_s_d, b1, d1, deltaU_b, deltaU_d, IsothermPar(type,13)] ;   
    
    %% Economic Parameters
    
    desired_flow                          = 100          ;   % Desired flow rate of flue gas per column [mol/s]
    electricity_cost                      = 0.07         ;   % Cost of electricity [$/kWh]
    hour_to_year_conversion               = 8000         ;   % total hours in a year, the remaining time is assumed to be down for maitenance  [hr/year]
    life_span_equipment                   = 20           ;   % Life Span of all equuipment besides adsorbent [years]
    life_span_adsorbent                   = 5            ;   % Life Span of adsorbent [years]
    CEPCI                                 = 536.4        ;   % CEPCI of present month year (Jan 2016).
    
    % change this according to the cycle to be simulated
    cycle_time = t_pres + t_ads + t_HR + t_CnCdepres + t_LR ;   % total time required for 1 cycle [s]
    
    EconomicParams    = zeros(6, 1)              ;
    EconomicParams(1) = desired_flow             ;
    EconomicParams(2) = electricity_cost         ;
    EconomicParams(3) = cycle_time               ;
    EconomicParams(4) = hour_to_year_conversion  ;
    EconomicParams(5) = life_span_equipment      ;
    EconomicParams(6) = life_span_adsorbent      ;
    EconomicParams(7) = CEPCI                    ;



    % Call PSA cycle functions
    CoCPressurization_fxn   = @(t, x) FuncCoCPressurization(t, x, Params, IsothermParams)   ;
    Adsorption_fxn          = @(t, x) FuncAdsorption(t, x, Params, IsothermParams)          ;
    HeavyReflux_fxn         = @(t, x) FuncHeavyReflux(t, x, Params, IsothermParams)         ;
    CnCDepressurization_fxn = @(t, x) FuncCnCDepressurization(t, x, Params, IsothermParams) ;
    
    % Retrieve times of PSA steps
    t_CoCPres   = Times(1) ;
    t_ads       = Times(2) ;
    t_HR        = Times(6) ;
    t_CnCDepres = Times(3) ;
    t_LR        = Times(4) ;
    
    % Dimensionless times 
    tau_CoCPres   = t_CoCPres*v_0/L   ;
    tau_ads       = t_ads*v_0/L       ;
    tau_HR        = t_HR*v_0/L        ;
    tau_CnCDepres = t_CnCDepres*v_0/L ;
    tau_LR        = t_LR*v_0/L        ;
    
   
    q                = Isotherm(y_0, P_l, 298.15, IsothermParams) ;
    x0               = zeros(5*N+10,1) ;
    x0(1:N+2)        = P_l/P_0         ;
    x0(N+3)          = y_0             ;
    x0(N+4:2*N+4)    = y_0             ;
    x0(2*N+5:3*N+6)  = q(1)/q_s0       ;
    x0(3*N+7:4*N+8)  = q(2)/q_s0       ;
    x0(4*N+9)        = 1               ;
    x0(4*N+10:5*N+10)= 298.15/T_0      ;

    
%     opts1 = odeset( 'RelTol', 1e-6) ;
%     opts2 = odeset( 'RelTol', 1e-6) ;
%     opts3 = odeset( 'RelTol', 1e-6) ;
%     opts4 = odeset( 'RelTol', 1e-6) ;
%     opts5 = odeset( 'RelTol', 1e-6) ;
    
    opts1 = odeset( 'JPattern', JacPressurization(N), 'RelTol', 1e-6)       ;
    opts2 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    opts3 = odeset( 'JPattern', JacAdsorption(N), 'RelTol', 1e-6)           ;
    opts4 = odeset( 'JPattern', Jac_CnCDepressurization(N), 'RelTol', 1e-6) ;
    opts5 = odeset( 'JPattern', Jac_LightReflux(N), 'RelTol', 1e-6)         ;
%   
%% Begin simulating PSA cycle. This is skipped if the first constraint is violated.
%   Run the simulation until the change in the temperature, gas mole fraction
%   and CO2 molar loading is less than 0.5%. For the molar loading and the
%   mole fraction, the simulation also stops if the absolute change in the
%   state variable is less than 5e-5 and 2.5e-4 respectively. 
    
    % initialize variables to store the desired data to be collected
    a_in  = [] ;
    b_in  = [] ;
    c_in  = [] ;
    d_in  = [] ;
    e_in  = [] ;
    a_fin = [] ;
    b_fin = [] ;
    c_fin = [] ;
    d_fin = [] ;
    e_fin = [] ;
   
    
    for i=1:500
        %disp(['Iteration for CSS condition number: ', num2str(i)])
        
        % Store initial conditions for CoCPressurization step - all iterations
        a_in = [a_in; x0'] ;
        
    %% 1. Simulate CoCPressurization step
        [t1, a] = ode15s(CoCPressurization_fxn, [0 tau_CoCPres], x0, opts1) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(a(:, 1) < a(:, 2))         ;  % P_1  < P_2
        a(idx ,1)       = a(idx, 2)                       ;  % P_1  = P_2
        a(idx, N+3)     = a(idx, N+4)                     ;  % y_1  = y_2
        a(idx, 4*N+9)   = a(idx, 4*N+10)                  ;  % T_1  = T_2
        a(:, 2*N+5)     = a(:, 2*N+6)                     ;  % x1_1 = x1_2
        a(:, 3*N+7)     = a(:, 3*N+8)                     ;  % x2_1 = x2_2
        a(:, 3*N+6)     = a(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        a(:, 4*N+8)     = a(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        a(:, N+3:2*N+4) = max(min(a(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % Store final conditions for CoCPressurization step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t1*L/v_0, a, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t1*L/v_0, a, 'LPEnd') ;
        a_fin = [a_fin; a(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;
        
        % Prepare initial conditions for Adsorption step
        x10         = a(end, :)' ;  % Final state of previous step is the
                                    % initial state for current step
        x10(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x10(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x10(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x10(2*N+4)  = x10(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x10(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x10(5*N+10) = x10(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
        
        % Store initial conditions for Adsorption step - all iterations
        b_in = [b_in; x10'] ;
        %disp('Pressurization done')
        
        % Initial conditions of states at first step of the PSA cycle
        statesIC = a(1, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%       
    %% 2. Simulate Adsorption step
        [t2, b] = ode15s(Adsorption_fxn, [0 tau_ads], x10, opts2) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(b(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        b(idx, N+2)     = b(idx, N+1)                     ;  % P_N+2 = P_N+1
        b(:, 2*N+5)     = b(:, 2*N+6)                     ;  % x1_1 = x1_2
        b(:, 3*N+7)     = b(:, 3*N+8)                     ;  % x2_1 = x2_2
        b(:, 3*N+6)     = b(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        b(:, 4*N+8)     = b(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        b(:, N+3:2*N+4) = max(min(b(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            %b = VelocityCorrection(b, ndot_0, 'HPEnd') ;
            b = velocitycleanup(b)                     ;
        end
        
        % Store final conditions for Adsorption step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t2*L/v_0, b, 'HPEnd') ;
        [totalEnd, CO2End, TEnd]  = StreamCompositionCalculator(t2*L/v_0, b, 'LPEnd') ;
        b_fin = [b_fin; b(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;
        
        % Add and update necessary parameters for Light Reflux step. These
        % are the composition and temperature going out of the adsorption 
        % (light product end), which are those for the inlet of the LR step
        y_LR        = CO2End/totalEnd ;
        T_LR        = TEnd            ;
        ndot_LR     = totalEnd/t_ads  ;
        Params(27)  = y_LR            ;
        Params(28)  = T_LR            ;
        Params(29)  = ndot_LR         ;
        
        % Call the function for Light Reflux step with updated parameters
        LightReflux_fxn = @(t, x) FuncLightReflux(t, x, Params, IsothermParams) ;
        
        % Prepare initial conditions for Heavy Reflux step
        x20         = b(end, :)' ;  % Final state of previous step is the
                                    % initial state for current step
        x20(1)      = P_inlet    ;  % BC z=0 P: P_1   = P_inlet
        x20(1)      = x20(2)     ;   
        x20(N+2)    = 1          ;  % BC z=1 P: P_N+2 = 1
        x20(N+3)    = y_HR       ;  % BC z=0 y: y_1   = y_HR
        x20(2*N+4)  = x20(2*N+3) ;  % BC z=1 y: y_N+2 = y_N+1
        x20(4*N+9)  = T_HR/T_0   ;  % BC z=0 T: T_1   = T_HR/T_0
        x20(5*N+10) = x20(5*N+9) ;  % BC z=1 T: T_N+2 = T_N+1
        
        % Store initial conditions for Heavy Reflux step - all iterations
        c_in = [c_in; x20'] ;
        %disp('Adsorption done')
%       
    %% 3. Simulate Heavy Reflux step
        [t3, c] = ode15s(HeavyReflux_fxn, [0 tau_HR], x20, opts3) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(c(:, N+1) < 1)             ;  % P_N+1 < 1 = P_N+2
        c(idx, N+2)     = c(idx, N+1)                     ;  % P_N+2 = P_N+1
        c(:, 2*N+5)     = c(:, 2*N+6)                     ;  % x1_1 = x1_2
        c(:, 3*N+7)     = c(:, 3*N+8)                     ;  % x2_1 = x2_2
        c(:, 3*N+6)     = c(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        c(:, 4*N+8)     = c(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        c(:, N+3:2*N+4) = max(min(c(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        if Params(end) == 0
            c = VelocityCorrection(c, ndot_HR, 'HPEnd') ;
            %c = velocitycleanup(c)                      ;
        end
        
        % Store final conditions for Heavy Reflux step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t3*L/v_0, c, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t3*L/v_0, c, 'LPEnd') ;
        c_fin = [c_fin; c(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;
        
        % Prepare initial conditions for CoCDepressurization step
        x30         = c(end,:)'   ;   % Final state of previous step is the
                                      % initial state for current step
        x30(1)      = x30(2)      ;   % BC z=0 P: P_1   = P_2
        x30(N+2)    = x30(N+1)    ;   % BC z=1 P: P_N+2 = P_N+1
        x30(N+3)    = x30(N+4)    ;   % BC z=0 y: y_1   = y_2
        x30(2*N+4)  = x30(2*N+3)  ;   % BC z=1 y: y_N+2 = y_N+1
        x30(4*N+9)  = x30(4*N+10) ;   % BC z=0 T: T_1   = T_2
        x30(5*N+10) = x30(5*N+9)  ;   % BC z=1 T: T_N+2 = T_N+1
        
        % Store initial conditions for CoCDepressurization step - all iterations
        d_in = [d_in; x30'] ;
        %disp('High Reflux done')
%       
    %% 4. Simulate CnCDepressurization step
        [t4, d] = ode15s(CnCDepressurization_fxn, [0 tau_CnCDepres], x30, opts4) ;
        
        % correct the output (clean up results from simulation)
        idx             = find(d(:, 2) < d(:, 1))         ;  % P_2  < P_1
        d(idx ,1)       = d(idx, 2)                       ;  % P_1  = P_2
        d(:, 2*N+5)     = d(:, 2*N+6)                     ;  % x1_1 = x1_2
        d(:, 3*N+7)     = d(:, 3*N+8)                     ;  % x2_1 = x2_2
        d(:, 3*N+6)     = d(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        d(:, 4*N+8)     = d(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        d(:, N+3:2*N+4) = max(min(d(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        % Store final donditions for CnCDepressurization step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, CO2Front, ~] = StreamCompositionCalculator(t4*L/v_0, d, 'HPEnd') ;
        [totalEnd, CO2End, ~]     = StreamCompositionCalculator(t4*L/v_0, d, 'LPEnd') ;
        d_fin = [d_fin; d(end, :), CO2Front, totalFront, CO2End, totalEnd]            ;
        
        % Prepare initial donditions for Light Reflux step
        x40         = d(end,:)'    ;  % Final state of previous step is the
                                      % initial state for durrent step
        x40(1)      = P_l/P_0      ;  % Bd z=0 P: P_1   = P_l/P_0
        %x40(N+2)    = 2*P_l/P_0    ;  % Bd z=1 P: P_N+2 = 2*P_l/P_0 % NOTE: be aware of this alpha here, on the other dode is just 2
        x40(N+3)    = x40(N+4)     ;  % Bd z=0 y: y_1   = y_2
        x40(2*N+4)  = y_LR         ;  % Bd z=1 y: y_N+2 = y_LR
        x40(4*N+9)  = x40(4*N+10)  ;  % Bd z=0 T: T_1   = T_2
        x40(5*N+10) = T_LR/T_0     ;  % Bd z=1 T: T_N+2 = T_LR/T_0
        
        % Store initial conditions for Light Reflux step - all iterations
        e_in = [e_in; x40'] ;
        %disp('Depressurization done')
%       
    %% 5. Simulate Light Reflux step
        [t5, e] = ode15s(LightReflux_fxn, [0 tau_LR], x40, opts5) ;
        
        % Correct the output (clean up results from simulation)
        idx             = find(e(:, 2) < e(:, 1))         ;  % P_2  < P_1
        e(idx ,1)       = e(idx, 2)                       ;  % P_1  = P_2
        e(:, 2*N+5)     = e(:, 2*N+6)                     ;  % x1_1 = x1_2
        e(:, 3*N+7)     = e(:, 3*N+8)                     ;  % x2_1 = x2_2
        e(:, 3*N+6)     = e(:, 3*N+5)                     ;  % x1_N+2 = x1_N+1
        e(:, 4*N+8)     = e(:, 4*N+7)                     ;  % x2_N+2 = x2_N+1
        e(:, N+3:2*N+4) = max(min(e(:, N+3:2*N+4), 1), 0) ;  % 0 <= y => 1
        
        e = VelocityCorrection(e, ndot_LR*alpha, 'LPEnd') ;
        %e = velocitycleanup(e)                            ;
        
        % Store final conditions for Light Reflux step - all iterations
        % and the CO2 and total moles at the Front and End of the column
        [totalFront, CO2Front, TFront]  = StreamCompositionCalculator(t5*L/v_0, e, 'HPEnd') ;
        [totalEnd, CO2End, ~]           = StreamCompositionCalculator(t5*L/v_0, e, 'LPEnd') ;
        e_fin = [e_fin; e(end, :), CO2Front, totalFront, CO2End, totalEnd]                  ;
        
        % Calculate necessary parameters for Heavy Reflux step
        y_HR       = CO2Front/totalFront    ;
        T_HR       = TFront                 ;
        ndot_HR    = totalFront.*beta/t_HR  ;
        Params(33) = y_HR                   ;
        Params(34) = T_HR                   ;
        Params(35) = ndot_HR                ;
        
        HeavyReflux_fxn = @(t, x) FuncHeavyReflux(t, x, Params, IsothermParams) ;
        
        % Prepare initial conditions for CoCPressurization step
        x0         = e(end, :)' ;  % Final state of previous step is the
                                   % initial state for current step
        x0(1)      = x0(2)      ;  % BC z=0 P: P_1   = P_2
        x0(N+2)    = x0(N+1)    ;  % BC z=1 P: P_N+2 = P_N+1
        x0(N+3)    = y_0        ;  % BC z=0 y: y_1   = y_0
        x0(2*N+4)  = x0(2*N+3)  ;  % BC z=1 y: y_N+2 = y_N+1
        x0(4*N+9)  = 1          ;  % BC z=0 T: T_1   = 1
        x0(5*N+10) = x0(5*N+9)  ;  % BC z=1 T: T_N+2 = T_N+1
        %disp('Light Reflux done')
        
        % Final conditions of states at lat step of the PSA cycle
        statesFC = e(end, [2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9]) ;
%       
    %% Check CCS condition
        
%         [cyclic_check, cyclic_display] = CCS_Check(a, b, c, d, e, t1, t2, t3, t4, t5) ;
%         
%         if strcmpi(it_disp, 'yes') == 1
%             display([i, cyclic_display]) ;
%         end
%         
%         if cyclic_check == 1
%             break
%         end 
        % CCS condition of states
        CSS_states = norm(statesIC-statesFC) ;
        % Mass balance condition
        [~, ~, massBalance] = ProcessEvaluation(a, b, c, d, e, t1, t2, t3, t4, t5) ;
        
        % Check if CCS has been acheived or not
        if CSS_states <= 0.001 && abs(massBalance-1) <= 0.005
            break
        end
%         % Condition to stop if the mass balance is not satisfied but also is
%         % not changing between ten consecutive iterations
%         mb(i) = massBalance ;
%         if i > 15
%             if CSS_states <= 1e-3 && abs(massBalance-1) > 0.005
%                 stateMB      = mb(end-15+1:end)                      ;
%                 %diffStatesMB = stateMB(end:-1:2)-stateMB(end-1:-1:1) ;
%                 diffStatesMB = stateMB(end-1:-1:1)-stateMB(end)      ;
%                 normStatesMB = norm(diffStatesMB)                    ;
%                 if normStatesMB < 1e-5
%                     break
%                 end 
%             end 
%         end     
%       
    end
    
    %% Process and Economic evaluation
    
    [purity, recovery, MB] = ProcessEvaluation(a, b, c, d, e, t1, t2, t3, t4, t5) ;
    desired_flow = EconomicParams(1) ;
    % cycle_time   = EconomicParams(3) ;
    cycle_time = t_CoCPres + t_ads + t_HR + t_CnCDepres + t_LR       ;
    
    % Calculate the amount of flue gas that is fed during the cycle
    [n_tot_pres, ~, ~] = StreamCompositionCalculator(t1*L/v_0, a, 'HPEnd') ;
    [n_tot_ads, ~, ~]  = StreamCompositionCalculator(t2*L/v_0, b, 'HPEnd') ;
    gas_fed = n_tot_pres + n_tot_ads ;  % mols/m^2
    
    % Calculate the required radius of the column to satisfy the molar flow rate
    radius_inner = sqrt((desired_flow.*cycle_time/gas_fed)/pi()) ;    % m
    r_in = radius_inner;                                          
    
    D = 2*(r_in+0.0508);  %Diameter of column (m with 2 inch thickness)
    Area = pi()*r_in^2;   %(Are of cross-section)m^2
    Inlet_molar_flow = ndot_0*Area;   %(mol/s)
    
    % Calculate the energy required during the Pressurization step
    E_pres  = CompressionEnergy(t1*L/v_0, a, 1e5)   ; % kWh
    
    % Calculate the energy required during the Feed step
    E_feed  = CompressionEnergy(t2*L/v_0, b, 1e5)   ; % kWh
    
    % Calculate the energy required during the Heavy Reflux Step
    E_HR    = CompressionEnergy(t3*L/v_0, c, 1e5)   ; % kWh
    
    % Calculate the energy required during the Counter Current Depressurization step
    E_bldwn = VacuumEnergy(t4*L/v_0, d, 1e5)        ; % kWh
    
    % Calculate the energy required during the Light Reflux Step
    E_evac  = VacuumEnergy(t5*L/v_0, e, 1e5)        ; % kWh
    
    % Calculate the total energy required
    energy_per_cycle = E_pres + E_feed + E_HR + E_bldwn + E_evac ; % [kWh per cycle]
    
    % Calculate the CO2 recovered during the cycle [ton CO_2 per cycle and mol/cycle]
    [~, n_CO2_CnCD, ~]   = StreamCompositionCalculator(t4*L/v_0, d, 'HPEnd') ;
    [~, n_CO2_LR, ~]     = StreamCompositionCalculator(t5*L/v_0, e,  'HPEnd') ;
    CO2_recovered_cycle  = (n_CO2_CnCD+(1-beta)*n_CO2_LR)*r_in^2*pi()*MW_CO2/1e3 ;        
    CO2_recovered_cycle2 = (n_CO2_CnCD+(1-beta)*n_CO2_LR)*r_in^2*pi() ;
    
    %Calculate the productivity of the column and the energy requirments
    mass_adsorbent     = L*pi()*r_in^2*(1-epsilon)*ro_s                    ;
    productivity       = CO2_recovered_cycle2./cycle_time./mass_adsorbent  ;
    energy_requirments = energy_per_cycle./CO2_recovered_cycle             ;
    
    Area_Heat_exchanger = 2 ;    %m^2
    
    
    %Direct Cost of equipments:
    Column_Cost = (exp((0.4148138*D) + (0.0738133*(L/D)) + (0.0231138*P_0/10^5) +  10.807870)) ; %Direct cost of adsorption column per ton CO2 captured
    
    t_depres = t_CnCdepres/3600;  %(Converting second to hour)
    molar_flow = (n_CO2_CnCD*(pi.*r_in^2))/t_depres;      %(mol/h)
    molar_flow_mol_s = (n_CO2_CnCD*(pi.*r_in^2))/t_CnCdepres ;  %(mol/s)
    Volumetric_flow = (molar_flow*8.314*313.15)/(P_0)  ;  %(m^3/h)  
  
    Vacuum_pump_Cost = ((423.9 * (Volumetric_flow).^0.653) + 30000) ;  %Direct cost of vacuum pump per ton CO2 captured
    
    Heat_exchanger_cost = 12003 * (Area_Heat_exchanger)^0.603 * 1.011187^(P_0/(10^5))  ;
    
    Adsorbent_required_kg = (pi()*r_in^2*L) * (1-epsilon) * SimParam(type,1) ;  %kg
    
    Adsorbent_cost = Adsorbent_required_kg * SimParam(type,4)  ;
    
    Fluegas_Inlet_molarflow = ndot_0 * Area  ;   %(mol/s)
    
    CO2capture = recovery * Fluegas_Inlet_molarflow * (t_pres+t_ads) * MW_CO2  ;   %(kg)
    
    TDC = Column_Cost + Vacuum_pump_Cost + Heat_exchanger_cost + Adsorbent_cost ; %Total Direct Cost (TDC)
    
    TDCPC = TDC + 0.15*TDC ;  %Total Direct Cost including Process Contingency
    
    TPC = TDCPC + 0.14*TDCPC + 0.2*TDCPC ; %Total Plant Cost
    
    TCR = TPC + 0.07*TDCPC  ;     %Total Captial requirement
    
 
%   
%% Filter stored data 
    % Prepare the collected data to store it in the output variables only
    % every 5 cycle but guaranteeing that the last cycle (CSS condition)
    % is stored
    
    if mod(i, 5) ==0
        idx_out = [1, 5:5:i]    ;
    else
        idx_out = [1, 5:5:i, i] ;
    end
    
%     a_fin = a_fin(idx_out, :)   ;
%     b_fin = b_fin(idx_out, :)   ;
%     c_fin = c_fin(idx_out, :)   ;
%     d_fin = d_fin(idx_out, :)   ;
%     e_fin = e_fin(idx_out, :)   ;
%     a_in  = a_in(idx_out, :)    ;
%     b_in  = b_in(idx_out, :)    ;
%     c_in  = c_in(idx_out, :)    ;
%     d_in  = d_in(idx_out, :)    ;
%     e_in  = e_in(idx_out, :)    ;
%   
%% Complementary Functions
   
    function [n_tot, n_CO2, Temp] = StreamCompositionCalculator(time, state_vars, ProductEnd)
    %MoleStreamCalculator: Calculate the composition of streams
    %   Calculate the composition (moles of CO2 and total moles), and 
    %   temperature of a stream at any end of the column end. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   ProductEnd: End of the column where the composition will be 
    %               calculated. OPTIONS: HPEnd and LPEnd. HPEnd stands for
    %               heavy product end, which is at the bottom of the column,
    %               and LPEnd stands for light product end, which is at the
    %               top of the column.
    %   
    %   Output:
    %   n_tot     : Average total number of moles at the desired end of the
    %               requested step 
    %   n_CO2     : Average number of moles of CO2 at the desired end of the
    %               requested step 
    %   Temp      : Average temperature at the desired end of the requested
    %               step  
    %   
    %% Check number of inputs
        % If no value is given for the ProductEnd, this is set up by default
        % as HPEnd, which is the heavy product end for any step.
        if nargin < 3
            ProductEnd = 'HPEnd';
        end
    %   
    %%  
        % Differential section length of the column
        dz = L/N ;
        
        % Dimensionalize all variables at the two ends of the columns
        % Collect pressure, temperature and mole fraction at the column end
        % of interest.
        if strcmpi(ProductEnd, 'HPEnd') == 1
        
            P = state_vars(:, 1:2)*P_0     ;
            y = state_vars(:, N+3)         ;
            T = state_vars(:, 4*N+9)*T_0   ;
        
            % Calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO2 + (1-y)*MW_N2).*P(:, 1)/R./T ;
        
            % calculate concentrations [mol/m^3]
            C_tot = P(:, 1)/R./T ;
            C_CO2 = C_tot.*y     ;
        
        elseif strcmpi(ProductEnd, 'LPEnd') == 1
        
            P = state_vars(:, N+1:N+2)*P_0 ;
            y = state_vars(:, 2*N+4)       ;
            T = state_vars(:, 5*N+10)*T_0  ;
        
            % calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO2 + (1-y)*MW_N2).*P(:, 2)/R./T ;
        
            % calculate concentrations [mol/m^3]
            C_tot = P(:, 2)/R./T ;
            C_CO2 = C_tot.*y     ;
        
        else
            error('Please specify in which end of the column the composition will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        % Calculate the pressure gradient at the edges [Pa/m]
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % calculate molar fluxes [mol/m^2/s]
        ndot_tot = abs(v.*C_tot) ;
        ndot_CO2 = abs(v.*C_CO2) ;
        
        % calculate total moles per column area [mol/m^2]
        n_tot = trapz(time, ndot_tot) ;

        n_CO2 = trapz(time, ndot_CO2) ;
        
        % calculate the average temperature of the gas [K]. Only important
        % if the emissions are being used in another step in the cycle (e.g.
        % light reflux, light product pressurization)
        energy_flux_tot = ndot_tot.*T                  ;
        energy_tot      = trapz(time, energy_flux_tot) ;
        Temp            = energy_tot/n_tot         ;
    %   
    end 
    
    function [purity, recovery, mass_balance] = ProcessEvaluation(varargin)
    %ProcessEvaluation: Calculate the process objectives
    %   Calculate the purity and the recovery of the heavy product and the
    %   overall mass balance of the heavy product for the cycle.
    %   
    %   Input:
    %   a-e         : The dimensionless state variables of each step for the
    %                 cycle
    %   t1-t5       : The time steps for each step of the cycle
    %   
    %   Output:
    %   purity      : Purity of the heavy product
    %   recovery    : Recovery of the heavy product
    %   mass_balance: Mass balance of the heavy product (ensure all heavy
    %                 product that enters, leaves. No accumulation)
    %   
    %%  
        step  = cell(nargin/2,1) ;
        tau   = cell(nargin/2,1) ;
        for st = 1:nargin/2
            step{st} = varargin{st}          ;
            tau{st}  = varargin{nargin/2+st} ;
        end
    %   
    %% Calculate total moles going in and out of column [mols/m^2]
        
        [~, n_CO2_CoCPres_HPEnd, ~]                       = StreamCompositionCalculator(tau{1}, step{1}, 'HPEnd') ;
        [~, n_CO2_ads_HPEnd, ~]                           = StreamCompositionCalculator(tau{2}, step{2}, 'HPEnd') ;
        [~, n_CO2_ads_LPEnd, ~]                           = StreamCompositionCalculator(tau{2}, step{2}, 'LPEnd') ;
        [~, n_CO2_HR_LPEnd, ~]                            = StreamCompositionCalculator(tau{3}, step{3}, 'LPEnd') ;
        [~, n_CO2_HR_HPEnd, ~]                            = StreamCompositionCalculator(tau{3}, step{3}, 'HPEnd') ;
        [n_tot_CnCDepres_HPEnd, n_CO2_CnCDepres_HPEnd, ~] = StreamCompositionCalculator(tau{4}, step{4}, 'HPEnd') ;
        [~, n_CO2_LR_LPEnd, ~]                            = StreamCompositionCalculator(tau{5}, step{5}, 'LPEnd') ;
        [n_tot_LR_HPEnd, n_CO2_LR_HPEnd, ~]               = StreamCompositionCalculator(tau{5}, step{5}, 'HPEnd') ;
	%   
    %% Calculate Purity, recovery and mass balance of the column
        
        purity       = (n_CO2_CnCDepres_HPEnd+(1-beta)*n_CO2_LR_HPEnd)/(n_tot_CnCDepres_HPEnd+(1-beta)*n_tot_LR_HPEnd) ;
        recovery     = (n_CO2_CnCDepres_HPEnd+(1-beta)*n_CO2_LR_HPEnd)/(n_CO2_CoCPres_HPEnd+n_CO2_ads_HPEnd)           ;
        
        mass_balance = (n_CO2_CnCDepres_HPEnd+n_CO2_ads_LPEnd+n_CO2_HR_LPEnd+n_CO2_LR_HPEnd)/... 
                       (n_CO2_CoCPres_HPEnd+n_CO2_ads_HPEnd+n_CO2_HR_HPEnd+n_CO2_LR_LPEnd)                             ;
    %   
    end 
    
    function energy = CompressionEnergy(time, state_vars, Patm)
    %CompressionEnergy: Calculate the  compression energy from entrance
    %   Calculate the total energy required by the compressor to increase
    %   the pressure to the desired value. Intended to be used with PSA
    %   dimensionless simulations. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   Patm      : Atmospheric pressure [Pa]. This value can also be changed
    %               if the flue gas is at a higher pressure. This value is
    %               the cutoff for the compressor. Above this value, energy
    %               is required. Below his value, no energy is required
    %   
    %   Output:
    %   energy    : Energy Requirments [kWh]  
    %   
    %% Check number of inputs
    
        % Differential section length of the column
        dz = L/N  ;
        
        % Compresor parameters
        adiabatic_index       = 1.4  ;
        compressor_efficiency = 0.72 ;
        
        % Calculate the pressure gradient at the edges [Pa/m]
        P = state_vars(:, 1:2)*P_0     ;
        y = state_vars(:, N+3)         ;
        T = state_vars(:, 4*N+9)*T_0   ;
        
        % Calculate the density of the gas [kg/m^3]
        ro_g = (y*MW_CO2 + (1-y)*MW_N2).*P(:, 1)/R./T ;
        
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((P(:, 1)/Patm).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                        ;
        
        %Calculate the total energy required by the compressor
        integral_term = abs(v.*P(:, 1).*ratio_term) ;
        
        energy = trapz(time, integral_term).*((adiabatic_index)./(adiabatic_index-1))./compressor_efficiency*pi()*r_in.^2 ;
        
        energy = energy/3.6e6 ;
    end
    
    function energy = VacuumEnergy(time, state_vars, Patm, ProductEnd)
    %VacuumEnergy: Calculate the vacuum energy at both ends of the column
    %   Calculate the total energy required by the vacuum pump to operate 
    %   at the desired pressure. Intended to be used with PSA dimensionless 
    %   simulations. 
    %   
    %   Input:
    %   time      : Dimensional time vector supplied by the ODE solver
    %   state_vars: Dimensionless state variable matrix supplied for the
    %               step from the ODE solver
    %   ProductEnd: End of the column where the composition will be 
    %               calculated. OPTIONS: HPEnd and LPEnd. HPEnd stands for
    %               heavy product end, which is at the bottom of the column,
    %               and LPEnd stands for light product end, which is at the
    %               top of the column.
    %   Patm      : Atmospheric pressure [Pa].This value is the cutoff for 
    %               the vacuum. Above this value, energy is not required.
    %               Below his value, energy is required
    %   
    %   Output:
    %   energy    : Energy Requirments [kWh]
    %   
    %% Check number of inputs
        % If no value is given for the ProductEnd, this is set up by default
        % as HPEnd, which is the heavy product end for any step.
        if nargin < 4
            ProductEnd = 'HPEnd';
        end
    %   
    %%  
        % Differential section length of the column
        dz = L/N ;
        
        % Vacuum parameters
        adiabatic_index   = 1.4  ;
        vacuum_efficiency = 0.72 ;
        
        % Dimensionalize all variables at the two ends of the columns
        % Collect pressure, temperature and mole fraction at the column end
        % of interest.
        if strcmpi(ProductEnd, 'HPEnd') == 1
        
            P = state_vars(:, 1:2)*P_0     ;
            y = state_vars(:, N+3)         ;
            T = state_vars(:, 4*N+9)*T_0   ;
        
            % Calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO2 + (1-y)*MW_N2).*P(:, 1)/R./T ;
            
            P_out = P(:, 1) ;
        
        elseif strcmpi(ProductEnd, 'LPEnd') == 1
        
            P = state_vars(:, N+1:N+2)*P_0 ;
            y = state_vars(:, 2*N+4)       ;
            T = state_vars(:, 5*N+10)*T_0  ;
        
            % calculate the density of the gas [kg/m^3]
            ro_g = (y*MW_CO2 + (1-y)*MW_N2).*P(:, 2)/R./T ;
            
            P_out = P(:, 2) ;
        
        else
            error('Please specify in which end of the column the composition will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        % Calculate the pressure gradient at the edges [Pa/m]
        dPdz = 2*(P(:, 2)-P(:, 1))/dz ;
        
        % calculate superficial velocity using ergun equation [m/s]
        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        % Calculate the compression ratio along with the impact it has on the
        % energy requirments
        ratio_term = ((Patm./P_out).^((adiabatic_index-1)/adiabatic_index)-1) ;
        ratio_term = max(ratio_term, 0)                                       ;
        
        integral_term = abs(v.*P_out.*ratio_term) ;
        
        %Calculate the total energy required by the compressor  
        energy=trapz(time, integral_term).*((adiabatic_index)./(adiabatic_index-1))./vacuum_efficiency*pi()*r_in.^2;
        
        energy = energy/3.6e6 ;
    %   
    end 
    
    function x_new = VelocityCorrection(x, n_hr, CorrectionEnd)
    %% Check number of inputs
        % If no value is given for the CorrectionEnd, this is set up by
        % default as HPEnd, which is the heavy product end for any step.
        if nargin < 3
            CorrectionEnd = 'HPEnd';
        end
    %   
    %%  
        x_new = x   ;
        
        % Differential section length of the column
        dz    = L/N ;
        
        % Dimensionalize all variables at the two ends of the columns
        if strcmpi(CorrectionEnd, 'HPEnd') == 1
        
            T = x(:, 4*N+9)*T_0 ;
            y = x(:, N+3)       ;
            P = x(:, 2)*P_0     ;
        
        elseif strcmpi(CorrectionEnd, 'LPEnd') == 1
        
            T = x(:, 5*N+10)*T_0 ;
            y = x(:, 2*N+4)      ;
            P = x(:, N+1)*P_0    ;
        
        else
            error('Please specify in which end of the column the velocity correction will be calculated. OPTIONS: HPEnd and LPEnd')
        end
        
        MW = MW_N2+(MW_CO2-MW_N2)*y ;
        
        a_1   = 150*mu*(1-epsilon)^2*dz/2/4/r_p^2/epsilon^3/R./T    ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz/2 ;
        a_2   = a_2_1/R./T*n_hr.*MW                                 ;
        
        a_a =  a_1+a_2  ;
        b_b =  P./T/R   ;
        c_c = -n_hr     ;
        
        vh = (-b_b+sqrt(b_b.^2-4.*a_a.*c_c))/2./a_a ;
        
        a_p = a_1.*T*R       ;
        b_p = a_2_1.*MW/R./T ;
        
        % Correction
        if strcmpi(CorrectionEnd, 'HPEnd') == 1
        
            x_new(:, 1)   = ((a_p.*vh+P)./(1-b_p.*vh.*vh))./P_0 ;
        
        elseif strcmpi(CorrectionEnd, 'LPEnd') == 1
        
            x_new(:, N+2) = ((a_p.*vh+P)./(1-b_p.*vh.*vh))./P_0 ;
        
        end
	%   
    end 
    
    function x_new = velocitycleanup(x)

        x_new=x;
        numb1=150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2;
        ro_gent=x(:,2).*P_0/R/T_0;
        numb2_ent=(ro_gent.*(MW_N2+(MW_CO2-MW_N2)*x(:,N+3)).*(1.75*(1-epsilon)/2/r_p/epsilon));


        x_new(:,1)=(numb1*v_0+numb2_ent*v_0*v_0)*L/P_0/2/N + x(:,2);
    end 
%   
%% Jacobian patterns functions
    
    function J_pres      = JacPressurization(N)
    %JacPressurization: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N     : Number of finite volumes used in the column
    %   Output:
    %       J_pres: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4 = ones(N+2, 4)                      ;
        A4 = full(spdiags(B4, -2:1, N+2, N+2)) ;
        
        % One band Jacobian scheme for adsorption/desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_pres = [ A4, A4, A1, A1, A4;... 
                   A4, A4, A1, A1, A4;... 
                   A1, A1, A1, A0, A1;... 
                   A1, A1, A0, A1, A1;... 
                   A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_pres(1, :) = 0 ;
        J_pres(1, 1) = 1 ;
        
        % Pressure Outlet
        J_pres(N+2, :) = J_pres(N+1, :) ;
        J_pres(:, N+2) = 0              ;
        
        % Mole Fraction Inlet
        J_pres(N+3, :) = 0 ;
        J_pres(:, N+3) = 0 ;
        
        % Mole Fraction Outlet
        J_pres(2*N+4, :) = J_pres(2*N+3, :) ;
        J_pres(:, 2*N+4) = 0                ;
        
        % Temperature Inlet
        J_pres(4*N+9, :) = 0 ;
        J_pres(:, 4*N+9) = 0 ;
        
        % Temperature Outlet
        J_pres(5*N+10, :) = J_pres(5*N+9, :) ;
        J_pres(:, 5*N+10) = 0                ;
        
        J_pres=sparse(J_pres);
	%   
    end 
    
    function J_ads       = JacAdsorption(N)
    %JacAdsorption: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N    : Number of finite volumes used in the column
    %   Output:
    %       J_ads: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4 = ones(N+2, 4)                      ;
        A4 = full(spdiags(B4, -2:1, N+2, N+2)) ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_ads = [ A4, A4, A1, A1, A4;... 
                  A4, A4, A1, A1, A4;... 
                  A1, A1, A1, A0, A1;... 
                  A1, A1, A0, A1, A1;... 
                  A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_ads(1, :) = 0 ;
        J_ads(:, 1) = 0 ;
        
        % Pressure Outlet
        J_ads(N+2, :) = 0 ;
        J_ads(:, N+2) = 0 ;
        
        % Mole Fraction Inlet
        J_ads(N+3, :) = 0 ;
        J_ads(:, N+3) = 0 ;
        
        % Mole Fraction Outlet
        J_ads(2*N+4, :) = J_ads(2*N+3, :) ;
        J_ads(:, 2*N+4) = 0               ;
        
        % Temperature Inlet
        J_ads(4*N+9, :) = 0 ;
        J_ads(:, 4*N+9) = 0 ;
        
        % Temperature Outlet
        J_ads(5*N+10, :) = J_ads(5*N+9, :) ;
        J_ads(:, 5*N+10) = 0               ;
        
        J_ads = sparse(J_ads) ;
	%   
    end 
    
    function J_CnCdepres = Jac_CnCDepressurization(N)
    %Jac_CnCDepressurization: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N          : Number of finite volumes used in the column
    %   Output:
    %       J_CnCdepres: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4         = ones(N+2, 4)                      ;
        A4         = full(spdiags(B4, -1:2, N+2, N+2)) ;
        A4(1, :)   = A4(2, :)                          ;
        A4(N+2, :) = A4(N+1, :)                        ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        A1(1, 2)     = 1                              ;
        A1(N+2, N+1) = 1                              ;
        
       % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_CnCdepres = [ A4, A4, A1, A1, A4;... 
                        A4, A4, A1, A1, A4;... 
                        A1, A1, A1, A0, A1;... 
                        A1, A1, A0, A1, A1;... 
                        A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_CnCdepres(1, :) = 0 ;
        J_CnCdepres(1, 1) = 1 ;
        
        % Pressure Outlet
        J_CnCdepres(N+2, :) = J_CnCdepres(N+1, :) ;
        
        % Mole Fraction Inlet
        J_CnCdepres(N+3, :) = J_CnCdepres(N+4, :) ;
        
        % Mole Fraction Outlet
        J_CnCdepres(2*N+4) = J_CnCdepres(2*N+3) ;
        
        % Molar Loading
        J_CnCdepres(2*N+5, :)       = 0 ;
        J_CnCdepres(3*N+6:3*N+7, :) = 0 ;
        J_CnCdepres(4*N+8, :)       = 0 ;
        
        % Temperature Inlet
        J_CnCdepres(4*N+9, :) = J_CnCdepres(4*N+10, :) ;
        
        %Temperature Outlet
        J_CnCdepres(5*N+10, :) = J_CnCdepres(5*N+9) ;
        
        J_CnCdepres = sparse(J_CnCdepres) ;
	%   
    end 
    
    function J_LR        = Jac_LightReflux(N)
    %Jac_LightReflux: Calculates a Jacobian pattern for the step
    %   This function calculates a jacobian pattern for the pressurization
    %   step. This is to be used with the ode solver to decrease the 
    %   computational time
    %   
    %   Input:
    %       N   : Number of finite volumes used in the column
    %   Output:
    %       J_LR: The sparse Jacobian pattern scheme
    %   
    %% Create individual segments
        % Four band Jacobian scheme for advection terms
        B4         = ones(N+2, 4)                      ;
        A4         = full(spdiags(B4, -1:2, N+2, N+2)) ;
        A4(1, :)   = A4(2, :)                          ;
        A4(N+2, :) = A4(N+1, :)                        ;
        
        % One band Jacobian scheme for adsorption/ desoprtion term
        B1           = ones(N+2, 1)                   ;
        A1           = full(spdiags(B1, 0, N+2, N+2)) ;
        A1(1, 1)     = 0                              ;
        A1(N+2, N+2) = 0                              ;
        A1(1, 2)     = 1                              ;
        A1(N+2, N+1) = 1                              ;
        
        % Zero band Jacobian Term
        A0 = zeros(N+2) ;
	%   
    %% Create Overall Jacobian based on individual segments
        J_LR = [ A4, A1, A1, A1, A4;... 
                 A4, A4, A1, A1, A4;... 
                 A1, A1, A1, A0, A1;... 
                 A1, A1, A0, A1, A1;... 
                 A4, A1, A1, A1, A4 ]  ;
	%   
    %% Modify based on boundary conditions
        % Pressure Inlet
        J_LR(1, :) = 0 ;
        J_LR(1, 1) = 1 ;
        
        % Pressure Outlet
        J_LR(N+2, :) = J_LR(N+1, :) ;
        
        % Mole Fraction Inlet
        J_LR(N+3, :) = J_LR(N+4, :) ;
        
        % Mole Fraction Outlet
        J_LR(2*N+4) = J_LR(2*N+3) ;
        
        % Molar Loading
        J_LR(2*N+5, :)       = 0 ;
        J_LR(3*N+6:3*N+7, :) = 0 ;
        J_LR(4*N+8, :)       = 0 ;
        
        % Temperature Inlet
        J_LR(4*N+9, :) = J_LR(4*N+10, :) ;
        
        % Temperature Outlet
        J_LR(5*N+10, :) = J_LR(5*N+9) ;
        
        J_LR = sparse(J_LR) ;
	%   
    end 
%   
%% Press and Depress Termination functions
    
    function [value, isterminal, direction] = PressurizationStop(~, y)
    %PressurizationStop: Function used with ode15s to stop the step
    %   Termination function for the depressurization step to stop the
    %   solution once the desired pressure is reached
    %
    %%  
        PressureOutlet = y(N+2)              ;
        value          = 0.99-PressureOutlet ;
        isterminal     = 1                   ;
        direction      = 0                   ;
	%   
    end 
    
    function [value, isterminal, direction] = CnCDepressurizationStop(~, y)
    %CnCDepressurizationStop: Function used with ode15s to stop the step
    %   Termination function for the CnC depressurization step to stop the
    %   solution once the desired pressure is reached
    %
    %%  
        PressureOutlet =  real(y(N+2)) ;
        value          =  PressureOutlet-1.9*P_l/P_0 ;
        isterminal     =  1                          ;
        direction      = -1                          ;
    %   
    end 
%   
end 
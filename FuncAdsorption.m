function derivatives = FuncAdsorption(~, state_vars, params, isotherm_params)

% Retrieve process parameters
    N				=	params(1)	;
    deltaU_1	    =	params(2)	;
    deltaU_2	    =	params(3)	;
    ro_s			=	params(4)	;
    T_0				=	params(5)	;
    epsilon			=	params(6)	;
    r_p				=	params(7)	;
    mu				=	params(8)	;
    R				=	params(9)	;
    v_0				=	params(10)	;
    q_s0			=	params(11)	;
    C_pg			=	params(12)	;
    C_pa			=	params(13)	;
    C_ps			=	params(14)	;
    D_m				=	params(15)	;
    K_z				=	params(16)	;
    P_0				=	params(17)	;
    L				=	params(18)	;
    MW_CO2			=	params(19)	;
    MW_N2			=	params(20)	;
    k_1_LDF		    =	params(21)	;
    k_2_LDF		    =	params(22)	;
    y_0				=	params(23)	;
    P_inlet			=	params(26)	;
    ndot_0          =   P_0/R/T_0*v_0   ;
%   
%Initialize state variables
    P  = zeros(N+2, 1) ;
    y  = zeros(N+2, 1) ;
    x1 = zeros(N+2, 1) ;
    x2 = zeros(N+2, 1) ;
    T  = zeros(N+2, 1) ;
    
    P(1:N+2)  = state_vars(1:N+2)               ;
    y(1:N+2)  = max(state_vars(N+3:2*N+4), 0)   ;
    x1(1:N+2) = max(state_vars(2*N+5:3*N+6), 0) ;
    x2(1:N+2) = state_vars(3*N+7:4*N+8)         ;
    T(1:N+2)  = state_vars(4*N+9:5*N+10)        ;
%   
% Initialize all variables used in the function
    % Temporal derivatives
    derivatives = zeros(5*N+10, 1) ;
    dPdt        = zeros(N+2, 1)    ;
    dPdt1       = zeros(N+2, 1)    ;
    dPdt2       = zeros(N+2, 1)    ;
    dPdt3       = zeros(N+2, 1)    ;
    dydt        = zeros(N+2, 1)    ;
    dydt1       = zeros(N+2, 1)    ;
    dydt2       = zeros(N+2, 1)    ;
    dydt3       = zeros(N+2, 1)    ;
    dx1dt       = zeros(N+2, 1)    ;
    dx2dt       = zeros(N+2, 1)    ;
    dTdt        = zeros(N+2, 1)    ;
    dTdt1       = zeros(N+2, 1)    ;
    dTdt2       = zeros(N+2, 1)    ;
    dTdt3       = zeros(N+2, 1)    ;
    % Spatial derivatives
    dpdz        = zeros(N+2, 1)    ;
    dpdzh       = zeros(N+1, 1)    ;
    dydz        = zeros(N+2, 1)    ;
    d2ydz2      = zeros(N+2, 1)    ;
    dTdz        = zeros(N+2, 1)    ;
    d2Tdz2      = zeros(N+2, 1)    ;

   
% Calculate all parameters used
    dz   = 1/N                                ;
    D_l  = 0.7*D_m + v_0*r_p                  ;
    Pe   = v_0*L/D_l                          ;
    phi  = R*T_0*q_s0*(1-epsilon)/epsilon/P_0 ;
    ro_g = P(1:N+2).*P_0/R./T(1:N+2)/T_0      ;
%   

    
    y(1) = y_0     ;
    T(1) = T_0/T_0 ;
    
    if params(end) == 1
        % Inlet pressure is specified as a constant
        P(1) = P_inlet ;
    elseif params(end) == 0
        % Inlet velocity is specified as a constant
        vh = zeros(N+1, 1) ;
        
        MW = MW_N2+(MW_CO2-MW_N2)*y(1) ;
    
        a_1   = 150*mu*(1-epsilon)^2*dz*L/2/4/r_p^2/epsilon^3/T(1)/T_0/R ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz*L/2    ;
        a_2   = a_2_1/R/T(1)/T_0*ndot_0*MW                               ;
    
        a =  a_1+a_2             ;
        b =  P(2)/T(1)*P_0/R/T_0 ;
        c = -ndot_0              ;
    
        vh(1) = (-b+sqrt(b^2-4*a*c))/2/a/v_0 ;
    
        a_p = a_1*T(1)*T_0*R      ;
        b_p = a_2_1*MW/R/T(1)/T_0 ;
    
        P(1)  = ((a_p*vh(1)*v_0+P(2)*P_0)./(1-b_p*vh(1)*v_0*vh(1)*v_0))/P_0 ;
		
        % to reproduce results of chapter 5 of Karson's thesis, the code
        % below has to be added, this is the issue i found in the Karson codes 
        viscous_term=150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2;
        vh(1) = 1 ;
        P(1)=((viscous_term.*vh(1)*v_0*dz/2*L./P_0)+P(2))./(1-(dz/2*L/R/T(1)/T_0)*MW*(1.75*(1-epsilon)/2/r_p/epsilon)*vh(1).^2.*v_0.^2);
        
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end

    y(N+2) = y(N+1)     ;
    T(N+2) = T(N+1)     ;
    if P(N+1) >= 1
        P(N+2) = 1      ;
    else
        P(N+2) = P(N+1) ;
    end
%   
% Spatial Derivative Calculations

    
    Ph          = WENO(P, 'upwind')      ;
    
    dpdz(2:N+1) = (Ph(2:N+1)-Ph(1:N))/dz ;
    dpdzh(2:N)  = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(1)    =  2*(P(2)-P(1))/dz      ;
    dpdzh(N+1)  =  2*(P(N+2)-P(N+1))/dz  ;
%   
%  
%   Mole Fraction: at the center of the volumes
    
    yh          = WENO(y, 'upwind')      ;
    
    dydz(2:N+1) = (yh(2:N+1)-yh(1:N))/dz ;
%   
%  
%   Temperature: at the center of the volumes
    
    Th          = WENO(T, 'upwind')      ;
    
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;
%   

%  
%   Mole Fraction
    d2ydz2(3:N) = (y(4:N+1)+y(2:N-1)-2*y(3:N))/dz/dz ;
    d2ydz2(2)   = (y(3)-y(2))/dz/dz                  ;
    d2ydz2(N+1) = (y(N)-y(N+1))/dz/dz                ;
%   
% 
%   Temperature
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;
%   
% Velocity Calculations
%   Calculates the interstitial velocity of the gas at the walls of the
%   volumes based off of the pressure gradients. NOTE: bear in mind that
%   the velocity in the Ergun's equation is the superficial velocity, and
%   not the interstitial, that's why in the denominator in the viscous term
%   it is found epsilon^2 and not epsilon^3, the same for kinetic terms,
%   where it is found epsilon and not epsilon^3. Superficial velocity is 
%   equal to interstitial velocity times the void fraction

    ro_gh          = (P_0/R/T_0)*Ph(1:N+1)./Th(1:N+1)               ;
    
    viscous_term   = 150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2         ;
    kinetic_term_h = (ro_gh.*(MW_N2+(MW_CO2-MW_N2).*yh)).*(1.75*(1-epsilon)...
                     /2/r_p/epsilon)                                         ;
                    
    % Velocities at walls of volumes 
    vh = -sign(dpdzh).*(-viscous_term+(abs(viscous_term^2+4*kinetic_term_h... 
               .*abs(dpdzh)*P_0/L)).^(.5))/2./kinetic_term_h/v_0             ;
%   
% Temporal Derivatives

%   1.1) Calculate equilibrium molar loading
    q   = Isotherm(y, P*P_0, T*T_0, isotherm_params) ;
    q_1 = q(:, 1)*ro_s                               ;
    q_2 = q(:, 2)*ro_s                               ;
%   
    k_1 = k_1_LDF*L/v_0 ;
    k_2 = k_2_LDF*L/v_0 ;

    dx1dt(2:N+1) = k_1*(q_1(2:N+1)/q_s0 - x1(2:N+1)) ;
    dx2dt(2:N+1) = k_2*(q_2(2:N+1)/q_s0 - x2(2:N+1)) ;

    % [J/m^3/K]
    sink_term = ((1-epsilon)*(ro_s*C_ps+q_s0*C_pa)+(epsilon.*ro_g(2:N+1).*C_pg)) ;

 

    transfer_term = K_z./v_0./L                             ;
    dTdt1(2:N+1)  = transfer_term.*d2Tdz2(2:N+1)./sink_term ;
%   

    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1) ;
    Pv           = Ph(1:N+1).*vh(1:N+1)            ;
    dTdt2(2:N+1) = -epsilon.*C_pg.*P_0./R./T_0.*((Pv(2:N+1)-Pv(1:N))- ... 
                    T(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz./sink_term      ;

    generation_term_1 = (1-epsilon).*q_s0.*(-(deltaU_1-R*T(2:N+1)*T_0))./T_0 ;
    generation_term_2 = (1-epsilon).*q_s0.*(-(deltaU_2-R*T(2:N+1)*T_0))./T_0 ;
    
    dTdt3(2:N+1)      = (generation_term_1.*dx1dt(2:N+1)+...
                         generation_term_2.*dx2dt(2:N+1))./sink_term  ;
%   
%   2.4) Total sum of all temperature derivatives 
    dTdt(2:N+1) = dTdt1(2:N+1) + dTdt2(2:N+1) + dTdt3(2:N+1) ;

    dPdt1(2:N+1) = -T(2:N+1).*(PvT(2:N+1)-PvT(1:N))./dz  ;
%   

    dPdt2(2:N+1) = -phi*T(2:N+1).*(dx1dt(2:N+1)+dx2dt(2:N+1)) ;

    dPdt3(2:N+1) = P(2:N+1).*dTdt(2:N+1)./T(2:N+1)            ;
%   

    dPdt(2:N+1) = dPdt1(2:N+1) + dPdt2(2:N+1) + dPdt3(2:N+1)  ;
%   

%   
 
    dydt1(2:N+1) = (1/Pe)*(d2ydz2(2:N+1)+(dydz(2:N+1).*dpdz(2:N+1)./P(2:N+1))... 
                  -(dydz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
%   

    ypvt         = yh(1:N+1).*Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;
    dydt2(2:N+1) = -(T(2:N+1)./P(2:N+1)).*((ypvt(2:N+1)-ypvt(1:N))... 
                   -y(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz             ;
 
    dydt3(2:N+1) = (phi*T(2:N+1)./P(2:N+1)).*((y(2:N+1)-1).*dx1dt(2:N+1)... 
                  + y(2:N+1).*dx2dt(2:N+1))                                ;

    dydt(2:N+1) = dydt1(2:N+1) + dydt2(2:N+1) + dydt3(2:N+1) ;
%   
%  Boundary Derivatives
    dPdt(1)    = 0         ;
    dPdt(N+2)  = 0         ;
    dydt(1)    = 0         ;
    dydt(N+2)  = dydt(N+1) ;
    dx1dt(1)   = 0         ;
    dx2dt(1)   = 0         ;
    dx1dt(N+2) = 0         ;
    dx2dt(N+2) = 0         ;
    dTdt(1)    = 0         ;
    dTdt(N+2)  = dTdt(N+1) ;
%   
%  Export derivatives to output
    derivatives(1:N+2)        = dPdt(1:N+2)  ;
    derivatives(N+3:2*N+4)    = dydt(1:N+2)  ;
    derivatives(2*N+5:3*N+6)  = dx1dt(1:N+2) ;
    derivatives(3*N+7:4*N+8)  = dx2dt(1:N+2) ;
    derivatives(4*N+9:5*N+10) = dTdt(1:N+2)  ;
%   
end 
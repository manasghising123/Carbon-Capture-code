function flux_w = WENO(flux_c, FlowDir)
    oo     = 10^-10              ;
    [N, m] = size(flux_c)        ;
    N      = N-2                 ;
    flux_w = zeros(N+1, m)       ;
    alpha0 = zeros(size(flux_c)) ;
    alpha1 = zeros(size(flux_c)) ;
    
    % Fluxes at boundaries of the domain
    flux_w(1, :)   = flux_c(1, :)   ;
    flux_w(N+1, :) = flux_c(N+2, :) ;
    
    if strcmpi(FlowDir, 'upwind') == 1
        
        alpha0(2:N, :) =(2/3)./((flux_c(3:N+1, :)-flux_c(2:N, :)+oo).^4) ;
        alpha1(3:N, :) =(1/3)./((flux_c(3:N, :)-flux_c(2:N-1, :)+oo).^4) ;
        alpha1(2, :)   =(1/3)./((2*(flux_c(2, :)-flux_c(1, :))+oo).^4)   ;
        
        flux_w(3:N, :) = (alpha0(3:N, :)./(alpha0(3:N, :)+alpha1(3:N, :)))...      
                       .*((flux_c(3:N, :)+flux_c(4:N+1, :))./2)+(alpha1(3:N, :)... 
                       ./(alpha0(3:N, :)+alpha1(3:N, :))).*(1.5*flux_c(3:N, :)...  
                        -.5*flux_c(2:N-1, :))                                     ;
        
        flux_w(2, :)   = (alpha0(2, :)./(alpha0(2, :)+alpha1(2, :)))...    
                       .*((flux_c(2, :)+flux_c(3, :))./2)+(alpha1(2, :)... 
                       ./(alpha0(2, :)+alpha1(2, :))).*(2*flux_c(2, :)...  
                         -flux_c(1, :))                                   ;
        
    elseif strcmpi(FlowDir, 'downwind') == 1
        
        alpha0(2:N, :)   = (2/3)./((flux_c(2:N, :)-flux_c(3:N+1, :)+oo).^4)   ;
        alpha1(2:N-1, :) = (1/3)./((flux_c(3:N, :)-flux_c(4:N+1, :)+oo).^4)   ;
        alpha1(N, :)     = (1/3)./((2*(flux_c(N+1, :)-flux_c(N+2, :))+oo).^4) ;
        
        flux_w(2:N-1, :) = (alpha0(2:N-1, :)./(alpha0(2:N-1, :)+alpha1(2:N-1, :)))...   
                         .*((flux_c(2:N-1, :)+flux_c(3:N, :))./2)+(alpha1(2:N-1, :)...  
                         ./(alpha0(2:N-1, :)+alpha1(2:N-1, :))).*(1.5*flux_c(3:N, :)... 
                          -.5*flux_c(4:N+1, :))                                        ;
                      
        flux_w(N, :)     = (alpha0(N, :)./(alpha0(N, :)+alpha1(N, :)))...      
                         .*((flux_c(N, :)+flux_c(N+1, :))./2)+(alpha1(N, :)... 
                         ./(alpha0(N, :)+alpha1(N, :))).*(2*flux_c(N+1, :)...  
                           -flux_c(N+2, :))                                   ;
    else
        error('Please specify the direction of flow. OPTIONS: upwind and downwind')
    end 
end 

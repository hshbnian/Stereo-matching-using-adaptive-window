function [sigma_squared_d , r_v , c_v, expanding_priority_flag]=...
    sigma_squared_calculator(r_v ,c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, expand_side_flage, expanding_priority_flag)
% HS 11July2019
% calculating sigma squared _d, represents the uncertainty of the estimation
% for the given window w

switch expand_side_flage
    case 'left'
        r_v(1) = r_v(1) - 1; % expanding in the top direction
        if r_v(1) == 0 
            r_v(1) = 1;
            expanding_priority_flag = 0; % otherwise we would be stuck in this location 
        else
            expanding_priority_flag = 1;
        end 

        
    case 'right'
        r_v(2) = r_v(2) + 1;
        if r_v(2) > m
            r_v(2) = m;
            expanding_priority_flag = 0;
        else
            expanding_priority_flag = 1;
        end 

        
    case 'top'
        c_v(1) = c_v(1) - 1;
        if c_v (1) == 0 
            c_v(1) = 1;
            expanding_priority_flag = 0; 
        else
            expanding_priority_flag = 1;
        end 

        
    case 'bottom'
        c_v(2) = c_v(2) + 1;
        if c_v(2) > n
            c_v(2) = n;
            expanding_priority_flag = 0; 
        else
            expanding_priority_flag = 1;
        end 

        
end   
   
if  expanding_priority_flag == 1  
    for eta = r_v(1) : r_v(2)
        for epsilon =  c_v(1) : c_v(2)
            d2 = sigma_two + (alfa_f * alfa_d * sqrt (eta ^2 + epsilon^2  ));
            %new_loc = epsilon + d_rc;
            new_loc = epsilon - d_rc;
            if new_loc <=0
                 new_loc =1;
            end
            sigma_squared_d_denominator= sigma_squared_d_denominator + (Deriv (eta , new_loc) ^ 2) /  d2 ; 

        end
    end
end
sigma_squared_d = 1 / double(sigma_squared_d_denominator);

% if (window_size_limit_r == 15) || (window_size_limit_c == 15 )
%      expanding_priority_flag = 0;
% end
end
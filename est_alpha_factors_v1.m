function [alpha_d, alpha_f] = est_alpha_factors_v1(disparity_map, d_rc, Deriv, loc_r, loc_c, r_v, c_v)
%Calculating alfa_d , representing the disparity fluctuations
%Calculating alfa_f, representing the intensity fluctuation
%MB July 23, 2019

r_top = r_v(1);
r_bottom = r_v(2);
c_top = c_v(1);
c_bottom = c_v(2);
%
% N_w is a number of samples within the window
N_w = (r_bottom - r_top + 1) * (c_bottom - c_top + 1);
%
alpha_d = 0;
alpha_f = 0;
for eta = r_top : r_bottom
    for epsilon =  c_top :  c_bottom
        %alpha_d = alpha_d +((disparity_map( eta , epsilon) -  d_rc )^ 2) / sqrt (eta^2 + epsilon ^2);
        %new_loc = epsilon + d_rc;
        alpha_d = alpha_d +((disparity_map( eta , epsilon) -  d_rc )^ 2) / sqrt ((eta - loc_r)^2 + (epsilon - loc_c)^2);
        new_loc = epsilon - d_rc;
        if (new_loc) < 1
            new_loc = 1;
        end
        alpha_f = alpha_f + (Deriv (eta, new_loc) ^ 2);
    end
end
alpha_d = double(alpha_d) / double(N_w);
alpha_f = double(alpha_f)/ double(N_w);
end
function [delta_d, sigma_squared_d]=...
    disparity_update_v4(left_image_mat, right_image_mat, loc_r, loc_c, disparity_map, window_r, window_c)
%HS 03July2019
%Calculating Sigma_Squared, which is an uncertainty value based on the
%local variation of the intensity and disparity variation within a window
%Takeo Kanade and et al's proposed method 1994
%ref: A stereo matching algorithm with an adaptive window: teory and experiment
%
%MB, July 17, 2019: Delta d to be calculated after the correct (maximum)
%window size is determined
%MB, July 23, 2019: v4 version - Updated formula based on f1(x, y) = f2(x - d, y)
%Note: paper uses f1(x, y) = f2(x + d, y)

%Parameters
sigma = sqrt(0.038); %1; 
max_window_size = 7;

% Parameters Initialization 
d_rc = disparity_map(loc_r , loc_c); % constant for a given (r,c))
sigma_two= 2 * (sigma ^ 2);
sigma_squared_d_denominator = 0;
%
[m, n] = size(rgb2gray(left_image_mat));
prohibition_flag_top = 0;
prohibition_flag_bottom = 0;
prohibition_flag_left = 0; 
prohibition_flag_right = 0;
expanding_priority_flag_top= 1;
expanding_priority_flag_bottom= 1;
expanding_priority_flag_left= 1; 
expanding_priority_flag_right = 1;
expanding_priority_flag_center= 1;
prohibition_flag_v =0;
% window_size_limit_r=0;
% window_size_limit_c=0;

%first order derivative of the right image in x direction
Deriv = diff(right_image_mat);
Deriv(end+1, :)=Deriv(end, :);% since the number of row of the derivative matrix is less in one row, I duplicated the last row

% window adjustment
r_top = loc_r - window_r; %loc_r - window_r + 1;
if r_top < 1
   r_top = loc_r ;
end   

r_bottom = loc_r + window_r; %loc_r + window_r - 1;
if r_bottom > m
    r_bottom = m;
end

c_top = loc_c - window_c;%loc_c - window_c + 1;
if c_top < 1
   c_top = loc_c;
end

c_bottom = loc_c + window_c; %loc_c + window_c - 1;
if c_bottom > n
    c_bottom = n;
end

% Coordinates of row and columns 
r_v = [ r_top , r_bottom ] ;
c_v = [ c_top , c_bottom ];

%Calculating alfa_d , representing the disparity fluctuations
%Calculating alfa_f, representing the intensity fluctuation
% [alpha_d, alpha_f] = est_alpha_factors(disparity_map, d_rc, Deriv, r_v, c_v);
 [alpha_d, alpha_f] = est_alpha_factors_v1(disparity_map, d_rc, Deriv, loc_r, loc_c, r_v, c_v); %MB July 23, 2019

% calculating sigma squared for the current window
%sigma_squared_d_current = sigma_squared_calculator (r_v ,c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'current', expanding_priority_flag_center);
%sigma_squared_d_current = sigma_squared_calculator_v1(loc_r, loc_c, r_v ,c_v, alpha_f, alpha_d, ...
%    sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'current', expanding_priority_flag_center); %MB July 23, 2019
[sigma_squared_d_current, ~, ~, window_size_v] = sigma_squared_calculator_v2(loc_r, loc_c, ...
    r_v ,c_v, disparity_map, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, ...
    m, n, 'current', expanding_priority_flag_center); %MB July 23, 2019

%while all prohibitation flags are zeros(if there is still a way to expand )
% all(prohibition_flag_v) means ==> it shows one if all directions got prohibited

stop_flag = 0;
%r_prior = r_v;
%c_prior = c_v;
% figure; imshow(disparity_map,[0,55]) ;colorbar; colormap jet;
while ~( stop_flag || all(prohibition_flag_v) )
 if   ( all(window_size_v >= max_window_size) )
     stop_flag = 1;
 end
 % we need to monitor expantion of the window based on its prior value in each side of the window 
%window_top_range =  (c_prior(1)- c_v(1)) + abs(c_prior(1) - c_prior(2));
%window_bottom_range = c_v(2)- c_prior(2)+ abs(c_prior(1) - c_prior(2)); 
%window_left_range =  r_prior(1)- r_v(1) + abs(c_prior(1) - c_prior(2));
%window_right_range = r_v(2) - r_prior(2) + abs(c_prior(1) - c_prior(2));
%window_size_limit_vec = [window_top_range, window_bottom_range, window_left_range, window_right_range ];


%fprintf(' window size==> Row:%d:%d Column:%d:%d ==> current uncertainty:%.4f \n ', r_v(1), r_v(2), c_v(1), c_v(2), sigma_squared_d_current );
%fprintf('current uncertainty:%.4f \n',  sigma_squared_d_current);

%calculating sigma squared for the expanded window in 4 directions
    if (prohibition_flag_top ~= 1 && expanding_priority_flag_top ==1)
        %[sigma_squared_d_top ,r_v_top , c_v_top, expanding_priority_flag_top] = sigma_squared_calculator (r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'top', expanding_priority_flag_top);
        %[sigma_squared_d_top ,r_v_top , c_v_top, expanding_priority_flag_top] = ...
        %    sigma_squared_calculator_v1(loc_r, loc_c, r_v, c_v, alpha_f, alpha_d, sigma_two, ...
        %    d_rc, Deriv, sigma_squared_d_denominator, m, n, 'top', expanding_priority_flag_top); %MB July 23, 2019
        [sigma_squared_d_top ,r_v_top , c_v_top, expanding_priority_flag_top, window_size_v] = ...
            sigma_squared_calculator_v2(loc_r, loc_c, r_v, c_v, disparity_map, sigma_two, ...
            d_rc, Deriv, sigma_squared_d_denominator, m, n, 'top', expanding_priority_flag_top); %MB July 23, 2019
        if sigma_squared_d_top > sigma_squared_d_current
             prohibition_flag_top = 1;
        end
        if ( window_size_v(1) >=  max_window_size )
            expanding_priority_flag_top = 0;
        end
    end

    if (prohibition_flag_bottom ~= 1 && expanding_priority_flag_bottom ==1)
        %[sigma_squared_d_bottom , r_v_bottom , c_v_bottom  , expanding_priority_flag_bottom] = sigma_squared_calculator (r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'bottom',expanding_priority_flag_bottom);
        %[sigma_squared_d_bottom, r_v_bottom, c_v_bottom, expanding_priority_flag_bottom] = ...
        %    sigma_squared_calculator_v1(loc_r, loc_c, r_v, c_v, alpha_f, alpha_d, sigma_two, ...
        %    d_rc, Deriv, sigma_squared_d_denominator, m, n, 'bottom',expanding_priority_flag_bottom); %MB July 23, 2019
        [sigma_squared_d_bottom, r_v_bottom, c_v_bottom, expanding_priority_flag_bottom, window_size_r] = ...
            sigma_squared_calculator_v2(loc_r, loc_c, r_v, c_v, disparity_map, sigma_two, ...
            d_rc, Deriv, sigma_squared_d_denominator, m, n, 'bottom',expanding_priority_flag_bottom); %MB July 23, 2019
        if sigma_squared_d_bottom > sigma_squared_d_current
            prohibition_flag_bottom = 1;
        end
        if (  window_size_v(2) >=  max_window_size )
            expanding_priority_flag_bottom = 0;
        end
        
    end

    if (prohibition_flag_left ~= 1 && expanding_priority_flag_left == 1 )
        %[sigma_squared_d_left, r_v_left , c_v_left , expanding_priority_flag_left]= sigma_squared_calculator (r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'left',  expanding_priority_flag_left);
        %[sigma_squared_d_left, r_v_left, c_v_left, expanding_priority_flag_left]= ...
        %    sigma_squared_calculator_v1(loc_r, loc_c, r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, ...
        %    sigma_squared_d_denominator, m, n, 'left',  expanding_priority_flag_left); %MB July 23, 2019
        [sigma_squared_d_left, r_v_left, c_v_left, expanding_priority_flag_left, window_size_v]= ...
            sigma_squared_calculator_v2(loc_r, loc_c, r_v, c_v, disparity_map, sigma_two, d_rc, Deriv, ...
            sigma_squared_d_denominator, m, n, 'left',  expanding_priority_flag_left); %MB July 23, 2019
        if sigma_squared_d_left > sigma_squared_d_current
            prohibition_flag_left = 1;
        end
        if ( window_size_v(3) >=  max_window_size )
            expanding_priority_flag_left = 0;
        end
    end

    if (prohibition_flag_right ~= 1 && expanding_priority_flag_right == 1 )
        %[sigma_squared_d_right, r_v_right , c_v_right , expanding_priority_flag_right ] = sigma_squared_calculator (r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'right', expanding_priority_flag_right);
        %[sigma_squared_d_right, r_v_right, c_v_right, expanding_priority_flag_right] = ...
        %    sigma_squared_calculator_v1(loc_r, loc_c, r_v, c_v, alpha_f, alpha_d, sigma_two, d_rc, Deriv, ...
        %    sigma_squared_d_denominator, m, n, 'right', expanding_priority_flag_right); %MB July 23, 2019
        [sigma_squared_d_right, r_v_right, c_v_right, expanding_priority_flag_right, window_size_v] = ...
            sigma_squared_calculator_v2(loc_r, loc_c, r_v, c_v, disparity_map, sigma_two, d_rc, Deriv, ...
            sigma_squared_d_denominator, m, n, 'right', expanding_priority_flag_right); %MB July 23, 2019
        if sigma_squared_d_right > sigma_squared_d_current
            prohibition_flag_right = 1;
        end
        if (  window_size_v(4) >=  max_window_size )
            expanding_priority_flag_right = 0;
        end
    end

    prohibition_flag_v = [prohibition_flag_top, prohibition_flag_bottom, prohibition_flag_left, prohibition_flag_right];
    % if one direction gets a zero as its expanding _priority _flag, its sigma squared value
    % needs to be removed from sigma_sq_vec. For example, assume that at the first
    % round, we got the sigma squared value for the top direction and since it
    % is at the border and we can not proceed with its expanstion, it
    % receives zero value for further expantion until the adaptive window reaches to the
    % middle regions for further expansions, in this case we cannot get
    % it's value in the sigma _sq_vec. 
    %Therefore, we need a min value of a direction that was not prohibited and has a expanding priority 
    %equal to one as well.
    expanding_priority_vec = [expanding_priority_flag_top, expanding_priority_flag_bottom, expanding_priority_flag_left, expanding_priority_flag_right];
    % if we get more than one direction to proceed, the one with the
    % minimum sigma squared value would be selected
    [x] = find ( ~prohibition_flag_v & expanding_priority_vec );
    sigma_sq_vec = [ sigma_squared_d_top , sigma_squared_d_bottom , sigma_squared_d_left , sigma_squared_d_right];
    %disp('sigma_sq_vec: top-bottom-left-right');
    %disp(sigma_sq_vec);
    [sigma_squared_d_current, sigma_ind] = min(sigma_sq_vec(x));
    direction_ind = x(sigma_ind);
    
    % updating the location of the window based on the picked direction( the minimum sigma squared value in 4 directions)
    [ r_v, c_v, stop_flag] = update_location(direction_ind, r_v, c_v, r_v_top, c_v_top, r_v_bottom, c_v_bottom, r_v_left, c_v_left, r_v_right, c_v_right);
    r_v_f2 = r_v;
    c_v_f2 = c_v - double(d_rc);
    
%     fig_handle_f1 = findobj('type', 'figure', 'name', 'Image 1');
%     figure(fig_handle_f1);
%     hold on;
%     rectangle('Position',[c_v(1), r_v(1), c_v(2)-c_v(1), r_v(2)-r_v(1)],...
%               'LineWidth',3,'LineStyle','-', 'FaceColor','none');
%     %
%     fig_handle_f2 = findobj('type', 'figure', 'name', 'Image 2');
%     figure(fig_handle_f2);
%     hold on;
%     rectangle('Position',[c_v_f2(1), r_v_f2(1), c_v_f2(2)-c_v_f2(1), r_v_f2(2)-r_v_f2(1)],...
%               'LineWidth',3,'LineStyle','-', 'FaceColor','none');
%     %
%     fig_handle_disparity = findobj('type', 'figure', 'name', 'Disparity Map');
%     figure(fig_handle_disparity);
%     hold on;
%     rectangle('Position',[c_v(1), r_v(1), c_v(2)-c_v(1), r_v(2)-r_v(1)],...
%               'LineWidth',3,'LineStyle','-', 'FaceColor','none');
          
    %[alpha_d, alpha_f] = est_alpha_factors(disparity_map, d_rc, Deriv, r_v, c_v);
    [alpha_d, alpha_f] = est_alpha_factors_v1(disparity_map, d_rc, Deriv, loc_r, loc_c, r_v, c_v);
end

sigma_squared_d = sigma_squared_d_current;

%calculating delta_d , represents incremental correlation of the estimate to be made
%[delta_d] = est_delta_disparity(left_image_mat, right_image_mat, Deriv, sigma_two, alpha_f, alpha_d, r_v, c_v, d_rc);
[delta_d] = est_delta_disparity_v1(left_image_mat, right_image_mat, Deriv, ...
    sigma_two, alpha_f, alpha_d, loc_r, loc_c, r_v, c_v, d_rc); %MB July 23, 2019
end

function [alpha_d, alpha_f] = est_alpha_factors(disparity_map, d_rc, Deriv, r_v, c_v)
%Calculating alfa_d , representing the disparity fluctuations
%Calculating alfa_f, representing the intensity fluctuation

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
        alpha_d = alpha_d +((disparity_map( eta , epsilon) -  d_rc )^ 2) / sqrt (eta^2 + epsilon ^2);
        new_loc = epsilon + d_rc;
        if (new_loc) < 1
            new_loc = 1;
        end
        alpha_f = alpha_f + (Deriv (eta, new_loc) ^ 2);
    end
end
alpha_d = double(alpha_d) / double(N_w);
alpha_f = double(alpha_f)/ double(N_w);
end

function [delta_d] = est_delta_disparity(left_image_mat, right_image_mat, ...
    Deriv, sigma_two, alpha_f, alpha_d, r_v, c_v, d_rc)
%calculating delta_d , represents incremental correlation of the estimate to be made

r_top = r_v(1);
r_bottom = r_v(2);
c_top = c_v(1);
c_bottom = c_v(2);
%
delta_d_numerator = 0;
delta_d_denominator = 0;

for eta = r_top : r_bottom
    for epsilon =  c_top :  c_bottom
        new_loc = epsilon + d_rc;
        if new_loc < 1
            new_loc =1;
        end

        n1 = (left_image_mat(eta, epsilon)- right_image_mat( eta ,new_loc ))* Deriv (eta, new_loc);
        t = round(sqrt (eta ^2 + epsilon^2));
        n2 = sigma_two + ((alpha_f) * (alpha_d) * t) ;
        delta_d_numerator = delta_d_numerator + (n1 /n2);
        delta_d_denominator = delta_d_denominator + ((Deriv ( eta , new_loc ) ^ 2)/ n2) ;
    end 
end
delta_d = double(delta_d_numerator) / double(delta_d_denominator);
end

function [delta_d] = est_delta_disparity_v1(left_image_mat, right_image_mat, ...
    Deriv, sigma_two, alpha_f, alpha_d, loc_r, loc_c, r_v, c_v, d_rc)
%calculating delta_d , represents incremental correlation of the estimate to be made
%MB July 23, 2019

r_top = r_v(1);
r_bottom = r_v(2);
c_top = c_v(1);
c_bottom = c_v(2);
%
delta_d_numerator = 0;
delta_d_denominator = 0;

for eta = r_top : r_bottom
    for epsilon =  c_top :  c_bottom
        %new_loc = epsilon + d_rc;
        new_loc = epsilon - d_rc; %MB July 23, 2019
        if new_loc < 1
            new_loc =1;
        end

        n1 = (left_image_mat(eta, epsilon)- right_image_mat(eta, new_loc))* Deriv(eta, new_loc);
        %t = round(sqrt (eta ^2 + epsilon^2));
        t = round(sqrt ((eta - loc_r)^2 + (epsilon - loc_c)^2)); %MB July 23, 2019
        n2 = sigma_two + ((alpha_f) * (alpha_d) * t) ;
        delta_d_numerator = delta_d_numerator + (n1 /n2);
        delta_d_denominator = delta_d_denominator + ((Deriv(eta , new_loc)^2)/ n2) ;
    end 
end
delta_d = double(delta_d_numerator) / double(delta_d_denominator);
end

function [ r_v, c_v, stop_flag] = update_location(direction_ind, r_v, c_v, r_v_top, c_v_top, r_v_bottom, c_v_bottom, r_v_left, c_v_left, r_v_right, c_v_right  )
stop_flag = 0;
if isempty(direction_ind) 
    direction_ind = 0;
    stop_flag = 1;
end
switch direction_ind
    case 1
        r_v = r_v_top;
        c_v = c_v_top;
        
    case 2
        r_v = r_v_bottom;
        c_v = c_v_bottom;
        
    case 3
        r_v = r_v_left;
        c_v = c_v_left;
        
    case 4
        r_v = r_v_right;
        c_v = c_v_right;
end
end

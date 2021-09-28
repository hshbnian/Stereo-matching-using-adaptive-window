function [delta_d, sigma_squared_d]=...
    disparity_update(left_image_mat, right_image_mat, loc_r, loc_c, disparity_map, window_r, window_c)
%HS 03July2019
%Calculating Sigma_Squared, which is an uncertainty value based on the
%local variation of the intensity and disparity variation within a window
%Takeo Kanade and et al's proposed method 1994
%ref: A stereo matching algorithm with an adaptive window: teory and experiment
% 
% if nargin == 0
%    load Cones_Left.mat
%    left_image_mat = Cones_Left;
%    load Cones_Right.mat
%    right_image_mat = Cones_Right;
%    load disparity_map.mat
%    disparity_map = disparity_map;
%    loc_r = 50;
%    loc_c = 50;
%    window_r = 5;
%    window_c = 4;
% end




% Parameters Initialization 
d_rc = disparity_map (loc_r , loc_c); % constant for a given (r,c))
sigma = 1; 
sigma_two= 2 * (sigma ^ 2);
alfa_f = 0;
alfa_d = 0;
sigma_squared_d_denominator = 0;
delta_d_numerator= 0;
delta_d_denominator= 0;

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
window_size_limit_r=0;
window_size_limit_c=0;

%first order derivative of the right image in x direction
Deriv = diff(right_image_mat);
Deriv(end+1, :)=Deriv(end, :);% since the number of row of the derivative matrix is less in one row, I duplicated the last row

% window adjustment
r_top = loc_r - window_r + 1;
if r_top < 1
   r_top = loc_r ;
end   

r_bottom = loc_r + window_r - 1;
if r_bottom > m
    r_bottom = m;
end

c_top = loc_c - window_c + 1;
if c_top < 1
   c_top = loc_c;
end

c_bottom = loc_c + window_c - 1;
if c_bottom > n
    c_bottom = n;
end

% N_w is a number of damples within the window
N_w = (r_bottom - r_top + 1) * (c_bottom - c_top + 1);


%Calculating alfa_d , representing the disparity fluctuations
%Calculating alfa_f, representing the intensity fluctuation
for eta = r_top : r_bottom
    for epsilon =  c_top :  c_bottom
        alfa_d = alfa_d +((disparity_map( eta , epsilon) -  d_rc )^ 2) / sqrt (eta^2 + epsilon ^2);
        new_loc = epsilon + d_rc;
        if (new_loc)==0
            new_loc = 1;
        end
        alfa_f = alfa_f + (Deriv (eta, new_loc) ^ 2);
    end
end
alfa_d= double(alfa_d) / double(N_w);
alfa_f = double(alfa_f )/ double(N_w);

%calculating delta_d , represents incremental correlation of the estimate to be made
for eta = r_top : r_bottom
    for epsilon =  c_top :  c_bottom
        new_loc = epsilon + d_rc;
        if new_loc == 0
            new_loc =1;
        end
        n1 = (left_image_mat(eta, epsilon)- right_image_mat( eta ,new_loc ))* Deriv (eta, new_loc);
        t = round(sqrt (eta ^2 + epsilon^2));
        n2 = sigma_two + ((alfa_f) * (alfa_d) * t) ;
        delta_d_numerator = delta_d_numerator + (n1 /n2);
        delta_d_denominator = delta_d_denominator + ((Deriv ( eta , new_loc ) ^ 2)/ n2) ;
    end 
end
delta_d = double(delta_d_numerator) / double(delta_d_denominator);

% Coordinates of row and columns 
r_v = [ r_top , r_bottom ] ;
c_v = [ c_top , c_bottom ];

% calculating sigma squared for the current window
sigma_squared_d_current = sigma_squared_calculator_v2(r_v ,c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'current', expanding_priority_flag_center, window_size_limit_r, window_size_limit_c);

%while all prohibitation flags are zeros(if there is still a way to expand )
% all(prohibition_flag_v) means ==> it shows one if all directions got prohibited
while ~(all(prohibition_flag_v) || ((window_size_limit_r >= 15) ||  (window_size_limit_c >= 15 )))
 window_size_limit_r= (r_v(2) - r_v(1));
 window_size_limit_c= (c_v(2) - c_v(1));    
%  fprintf(" window size==> Row:%d:%d Column:%d:%d ==> current uncertainty:%.4f \n ", r_v(1), r_v(2), c_v(1), c_v(2), sigma_squared_d_current );
% fprintf('current uncertainty:%.4f \n',  sigma_squared_d_current);

%calculating sigma squared for the expanded window in 4 directions
    if (prohibition_flag_top ~= 1 && expanding_priority_flag_top ==1)
        [sigma_squared_d_top ,r_v , c_v, expanding_priority_flag_top] = sigma_squared_calculator (r_v, c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'top', expanding_priority_flag_top, window_size_limit_r, window_size_limit_c);
        if sigma_squared_d_top > sigma_squared_d_current
             prohibition_flag_top = 1;
        end
    end

    if (prohibition_flag_bottom ~= 1 && expanding_priority_flag_bottom ==1)
        [sigma_squared_d_bottom ,r_v , c_v , expanding_priority_flag_bottom] = sigma_squared_calculator (r_v, c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'bottom',expanding_priority_flag_bottom, window_size_limit_r, window_size_limit_c);
        if sigma_squared_d_bottom > sigma_squared_d_current
            prohibition_flag_bottom = 1;
        end
    end

    if (prohibition_flag_left ~= 1 && expanding_priority_flag_left ==1 )
        [sigma_squared_d_left, r_v , c_v , expanding_priority_flag_left]= sigma_squared_calculator (r_v, c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'left',  expanding_priority_flag_left, window_size_limit_r, window_size_limit_c);
        if sigma_squared_d_left > sigma_squared_d_current
            prohibition_flag_left = 1;
        end
    end

    if (prohibition_flag_right ~= 1 && expanding_priority_flag_right ==1 )
        [sigma_squared_d_right, r_v , c_v , expanding_priority_flag_right ] = sigma_squared_calculator (r_v, c_v, alfa_f, alfa_d, sigma_two, d_rc, Deriv, sigma_squared_d_denominator, m, n, 'right', expanding_priority_flag_right, window_size_limit_r, window_size_limit_c);
        if sigma_squared_d_right > sigma_squared_d_current
            prohibition_flag_right = 1;
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
    [sigma_squared_d_current, sigma_ind] = min(sigma_sq_vec(x));
    %direction_ind = x(sigma_ind);
end
sigma_squared_d = sigma_squared_d_current;
end
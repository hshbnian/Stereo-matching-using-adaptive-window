function test_disparity_map_using_updating_window

%% Loading left, right, and disparity map matrices
load('Cones_Left.mat', 'Cones_Left' );
left_image_mat = Cones_Left;
load('Cones_Right.mat', 'Cones_Right');
right_image_mat = Cones_Right;
load('disparity_map.mat', 'disparity_map')
disparity_map=disparity_map;

% Debug: input images
figure('name', 'Image 1');
imshow(left_image_mat, []); colorbar;
impixelinfo;
%
figure('name', 'Image 2');
imshow(right_image_mat, []); colorbar;
impixelinfo;
%% Initialization
[m, n] = size(disparity_map);
no_of_itn= 1;

%% calling disparity_update funcion to calculate delta_d value
figure('name', 'Disparity Map');
no_of_figs = no_of_itn + 1;
no_of_fig_col = 1;
no_of_fig_row = ceil(no_of_figs / no_of_fig_col);
%
subplot(no_of_fig_row, no_of_fig_col, 1);
imshow(disparity_map,[0 55]);colorbar; colormap jet; title('Initial Disparity');
impixelinfo;
%
for iter = 1 : no_of_itn
    for loc_r = 1 : m  %208%132%78%214 : 215%
        fprintf('Itn: %d; row: %d / %d\n', iter, loc_r, m); 
        for loc_c = 1 : n %194%320%218%77 : 78%
           %fprintf('Itn: %d; row: %d / %d col: %d / %d\n', iter, loc_r, m, loc_c, n); 
           [delta_d, sigma_squared_d]= disparity_update_v4(left_image_mat, right_image_mat, loc_r, loc_c, disparity_map, 1, 1);
           delta_d_mat(loc_r,loc_c)= delta_d;
        end 
    end
    %fprintf('Out of loop\n');
    %%  Updating the initial disparity map using the calculated delta_d matrix
    %disparity_map_updated= double(disparity_map) + double(delta_d_mat);
    disparity_map_updated= double(disparity_map) - double(delta_d_mat); %MB July 23, 2019
    %
    fig_handle = findobj('type', 'figure', 'name', 'Disparity Map');
    figure(fig_handle);
    subplot(no_of_fig_row, no_of_fig_col, iter + 1);
    imshow(disparity_map_updated,[0 55]);colorbar; colormap jet; title(sprintf('Iteration %d', iter));
    disparity_map = uint8 (disparity_map_updated);
    %fprintf('Updated figure; itn: %d\n', iter);
end

end
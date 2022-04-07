function [sigma2_vect,theta_vect] = generate_random_angle_rotations(time_index_index,store_registration)

% this outputs the goodness of alignment measure for each of the random
% rotations (stored in theta_vect)

% 100 can be some variable instead.

sigma2_vect = zeros(100, 1);
theta_vect = zeros(100, 3);

for which_rot = 1:100
    % store this time index
    time_index = valid_time_indices(time_index_index);
    
    % store next in series
    time_index_plus_1 = valid_time_indices(time_index_index+1);
    
    % store combined image for both.
    A = imread([name_of_embryo,num2str(time_index,'%05.5d'),suffix_for_embryo],1);
    tiff_info = imfinfo([name_of_embryo,num2str(time_index,'%05.5d'),suffix_for_embryo]);
    % combine all tiff stacks into 1 3D image.
    combined_image = zeros(size(A,1), size(A,2), size(tiff_info, 1));
    for j = 1:size(tiff_info, 1)
        A = imread([name_of_embryo,num2str(time_index,'%05.5d'),suffix_for_embryo],j);
        combined_image(:,:,j) = A(:,:,1);
    end
    combined_image1 = combined_image;
    
    resXY = 0.208;
    resZ = 2.0;
    reduceRatio = 1/4;
    combined_image1 = isotropicSample_nearest(double(combined_image1), resXY, resZ, reduceRatio);
    
    A = imread([name_of_embryo,num2str(time_index_plus_1,'%05.5d'),suffix_for_embryo],1);
    tiff_info = imfinfo([name_of_embryo,num2str(time_index_plus_1,'%05.5d'),suffix_for_embryo]);
    % combine all tiff stacks into 1 3D image.
    combined_image = zeros(size(A,1), size(A,2), size(tiff_info, 1));
    for j = 1:size(tiff_info, 1)
        A = imread([name_of_embryo,num2str(time_index_plus_1,'%05.5d'),suffix_for_embryo],j);
        combined_image(:,:,j) = A(:,:,1);
    end
    combined_image2 = combined_image;
    
    resXY = 0.208;
    resZ = 2.0;
    reduceRatio = 1/4;
    combined_image2 = isotropicSample_nearest(double(combined_image2), resXY, resZ, reduceRatio);
    
    % STORE MESHGRID
    [X, Y, Z] = meshgrid(1:size(combined_image1, 2), 1:size(combined_image1, 1), 1:size(combined_image1, 3));
    
    % FRACTION OF POINTS (DOWNSAMPLING)
    fraction_of_selected_points = 1/10;
    %
    find1 = find(combined_image1(:)~=0);
    number_of_points = length(find1);
    
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find1 = find1(p);
    
    ptCloud1 = [X(find1), Y(find1), Z(find1)];
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
    %
    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);
    
    ptCloud2 = [X(find2), Y(find2), Z(find2)];
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [mean(X(find2)), mean(Y(find2)), mean(Z(find2))];
    
    ptCloud2 = pointCloud(ptCloud2);
    theta1 =rand*360;
    rot1 = [ cosd(theta1) -sind(theta1) 0; ...
        sind(theta1) cosd(theta1) 0; ...
        0 0 1];
    theta2 =rand*360;
    rot2 = [ 1 0 0; ...
        0 cosd(theta2) -sind(theta2); ...
        0 sind(theta2) cosd(theta2)];
    theta3 =rand*360;
    rot3 = [ cosd(theta3) 0 sind(theta3); ...
        0 1 0; ...
        -sind(theta3) 0 cosd(theta3)];
    tform = rigid3d(rot1*rot3*rot2,[0,0,0]);
    ptCloud2 = pctransform(ptCloud2,tform);
    ptCloud2 = ptCloud2.Location;
    theta_vect(which_rot, 1) = theta1;
    theta_vect(which_rot, 2) = theta2;
    theta_vect(which_rot, 3) = theta3;
    
    % Example 3. 3D Rigid CPD point-set registration. Full options intialization.
    %  3D face point-set.
    
    %         % Init full set of options %%%%%%%%%%
    %         opt.method='nonrigid'; % use nonrigid registration
    %         opt.beta=2;            % the width of Gaussian kernel (smoothness)
    %         opt.lambda=3;          % regularization weight
    %
    %         opt.viz=1;              % show every iteration
    %         opt.outliers=0;         % don't account for outliers
    %         opt.fgt=0;              % do not use FGT (default)
    %         opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
    %         opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
    %
    %         opt.max_it=100;         % max number of iterations
    %         opt.tol=1e-8;           % tolerance
    %
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Set the options
    opt.method='rigid'; % use rigid registration
    opt.viz=1;          % show every iteration
    opt.outliers=0;     % do not assume any noise
    
    opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
    opt.scale=0;        % estimate global scaling too (default)
    opt.rot=1;          % estimate strictly rotational matrix (default)
    opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
    
    opt.max_it=200;     % max number of iterations
    opt.tol=1e-5;       % tolerance
    
    
    % registering Y to X
    [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2,opt);
    sigma2_vect(which_rot) = sigma2;
    store_registration{time_index_index, 1} = Transform;
    
    %figure; hold all; title('Before'); cpd_plot_iter(ptCloud1, ptCloud2);
    figure; hold all;
    %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
    cpd_plot_iter(ptCloud1, ptCloud2);
    figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
    
    close all;
    
    [iou_matrix, M, corresponding_ious_for_matches, ...
        cell_labels_I_care_about1, cell_labels_I_care_about2, ...
        center_point_for_each_label1, center_point_for_each_label2, ...
        match_based_on_nearest_neighbors, ~, ~, ...
        alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,ptCloud1,...
        combined_image1,combined_image2,find1,find2);
    
    % DISPLAY CURRENT TIME INDEX (just to make sure that calculation is not stalled).
    disp(time_index);
    disp(which_rot);
end


end


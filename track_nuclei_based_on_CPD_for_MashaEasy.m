
function [] = track_nuclei_based_on_CPD_for_Maddy(  )
% Really barebones code for Maddy to learn image registration...
% Loop through and convince yourself that you are happy with the registrations it is giving you

% see https://sites.google.com/site/myronenko/research/cpd for this
% specific implementation of CPD

% see https://ieeexplore.ieee.org/document/5432191 for explanation of what
% code is doing

% Maddy: Change these paths 
% (make sure that these are right for your local environment)
addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/core'));
addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/data'));

%%

% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% Which image indices to run over...
which_number_vect = 0:200;


% What is the prefix for the embryo names?
%name_of_embryo = '/Users/hnunley/Desktop/test_for_maddy/OG_st0_NANOGGATA6_2105_better_model/Stardist3D_klbOut_Cam_Long_';
name_of_embryo = '/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/OG_st0_NANOGGATA6/Stardist3D_klbOut_Cam_Long_';

% Suffix: yours is probably '.lux.tif'
suffix_for_embryo = '.tif';

% Initialize empty graph and cell array for storing registration
valid_time_indices = which_number_vect;
store_registration = cell((length(valid_time_indices)-1), 1);
sigma2_vect_saved = zeros((length(valid_time_indices)-1), 1);
% also, check the alignment of this one with the time frame after
%for time_index_index = 146:(length(valid_time_indices)-1)
for time_index_index = 132:133
    
    % store this time index
    time_index = valid_time_indices(time_index_index)
    
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
    fraction_of_selected_points =  1/10;  % slow to run at full scale - but make full res points and xform?
    %
    find1 = find(combined_image1(:)~=0);  % this is the indices into combined_image1 to get indices into (X,Y,Z) to the full set of point
    number_of_points = length(find1);
    
    % why random points - why not just subsample by 10 ?
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    full_ptCloud1 =  [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
    find1 = find1(p);
    
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
    %
    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);
    
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [mean(X(find2)), mean(Y(find2)), mean(Z(find2))];
    ptCloud2 = pointCloud(ptCloud2);
    
    tform = rigid3d(eye(3), [0,0,0]);
    
    ptCloud2 = pctransform(ptCloud2,tform);
    ptCloud2 = ptCloud2.Location;
    
    % Example 3. 3D Rigid CPD point-set registration. Full options intialization.
    %  3D face point-set.
    
    % Set the options
    opt.method='rigid'; % use rigid registration
    opt.viz=1;          % show every iteration
    opt.outliers=0;     % do not assume any noise
    
    opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
    opt.scale=0;        % estimate global scaling too (default)
    opt.rot=1;          % estimate strictly rotational matrix (default)
    opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
    
    opt.max_it=100;     % max number of iterations
    opt.tol=1e-5;       % tolerance
      
    % registering Y to X
   % [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2,opt);
    [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2,opt);
    % run transform on full_ptCloud1
    
    store_registration{time_index_index, 1} = Transform;
    % sigma2 is a measure of goodness-of-alignment between the two point clouds
    % see 
    sigma2_vect_saved(time_index_index, 1) = sigma2;
    
    figure; hold all;
    cpd_plot_iter(ptCloud1, ptCloud2);
    figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
    
    pause;
    close all;
    
    [iou_matrix, M, corresponding_ious_for_matches, ...
        cell_labels_I_care_about1, cell_labels_I_care_about2, ...
        center_point_for_each_label1, center_point_for_each_label2, ...
        match_based_on_nearest_neighbors, ~, ~, ...
        alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,ptCloud1,...
        combined_image1,combined_image2,find1,find2);
    store_matches{time_index_index, 1} = M
    store_iou_table{time_index_index, 1} = iou_matrix;
        
    % loop through all matches
    % LB: M is nMatches x 2 (label at time t, label at time t+1)  

    for i = 1:length(M)
        find_non_zero_please = find(iou_matrix(M(i,1),:));
        if (length(find_non_zero_please) > 1)  % more than one iou match
            figure; hold all; cpd_plot_iter(ptCloud1, Transform.Y);
            %hold all; plot(alpha_shape_for_each_label1{M(i,1),1},'FaceColor','red','FaceAlpha',0.5);
            for j = 1:length(find_non_zero_please)
                if (find_non_zero_please(j) == M(i,2)) % this is the best match?
                    hold all; plot(alpha_shape_for_each_label2{M(i,2),1},'FaceColor','green','FaceAlpha',0.5);
                else
                    hold all; plot(alpha_shape_for_each_label2{find_non_zero_please(j),1},'FaceColor','black','FaceAlpha',0.5);
                end
                
            end
            
            title([num2str(corresponding_ious_for_matches(i)),';',num2str(i)]);
            pause;
            close all;
        end
    end
    
    % DISPLAY CURRENT TIME INDEX (just to make sure that calculation is not stalled).
    disp('time index');
    disp(time_index);
end

% Save vector of transformations...
save('transform_labels_pCloud.mat', 'store_registration');
save('matches.mat','store_matches');
save('iou_table.mat','store_iou_table');


end


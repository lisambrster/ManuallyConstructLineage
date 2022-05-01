
%function [] = PrecomputeRegistrationTransforms(  )
% uses already tracked embryo that was easy that Masha gave to me
% performs 'easy' registration
% creates graph
% shows some visualizations..
%% User Inputs
addpath(genpath('oldfiles'));
addpath(genpath('CPD2'));

verbosemode = 1; 

% What is the prefix for the embryo names?
name_of_embryo = '/home/jotsund/Desktop/tracking_test/220316_out/A/st9/klbOut_Cam_Long_';
% Suffix: yours is probably '.lux.tif'
suffix_for_embryo = '.lux.label.tif';

RegistrationFileName = 'transformsFirstFive.mat';

G_based_on_nn = graph;
% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% how many random orientations do you want - minimum.
maxItr = 100;

% Which image indices to run over...
which_number_vect = 1:100; %what is the valid time range.

firstTime = 1;
lastTime = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SCRIPT BEGINS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valid_time_indices = which_number_vect;
% Initialize empty graph and cell array for storing registration
store_registration = cell((length(valid_time_indices)-1), 1);
sigma2_vect_saved = zeros((length(valid_time_indices)-1), 1);



% Set the options for CPD - these are always the same.
opt.method='rigid'; % use rigid registration
opt.viz=0;          % show every iteration
opt.outliers=0;     % do not assume any noise

opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
opt.scale=0;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

opt.max_it=200;     % max number of iterations
opt.tol=1e-5;       % tolerance
opt.ftg = 1; % make faster

%% Note: last time point will look at this time point and the next one

for time_index_index = firstTime:lastTime  % time 78 takes a LONG time
    sigma2tests = zeros(maxItr,1)*nan;  % these are the 100 tests for one pair of images
    transforms = cell(maxItr,1);
    tic
    % store this time index
    time_index = valid_time_indices(time_index_index);
    disp(time_index)

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
    rng(1)
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    %full_ptCloud1 =  [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
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


    % Example 3. 3D Rigid CPD point-set registration. Full options intialization.
    %  3D face point-set.
    bStoredReg = false;
    if bStoredReg
        Transform = xforms.store_registration{time_index_index,1};
    else

        %         nobody is using this, it's been hardcoded to off, so we'll remove it
                % bEasy = false;
        %         if (bEasy)
        %             ptCloud2 = pctransform(ptCloud2,tform); % this makes ptCloud a pointCloud
        %             ptCloud2 = ptCloud2.Location;
        %             % Set the options
        %             opt.method='rigid'; % use rigid registration
        %             opt.viz=0;          % show every iteration
        %             opt.outliers=0;     % do not assume any noise
        %
        %             opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
        %             opt.scale=0;        % estimate global scaling too (default)
        %             opt.rot=1;          % estimate strictly rotational matrix (default)
        %             opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
        %
        %             opt.max_it=100;     % max number of iterations
        %             opt.tol=1e-5;       % tolerance
        %             opt.fgt =1; % should make it faster but maybe I'm not set up right?
        %
        %             % registering Y to X
        %             % [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2,opt);
        %             [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2,opt);
        %             % run transform on full_ptCloud1
        %
        %             store_registration{time_index_index, 1} = Transform;
        %             % sigma2 is a measure of goodness-of-alignment between the two point clouds
        %             % see
        %             sigma2_vect_saved(time_index_index, 1) = sigma2;
        %
        %
        %             figure;
        %             hold on;
        %             cpd_plot_iter(ptCloud1, ptCloud2);
        %             hold off;
        %
        %             figure;
        %             hold on;
        %             title('After registering Y to X.');
        %             cpd_plot_iter(ptCloud1, Transform.Y);
        %             hold off;
        %
        %
        %
        %         else


        sigma2_vect = zeros(100, 1);
        %theta_vect = zeros(100, 3);

        which_rot = 1;
        min_sigma2 = 100;
        counter = 0;
    
        while ((min_sigma2 > 10) && (counter < maxItr))
            fprintf('.')
            thetaHoldler = {};
            for whichrot = 1:gcp().NumWorkers
                rng((1+counter)*whichrot)%reseeding RNG call to ensure each rand is different below... 
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
                ptCloud3 = pctransform(ptCloud2,tform);
                ptCloud2Loc = ptCloud3.Location;
                %thetaz = [theta1,theta2,theta3];
                %thetaHoldler{whichrot} = thetaz;

                % registering Y to X
                [Transform,~, sigma2]=cpd_register(ptCloud1,ptCloud2Loc,opt);
                %sigma2_vect(which_rot) = sigma2;
                %store_registration{counter+time_index_index, 1} = Transform;
                transforms{counter+whichrot, 1} = Transform;
                sigma2tests(counter+whichrot) = sigma2;
                %which_rot = which_rot + 1;

                %figure; hold all; title('Before'); cpd_plot_iter(ptCloud1, ptCloud2);
                % figure; hold all;
                %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
                %cpd_plot_iter(ptCloud1, ptCloud2Loc);
                %figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
            end
%%%turns out we didn't need the thetas after all
%             for aye = 1:11 %resynchronize thetas - this could be simplified and it would save some millisecond/loop. 
%                 theta_vect(counter+aye, 1) = thetaHoldler(aye,1);
%                 theta_vect(counter+aye, 2) = thetaHoldler(aye,2);
%                 theta_vect(counter+aye, 3) = thetaHoldler(aye,3);
%             end
            %end
            counter = counter + gcp().NumWorkers;%which_rot;
            if which_rot > 99
                disp('did not find transformation with sigma2 < 10');
            end
            % get the best one we found this loop, 
            [min_sigma2, min_ind] = nanmin(sigma2tests);
            Transform = transforms{min_ind,1};
            if verbosemode
                disp(min_sigma2); disp(min_ind);
            end
            store_registration{time_index_index, 1} = Transform;

        end
        % put a check point here

        % need to pick the one with the lowest sigma
        %index_min = find(sigma2_vect ==min(sigma2_vect));
        %disp(['min sigma ',sigma2_vect(index_min)]);
        %close all;
        %figure; hold all; title('Before registration'); cpd_plot_iter(ptCloud1, ptCloud2Loc);

        %figure; hold all;
        %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
        %cpd_plot_iter(ptCloud1, ptCloud2Loc);

        %figure;
        %hold on;
        %title('After registration');
        %cpd_plot_iter(ptCloud1, Transform.Y);
        %hold off;

    end
    disp(' ')
    toc
end

% Save vector of transformations...
save(RegistrationFileName, 'store_registration');
%save('matches.mat','store_matches');
%save('iou_table.mat','store_iou_table');
%save('graph.mat1to3','G_based_on_nn');

%% example to add node and edge after split
%G_based_on_nn = G_based_on_nn.addnode('078_016');
%G_based_on_nn = G_based_on_nn.addedge('077_010','078_016')
%end


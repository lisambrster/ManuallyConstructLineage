
function [] = track_nuclei_based_on_CPD_EasyFromMasha_make_graph(  )
% uses already tracked embryo that was easy that Masha gave to me
% performs 'easy' registration
% creates graph
% shows some visualizations..

addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/core'));
addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/data'));

%%

G_based_on_nn = graph;

% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% Which image indices to run over...
which_number_vect = 26:100;


% What is the prefix for the embryo names?
%name_of_embryo = '/Users/hnunley/Desktop/test_for_maddy/OG_st0_NANOGGATA6_2105_better_model/Stardist3D_klbOut_Cam_Long_';
%  stack 1 in 211018. 
% from here: /mnt/home/mavdeeva/ceph/mouse/data/segmentation_out/211018. /stack1_channel_2_obj_left/nuclear
%name_of_embryo = '/mnt/ceph/users/mavdeeva/mouse/data/segmentation_out/211018/stack_1_channel_2_obj_left/nuclear/Stardist3D_klbOut_Cam_Long_';
name_of_embryo = '/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/211018_stack_1/Stardist3D_klbOut_Cam_Long_';
name_of_embryo = '/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/211018_stack_1/Stardist3D_klbOut_Cam_Long_';


% Suffix: yours is probably '.lux.tif'
suffix_for_embryo = '.lux.tif';

% Initialize empty graph and cell array for storing registration
valid_time_indices = which_number_vect;
store_registration = cell((length(valid_time_indices)-1), 1);
sigma2_vect_saved = zeros((length(valid_time_indices)-1), 1);
sigma2tests = zeros(100,1);  % these are the 100 tests for one pair of images
transforms = cell(100,1);
% also, check the alignment of this one with the time frame after
%for time_index_index = 146:(length(valid_time_indices)-1)
for time_index_index = 22:100
     
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
    
    bEasy = false;
    if (bEasy)
        ptCloud2 = pctransform(ptCloud2,tform); % this makes ptCloud a pointCloud
        ptCloud2 = ptCloud2.Location;
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


    
    else
        
        
        sigma2_vect = zeros(100, 1);
        theta_vect = zeros(100, 3);

        which_rot = 1;
        sigma2 = 100;
        while ((sigma2 > 10) && (which_rot < 101))    
        
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
            ptCloud2Loc = ptCloud2.Location;
            theta_vect(which_rot, 1) = theta1;
            theta_vect(which_rot, 2) = theta2;
            theta_vect(which_rot, 3) = theta3;
     
            % Set the options
            opt.method='rigid'; % use rigid registration
            opt.viz=0;          % show every iteration
            opt.outliers=0;     % do not assume any noise

            opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
            opt.scale=0;        % estimate global scaling too (default)
            opt.rot=1;          % estimate strictly rotational matrix (default)
            opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

            opt.max_it=200;     % max number of iterations
            opt.tol=1e-5;       % tolerance


            % registering Y to X
            [Transform, ~, sigma2]=cpd_register(ptCloud1,ptCloud2Loc,opt);
            %sigma2_vect(which_rot) = sigma2;
            store_registration{time_index_index, 1} = Transform;
            transforms{which_rot, 1} = Transform;
            sigma2tests(which_rot) = sigma2;
            which_rot = which_rot + 1;

            %figure; hold all; title('Before'); cpd_plot_iter(ptCloud1, ptCloud2);
           % figure; hold all;
            %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
            %cpd_plot_iter(ptCloud1, ptCloud2Loc);
            %figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
        end
        
        if which_rot > 99
            disp('did not find transformation with sigma2 < 10');
            % get the best one we found
            [min_sigma2, min_ind] = min(sigma2tests) 
            Transform = transforms{min_ind,1};
        end
            % put a check point here
        
        % need to pick the one with the lowest sigma
        %index_min = find(sigma2_vect ==min(sigma2_vect));
        %disp(['min sigma ',sigma2_vect(index_min)]);
        %figure; hold all; title('Before'); cpd_plot_iter(ptCloud1, ptCloud2);
        close all;
        %figure; hold all;
        %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
        %cpd_plot_iter(ptCloud1, ptCloud2Loc);
        figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
   
        
    end
    
    [iou_matrix, M, corresponding_ious_for_matches, ...
            cell_labels_I_care_about1, cell_labels_I_care_about2, ...
            center_point_for_each_label1, center_point_for_each_label2, ...
            match_based_on_nearest_neighbors, ~, ~, ...
            alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,ptCloud1,...
            combined_image1,combined_image2,find1,find2);
     store_matches{time_index_index, 1} = M;
     store_iou_table{time_index_index, 1} = iou_matrix;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make the graph..
        
    [nn nd]=kNearestNeighbors(center_point_for_each_label1, center_point_for_each_label2,min(3,length(center_point_for_each_label2)));
    if length(nn(:,1))~=length(unique(nn(:,1))) % Reject duplicate nearest neighbors
        dup=find_duplicates(nn(:,1));
        for lvd=1:size(dup,1)
            [ic,ia,ib]=intersect(nn(dup(lvd).ind,2),setdiff(1:size(nn,1),nn(:,1)));
            if ~isempty(ia)
                nn(dup(lvd).ind(ia),1)=nn(dup(lvd).ind(ia),2);
                nd(dup(lvd).ind(ia),1)=nd(dup(lvd).ind(ia),2);
                loi=dup(lvd).ind(setdiff(1:size(dup(lvd).ind,1),ia)); % treat triple and more entries
                if length(loi)>1
                    [mv,mi]=min(nd(loi,1));
                    loi1=setdiff(1:length(loi),mi);
                    nn(loi(loi1),1)=NaN;
                end
            else
                [mv mi]=min(sum(nd(dup(lvd).ind,1),2));
                loi=setdiff(dup(lvd).ind,dup(lvd).ind(mi));
                nn(loi,1)=NaN;
            end
        end
    end
    
    nn=nn(:,1);
    nd=nd(:,1);
    
    sample_graph = graph;
    for iind = 1:length(cell_labels_I_care_about1)
        this_label = cell_labels_I_care_about1(iind);
        
        
        % store node props table... so that node can be added with volume
        NodePropsTable = table({[num2str(time_index,'%05.3d'),'_', num2str(this_label,'%05.3d')]}, center_point_for_each_label1(iind, 1), center_point_for_each_label1(iind, 2), center_point_for_each_label1(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph = addnode(sample_graph, NodePropsTable);
        
    end
    
    for iind = 1:length(cell_labels_I_care_about2)
        this_label = cell_labels_I_care_about2(iind);
        
        
        % store node props table... so that node can be added with volume
        NodePropsTable = table({[num2str(time_index_plus_1,'%05.3d'),'_', num2str(this_label,'%05.3d')]}, center_point_for_each_label2(iind, 1), center_point_for_each_label2(iind, 2), center_point_for_each_label2(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph = addnode(sample_graph, NodePropsTable);
        
    end
    
    for point_index = 1:length(nn)
        
        if (~isnan(nn(point_index)))
            % make directed edges (in time) between matches + store iou for the match as a graph weight
            sample_graph = addedge(sample_graph, [num2str(time_index,'%05.3d'),'_', num2str(cell_labels_I_care_about1(nn(point_index)),'%05.3d')],...
                [num2str(time_index_plus_1,'%05.3d'),'_', num2str(cell_labels_I_care_about2(point_index),'%05.3d')]);
        end
        if (~isnan(nn(point_index)))
            
            % make directed edges (in time) between matches + store iou for the match as a graph weight
            G_based_on_nn = addedge(G_based_on_nn, [num2str(time_index,'%05.3d'),'_', num2str(cell_labels_I_care_about1(nn(point_index)),'%05.3d')],...
                [num2str(time_index_plus_1,'%05.3d'),'_', num2str(cell_labels_I_care_about2(point_index),'%05.3d')]);
            
        end
        
    end
    % visualization for checking if everything is correct
    hold all; plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0);

    figure; h1 = plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name);
    disp(time_index);


        
    % loop through all matches
    % LB: M is nMatches x 2 (label at time t, label at time t+1)  

    for i = 1:length(M)
        find_non_zero_please = find(iou_matrix(M(i,1),:));
        if (length(find_non_zero_please) > 1)  % more than one iou match
            figure; hold all; cpd_plot_iter(ptCloud1, Transform.Y);
            hold all; plot(alpha_shape_for_each_label1{M(i,1),1},'FaceColor','red','FaceAlpha',0.5);
            for j = 1:length(find_non_zero_please)
                if (find_non_zero_please(j) == M(i,2)) % this is the best match?
                    hold all; plot(alpha_shape_for_each_label2{M(i,2),1},'FaceColor','green','FaceAlpha',0.5);
                else
                    hold all; plot(alpha_shape_for_each_label2{find_non_zero_please(j),1},'FaceColor','black','FaceAlpha',0.5);
                end
                
            end
            
            title([num2str(corresponding_ious_for_matches(i)),';',num2str(i)]);
            %pause;
            %close all;
        end
    end
    

    %% this is good place for breakpoint to check that all nuclei are matched
    figure;
    plot(G_based_on_nn,'layout','layered')
    disp('time index');
    disp(time_index);
    %pause;
    close all;
end

% Save vector of transformations...
save('transforms.mat', 'store_registration');
save('matches.mat','store_matches');
save('iou_table.mat','store_iou_table');
save('graph.mat','G_based_on_nn');

%% example to add node and edge after split
%G_based_on_nn = G_based_on_nn.addnode('078_015');
%G_based_on_nn = G_based_on_nn.addedge('077_010','078_015')
end


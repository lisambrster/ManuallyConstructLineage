
% testing out Hayden's code

%G_based_on_nn
name_of_embryo = 
suffix_for_embryo
store_registration


% read in two label images

% store this time index
time_index = 1; %valid_time_indices(time_index_index);

% store next in series
time_index_plus_1 = 2; %valid_time_indices(time_index_index+1);

% store combined image for both.
suffix_for_embryo1 = 
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

tform = rigid3d(eye(3), [0,0,0]);


ptCloud2 = pctransform(ptCloud2,tform);
ptCloud2 = ptCloud2.Location;

%
ptCloud2 = cpd_transform(ptCloud2, store_registration{time_index_index,1});
%ptCloud2 = ptCloud2;

% % Set the options
% opt.method='rigid'; % use rigid registration
% opt.viz=1;          % show every iteration
% opt.outliers=0;     % do not assume any noise
% 
% opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
% opt.scale=0;        % estimate global scaling too (default)
% opt.rot=1;          % estimate strictly rotational matrix (default)
% opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)
% 
% opt.max_it=100;     % max number of iterations
% opt.tol=1e-6;       % tolerance
% 
% % registering Y to X
% [Transform, ~]=cpd_register(ptCloud1,ptCloud2,opt);
% 
% %
% figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
% 
% %pause;
% 
% [iou_matrix, M, corresponding_ious_for_matches, ...
%     cell_labels_I_care_about1, cell_labels_I_care_about2, ...
%     center_point_for_each_label1, center_point_for_each_label2, ...
%     match_based_on_nearest_neighbors, ~, ~, ...
%     alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,...
%     ptCloud1,...
%     combined_image1,combined_image2,find1,find2);


figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, ptCloud2);

%pause;

[iou_matrix, M, corresponding_ious_for_matches, ...
    cell_labels_I_care_about1, cell_labels_I_care_about2, ...
    center_point_for_each_label1, center_point_for_each_label2, ...
    match_based_on_nearest_neighbors, ~, ~, ...
    alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(ptCloud2,...
    ptCloud1,...
    combined_image1,combined_image2,find1,find2);


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

for iind1 = 1:length(cell_labels_I_care_about1)
    
    for iind2 = 1:length(cell_labels_I_care_about2)
        
        this_label = cell_labels_I_care_about1(iind1);
        nodeID1 = [num2str(time_index,'%05.3d'),'_', num2str(this_label,'%05.3d')];
        this_label = cell_labels_I_care_about2(iind2);
        nodeID2 = [num2str(time_index_plus_1,'%05.3d'),'_', num2str(this_label,'%05.3d')];
        
        if (findnode(G_based_on_nn, nodeID1) ~= 0 && findnode(G_based_on_nn, nodeID2) ~= 0)
            if (findedge(G_based_on_nn,nodeID1,nodeID2) ~= 0)
                % make directed edges (in time) between matches + store iou for the match as a graph weight
                sample_graph = addedge(sample_graph, nodeID1,nodeID2);
            end
        end
    end
    
end

hold all; plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0);

figure; h1 = plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name);
disp(time_index);


end


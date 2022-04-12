function [outputArg1,outputArg2] = compute_G_based_on_nn_given_rotations(inputArg1,inputArg2)


G_based_on_nn = graph;
% also, check the alignment of this one with the time frame after
for time_index_index = 1:(length(valid_time_indices)-1)
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
    
    tform = rigid3d(eye(3), [0,0,0]);
    
%     if (time_index_index == 145)
%         theta1 = 60.1580;
%         rot1 = [ cosd(theta1) -sind(theta1) 0; ...
%             sind(theta1) cosd(theta1) 0; ...
%             0 0 1];
%         
%         theta2 = 98.1678;
%         rot2 = [ 1 0 0; ...
%             0 cosd(theta2) -sind(theta2); ...
%             0 sind(theta2) cosd(theta2)];
%         
%         theta3 =103.0349;
%         rot3 = [ cosd(theta3) 0 sind(theta3); ...
%             0 1 0; ...
%             -sind(theta3) 0 cosd(theta3)];
%         tform = rigid3d(rot1*rot3*rot2,[0,0,0]);
%     end
    
    ptCloud2 = pctransform(ptCloud2,tform);
    ptCloud2 = ptCloud2.Location;
    
    %Transform = store_registration{time_index_index,1};
    [iou_matrix, M, corresponding_ious_for_matches, ...
        cell_labels_I_care_about1, cell_labels_I_care_about2, ...
        center_point_for_each_label1, center_point_for_each_label2, ...
        match_based_on_nearest_neighbors, ~, ~, ...
        alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(cpd_transform(ptCloud2, store_registration{time_index_index,1}),...
        ptCloud1,...
        combined_image1,combined_image2,find1,find2);
    
    
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
    disp(time_index);
end



end


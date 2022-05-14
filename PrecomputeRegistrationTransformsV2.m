
%function [] = PrecomputeRegistrationTransforms(  )
% uses already tracked embryo that was easy that Masha gave to me
% performs 'easy' registration
% creates graph
% shows some visualizations..
%% User Inputs

verbosemode = 1;  %periodically output the minimum sigma2 value.
bRerun = false;
poolsize = 30; % set the number of workers in your parallel pool


addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/core'));
addpath(genpath('/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/CPD2/data'));

% What is the prefix for the embryo names?
name_of_embryo = '/Users/lbrown/Documents/PosfaiLab/Registration/HaydensReg2022/211018_stack_1/Stardist3D_klbOut_Cam_Long_';
% Suffix: yours is probably '.lux.tif'
suffix_for_embryo = '.lux.tif';


name_of_embryo = '/Users/lbrown/Documents/PosfaiLab/3DStardist/GataNanog/HaydenJan22Set/Stardist3D_klbOut_Cam_Long_';
% Stardist3D_klbOut_Cam_Long_00084.tif
suffix_for_embryo = '.tif';

data_path = '/Users/lbrown/Documents/PosfaiLab/3DStardist/GataNanog/HaydenJan22Set/';
RegistrationFileName = strcat(data_path,'transforms100_150Match.mat');
%InitFileName = strcat(data_path,'init1_10.mat')du
firstTime = 100;
lastTime =  150;

G_based_on_nn = graph;
% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% how many random orientations do you want - minimum.
maxItr = 8;

% Which image indices to run over...
which_number_vect = 1:200; %what is the valid time range.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SCRIPT BEGINS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolsize = 8;


valid_time_indices = which_number_vect;
% Initialize empty graph and cell array for storing registration
store_registration = cell((length(valid_time_indices)-1), 1);
store_init_transform = cell((length(valid_time_indices)-1), 1);
store_sigma = cell((length(valid_time_indices))-1,1);
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

%for time_index_index = firstTime:lastTime  % time 78 takes a LONG time
time_index_index = firstTime;
while time_index_index <= lastTime
    
    sigma2tests = zeros(maxItr,1)*nan;  % these are the 100 tests for one pair of images
    transforms = cell(maxItr,1);
    tic
    % store this time index
    time_index = valid_time_indices(time_index_index);
    disp(time_index)

    % store next in series
    time_index_plus_1 = valid_time_indices(time_index_index+1);

    if (~bRerun)
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
    end

    % STORE MESHGRID
    [X, Y, Z] = meshgrid(1:size(combined_image1, 2), 1:size(combined_image1, 1), 1:size(combined_image1, 3));

    % FRACTION OF POINTS (DOWNSAMPLING)
    fraction_of_selected_points =  1/10;  % slow to run at full scale - but make full res points and xform?
    %
    find1 = find(combined_image1(:)~=0);  % this is the indices into combined_image1 to get indices into (X,Y,Z) to the full set of point
    number_of_points = length(find1);
    
    % compute mean BEFORE random subsampling of points - so constant per
    % iteration
    meanX1 = mean(X(find1));
    meanY1 = mean(Y(find1));
    meanZ1 = mean(Z(find1));

    % why random points - why not just subsample by 10 ?
    rng(1)
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    %full_ptCloud1 =  [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
    find1 = find1(p);
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [meanX1, meanY1, meanZ1];
    
    
    %
    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    % compute mean BEFORE random subsampling of points - so constant per
    % iteration
    meanX2 = mean(X(find2));
    meanY2 = mean(Y(find2));
    meanZ2 = mean(Z(find2));

    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);   
    %sub_meanX2 = mean(X(find2));  %% these are a little different..
    %sub_meanY2 = mean(Y(find2));
    %sub_meanZ2 = mean(Z(find2));
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [meanX2, meanY2, meanZ2];
    ptCloud2 = pointCloud(ptCloud2);

    tform = rigid3d(eye(3), [0,0,0]);
    
    % get 100 random numbers
    for i=1:100
        store_rand(i) = rand;
    end
  

    % Example 3. 3D Rigid CPD point-set registration. Full options intialization.
    %  3D face point-set.
    bStoredReg = false;
    if bStoredReg % get from file
        Transform = xforms.store_registration{time_index_index,1};
    else

    sigma2_vect = zeros(100, 1);
    %theta_vect = zeros(100, 3);

    which_rot = 1;
    min_sigma2 = 100;
    counter = 0;

        while ((min_sigma2 > 10) && (counter < maxItr))
            fprintf('.')
            thetaHoldler = {};
            parfor whichrot = 1:gcp().NumWorkers

                newRotation = RandomRotationMatrix(counter+whichrot);
                ptCloud2Loc = ptCloud2.Location*newRotation;


                % registering Y to X
                [Transform,~, sigma2]=cpd_register(ptCloud1,ptCloud2Loc,opt);
                %sigma2_vect(which_rot) = sigma2;
                %store_registration{counter+time_index_index, 1} = Transform;
                transforms{counter+whichrot, 1} = Transform;
                sigma2tests(counter+whichrot) = sigma2;
                init_transforms{counter+whichrot, 1} = newRotation;
                %which_rot = which_rot + 1;

                %figure; hold all; title('Before'); cpd_plot_iter(ptCloud1, ptCloud2);
                % figure; hold all;
                %title([num2str(theta1),';',num2str(theta2),';',num2str(theta3)]);
                
                %cpd_plot_iter(ptCloud1, ptCloud2Loc);
                %figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
            end

            counter = counter + gcp().NumWorkers;%which_rot;
            if counter > 99
                disp('did not find transformation with sigma2 < 10');
            end
            % get the best one we found this loop, 
            [tmp_min_sigma2, min_ind] = nanmin(sigma2tests);
            min_sigma2 = tmp_min_sigma2;
            Transform = transforms{min_ind,1};
            %if verbosemode
            disp(['min sigma, min index ',min_sigma2,min_ind]);
            %end
            store_registration{time_index_index, 1} = Transform;
            store_sigma{time_index_index, 1} = min_sigma2;
            store_init_transform{time_index_index, 1} = init_transforms{min_ind,1};
            
            %close all;
            if verbosemode
                %figure; hold all; title('After registering Y to X.'); cpd_plot_iter(ptCloud1, Transform.Y);
            end
            Rinit = store_init_transform{time_index_index,1};
%             tform = rigid3d(Rinit,[0,0,0]);
%             Y = pctransform(ptCloud2,tform);
%             Y = Y.Location;
%             [M, D]=size(Y);
%             newY = Transform.s * Y * Transform.R'  + repmat(Transform.t',[M 1]);
%             if verbosemode
%                 figure; hold all; title('LB NewY.');cpd_plot_iter(newY, Transform.Y);
%             end
            
            
            % match Transform.Y before Transform to ptCloud2 after Rinit            
%             Rinit = store_init_transform{time_index_index,1};
%             tform = rigid3d(Rinit,[0,0,0]);
%             Y = pctransform(ptCloud2,tform);
%             [M, D] = size(Transform.Y);
%             newY = (Transform.Y - repmat(Transform.t',[M 1]) ) * Transform.R;
%             if verbosemode
%                 figure; hold all; title('LB Y after Rinit to Transform.Y before transform.');cpd_plot_iter(newY, Y.Location);
%             end
             
             % match Transform.Y before Transform and Rinit to ptCloud2
%             Rinit = store_init_transform{time_index_index,1};
%             tform = rigid3d(Rinit',[0,0,0]);
%             [M, D] = size(Transform.Y);
%             newY = (Transform.Y - repmat(Transform.t',[M 1]) ) * Transform.R;
%             newY = newY * Rinit';

%             CombinedR = Transform.R * Rinit';
%             newY = (Transform.Y - repmat(Transform.t',[M 1]) ) * CombinedR;
%             if verbosemode
%                 figure; hold all; title('LB Y before Rinit to Transform Y before all');cpd_plot_iter(newY, ptCloud2.Location);
%             end

            X = ptCloud1;
            [M, D] = size(X);
%             %Xtemp = pctransform(pointCloud(X),tform);
%             newX = ( (1/Transform.s) * ( X  - repmat(Transform.t',[M 1]) ) * Transform.R ) * Rinit' ;
%             tform = rigid3d(Rinit',[0,0,0]);
%             %newXp = newX*R'; % no good
%             %newX = pctransform(pointCloud(newX),tform);
            
            % combined transform
            final_translation = -repmat(Transform.t',[M,1])*Transform.R*Rinit';
            final_rotation = Transform.R*Rinit';
            newX = X*final_rotation + final_translation;
            if verbosemode
                %close all;
                figure; hold all; title('LB New X.'); cpd_plot_iter(newX, ptCloud2.Location);
            end
            %FullTransform = rigid3d(final_rotation,final_translation(1,:));
            s.Rotation = final_rotation;
            s.Translation = final_translation;
            s.Centroids1 = [meanX1,meanY1,meanZ1];
            s.Centroids2 = [meanX2,meanY2,meanZ2];
            s.NumberTrials = counter;
            s.minSigma = min_sigma2;
            store_registration{time_index_index, 1} = s;

        end
    end
    
    %% compute the matches 
    [iou_matrix, M, corresponding_ious_for_matches, ...
            cell_labels_I_care_about1, cell_labels_I_care_about2, ...
            center_point_for_each_label1, center_point_for_each_label2, ...
            match_based_on_nearest_neighbors, ~, ~, ...
            alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,ptCloud1,...
            combined_image1,combined_image2,find1,find2);
     store_matches{time_index_index, 1} = M;
     store_iou_table{time_index_index, 1} = iou_matrix;
     
    %% make the graph..
        
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
    % visualization for checking if everything is correct - 3d plot of
    % edges and nodes
    if verbosemode
        hold all; plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0);
        figure; h1 = plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name);
    end
    %disp(time_index);
   
    % loop through all matches
    % LB: M is nMatches x 2 (label at time t, label at time t+1)  

    for i = 1:length(M)
        find_non_zero_please = find(iou_matrix(M(i,1),:));
        if (length(find_non_zero_please) > 1)  % more than one iou match
            if verbosemode
                figure; hold all; cpd_plot_iter(ptCloud1, Transform.Y);
                hold all; plot(alpha_shape_for_each_label1{M(i,1),1},'FaceColor','red','FaceAlpha',0.5);
                for j = 1:length(find_non_zero_please)
                    if (find_non_zero_please(j) == M(i,2)) % this is the best match?
                        hold all; plot(alpha_shape_for_each_label2{M(i,2),1},'FaceColor','green','FaceAlpha',0.5);
                    else
                        hold all; plot(alpha_shape_for_each_label2{find_non_zero_please(j),1},'FaceColor','black','FaceAlpha',0.5);
                    end
                end                
            end
            
            title([num2str(corresponding_ious_for_matches(i)),';',num2str(i)]);
        end
    end
    
    
    %% compute size of each nucleus
    % if reg failed - remove two smallest
    if min_sigma2 > 10
        %% find nuclei with no match
        for iImage = 1:2
            iNoInd = 1;
            NoM = [];
            if iImage == 1
                nlabels = size( unique(combined_image1), 1) - 1;
            else
                 nlabels = size( unique(combined_image2), 1) - 1;
            end
            for itest = 1: nlabels
                bFound = false;
                for imatch= 1:length(M)
                    if (M(imatch,iImage) == itest)
                        bFound = true;
                    end
                end
                if ~bFound
                    NoM(iNoInd) = itest;
                    iNoInd = iNoInd + 1;
                end
            end
            %disp(length(NoM));
            %disp(NoM);
            % remove these nuclei from image
            for ind =1:length(NoM)
                iLabel = NoM(ind);
                if iImage == 1
                    nucleus_pts = find(combined_image1(:)==iLabel); 
                    combined_image1(nucleus_pts) = 0;
                else
                    nucleus_pts = find(combined_image2(:)==iLabel); 
                    combined_image2(nucleus_pts) = 0;
                end
            end
        end   
     
        bRemoveSmall = false;
        if bRemoveSmall
            % get nuclei labels
            nlabels = size( unique(combined_image2), 1) - 1;
            vol = [];
            for ilabel = 1:nlabels
                 nucleus_pts = find(combined_image2(:)==ilabel);  % this is the indices into combined_image1 to get indices into (X,Y,Z) to the full set of point
                vol(ilabel) = length(nucleus_pts);
            end
            if verbosemode
                figure;
                plot(1:nlabels,vol,'LineWidth',5);
                xlabel('Nucleus Label ID');
                ylabel('Nucleus Volume');
            end
            % order by size
            [sort_vol,sort_idx] = sort(vol);
             % remove two smallest (if smaller than 1000 ?? depends on time)
            smallest_id = sort_idx(1);
            nucleus_pts = find(combined_image2(:)==smallest_id); 
            combined_image2(nucleus_pts) = 0;
            second_smallest_id = sort_idx(2);
            nucleus_pts = find(combined_image2(:)==second_smallest_id); 
            combined_image2(nucleus_pts) = 0;
        end
        if ~bRerun
            bRerun = true;
        else
            bRerun = false;
            time_index_index = time_index_index + 1;
        end
    else
        bRerun = false;
        time_index_index = time_index_index + 1;
    end
    close all;
    disp(' ')
    toc
end

% Save vector of transformations...
save(RegistrationFileName, 'store_registration');
%save(InitFileName,'store_init_transform');
save(strcat(data_path,'matches.mat','store_matches'));
save(strcat(data_path,'iou_table.mat','store_iou_table'));
save(strcat(data_path,'graph.mat','G_based_on_nn'));


%% example to add node and edge after split
%G_based_on_nn = G_based_on_nn.addnode('078_016');
%G_based_on_nn = G_based_on_nn.addedge('077_010','078_016')
%end
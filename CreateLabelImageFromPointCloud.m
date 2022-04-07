
% given transformed point cloud - make label image
% needs original label number for each point
% also too downsampled - so only sparse set of points..
addpath('/Users/lbrown/Downloads/bfmatlab');
T = load('transform_labels_pCloud.mat');
for itime = 131:131
    s = T.store_registration{itime};
    pts = s.Y;
    npts = size(pts,1);
    label_img = zeros(160,160,160,'uint8');
    for ipt = 1:npts
        p = pts(ipt,:);
        label_img(round(p(1)+80),round(p(2)+80),round(p(3)+80)) = 255;
    end

    %plane = zeros(64, 64, 1, 2, 2, 'uint8');
    bfsave(label_img, 'TransformedLabels.ome.tiff');


end

        
    
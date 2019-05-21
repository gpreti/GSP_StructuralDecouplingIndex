%% This function uses routines formerly implemented in MIP:Lab (Jonas
% Richiardi) to plot brain plots
%
% Inputs:
% - CM is a n_ROI x n_ROI matrix denoting the strength of connections (if
% wished in the plot). One can enter a matrix full of zeros if connections
% should not be displaye
% - CC is a vector representing the size of the nodes in the plot
% - CodeBook contains the information about where to locate the nodes of
% the atlas, and how to label them
% - T_conn is the treshold above which connections will be represented. For
% example, if we have a FC matrix and T_conn=0.5, only the strongest
% positive connections will be displayed
% - Factor_SphereSize rescales the size of the nodes; recommanded value is
% the max of the CC vector or something similar
% - Same for Factor_Col, but controling the color of the nodes
% - Exp_Sphere controls the linearity of the sphere sizes (quadratically
% increasing in the display if Exp_Sphere=2)
% - View is an integer from 1 to 3 depicting the type of slices to show
% - Colormap_nodes is the matlab colormap to use to represent nodes (e.g.
% 'hot' or 'autumn'
% - Colormap_edges is the matlab colormap to use to represent edges
% - Gamma is a factor controling the linear of edge thickness when plotting
% connections
% - LinearWeight is a rescaling factor for connections
% - CA can be set to [] to have dynamic ranges for the connection colormap;
% if set to something specific, that range is used for the plotting
function [] = PlotBrainGraph(CM,CC,CC2,CodeBook,T_conn,Factor_SphereSize,...
    Factor_Col,Exp_Sphere,View,Colormap_nodes,Colormap_edges,Gamma,...
    LinearWeight,CA)

    % edges
    pp.doPlotConnections=1;
    pp.pruneTinyConnections=1;
    pp.pruningWeight=T_conn;
    pp.linearWeight=LinearWeight;
    pp.gamma=Gamma;
    pp.doColourCode=true;
    pp.cmapchoiceE=Colormap_edges;

    % nodes
    pp.doPlotDegreeSpheres=true;%**
    pp.myDSplotT=0;      % this * 1/N  is the threshold
    pp.showLabels=0;
    pp.cmapchoiceV=Colormap_nodes; % gray
    
    % Controls the exponent for the display of the spheres
    pp.e_matD_to_sphS=Exp_Sphere;
    
    
    pp.f_matD_to_sphS=1;

    % brain outline
    pp.bmesh_faceAlpha=0.2; % transparency

    % node size and color
    pp.FctForColSph=CC2/Factor_Col;
    pp.FctForSizSph=CC/Factor_SphereSize;
    
    % Edge transparency
    pp.doAlpha = true;
    pp.alphagamma = 0.05;

    figure
    
    if ~isempty(CA)
    
      [~]=show_cm_extended(CM,CC,CC2,CodeBook,3,'gamma',pp.gamma,'FctForColSph',pp.FctForColSph,...
        'linearWeight',pp.linearWeight,'showLabels',pp.showLabels,...
        'createFigure',false,'doColourCode',pp.doColourCode,'doAlpha',pp.doAlpha,'alphagamma',pp.alphagamma,...
        'cmapchoiceV',pp.cmapchoiceV,'cmapchoiceE',pp.cmapchoiceE,...
        'doPlotConnections',pp.doPlotConnections,'doPlotDegreeSpheres',pp.doPlotDegreeSpheres,...
        'bmesh_show',1,'bmesh_hemi','lr','bmesh_useConvexHull',false,...
        'bmesh_convexHullOutlinePlanes','','FctForSizSph',pp.FctForSizSph,...
        'bmesh_faceAlpha',pp.bmesh_faceAlpha,'pruneTinyConnections',pp.pruneTinyConnections,...
        'pruningWeight',pp.pruningWeight,'e_matD_to_sphS',pp.e_matD_to_sphS,'f_matD_to_sphS', pp.f_matD_to_sphS,...
        'plotDegreeSpheresOnlyAboveThreshold',pp.myDSplotT,...
        'showTopNlabels',0,'txtPosOffset',5*[1 1 1],'ColAxis',CA,...
        'showTopNlabelsBy','vertexW','useShortNames',0);
    
    else
        [~]=show_cm_extended(CM,CC,CC2,CodeBook.full,3,'gamma',pp.gamma,'FctForColSph',pp.FctForColSph,...
        'linearWeight',pp.linearWeight,'showLabels',pp.showLabels,...
        'createFigure',false,'doColourCode',pp.doColourCode,'doAlpha',pp.doAlpha,'alphagamma',pp.alphagamma,...
        'cmapchoiceV',pp.cmapchoiceV,'cmapchoiceE',pp.cmapchoiceE,...
        'doPlotConnections',pp.doPlotConnections,'doPlotDegreeSpheres',pp.doPlotDegreeSpheres,...
        'bmesh_show',1,'bmesh_hemi','lr','bmesh_useConvexHull',false,...
        'bmesh_convexHullOutlinePlanes','','FctForSizSph',pp.FctForSizSph,...
        'bmesh_faceAlpha',pp.bmesh_faceAlpha,'pruneTinyConnections',pp.pruneTinyConnections,...
        'pruningWeight',pp.pruningWeight,'e_matD_to_sphS',pp.e_matD_to_sphS,'f_matD_to_sphS', pp.f_matD_to_sphS,...
        'plotDegreeSpheresOnlyAboveThreshold',pp.myDSplotT,...
        'showTopNlabels',0,'txtPosOffset',5*[1 1 1],...
        'showTopNlabelsBy','vertexW','useShortNames',0);
    end


    %colorbar;
    axis off
    axis tight;

    if View==1, az=0; el=90; cam='right'; vview='ax';
    elseif View==2, az=90; el=0; cam='right'; vview='sag';
    else az=0; el=0; cam='left'; vview='cor';
    end
    view(az,el);
    camlight(cam);

    disp('done');

end
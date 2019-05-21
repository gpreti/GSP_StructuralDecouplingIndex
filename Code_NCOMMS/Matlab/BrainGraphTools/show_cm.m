function [myH,varargout]=show_cm(CM,cod,mode,varargin)
% SHOW_CM show a brain connectivity graph with many customisable options
%
% EXAMPLE:
%   show_cm(myCM,myCB,3);
%
% INPUT:
%  CM: nRegions x nRegions adjacency matrix. Entries indicate graph edge
%       weights.
%  cod: a struct for a codebook with fields
%       .center: 1 x nRegions cell array of 3x1 vectors - MNI coordinates of each region
%       .id: nRegions x 1 array of integers - code values corresponding to
%           atlas values (identify each region)
%       .name: nRegions x 1 cell array of strings - region names
%       .sname: nRegions x 1 cell array of strings - SHORT region names
%       .num: a scalar containing the number of regions defined in the
%               codebook
%  mode: scalar -
%           1 displays adjacency matrix
%           2 display 2D projections of graph
%           3 shows graph in 3D space
% 
% ** NAME/VALUE PAIRS (incomplete list)
%
% * OPTIONS TO SELECT WHAT TO PLOT
%
%  'doPlotConnections': boolean - whether to plot graph edges (true) or
%        not (default: true)
%  'doPlotDegreeSpheres': boolean - whether to plot spheres on vertex
%       centers representing the degree of the node or a given function (default false)
%  'plotDegreeSpheresOnlyAboveThreshold': factor by which degree sphere
%       must be larger than uniform weight (1/nVertices) in order to be
%       plotted. 0 - plot all degree spheres (default: 0) / if used together
%       with FctForColSph, sphere plotted for abs(fct) > threshold
%  'nodesXonly': scalar or row vector - show only connections of the specified node(s)
%       (default: [])
%  'pruneTinyConnections': boolean - don't plot connections with very small weight
%  'pruningWeight': scalar - weight below which connections should not be
%   plotted
%
%
% * OPTIONS REGARDING TEXT LABELS
%
%  'showLabels': logical - show region labels or not (default: true)
%  'showWeights': logical - show edge weights or not (defalt: false)
%  'myFontSize': scalar between 6 and 50 - size of font for region labels
%        (default: 12)
%  'showTopNlabels': scalar or empty - if scalar, choice of how many
%        region labels to display - these will be for top-weighted edges. If
%        empty, label display is either all (showLabels==true) or none.
%        (default: [])
%   'showTopNlabelsBy': whether to take into account edge weights 'edgeW' or vertex
%   weights 'vertexW'. Default: 'edgeW';
%  'useShortNames': logical - whether to use short names for labels or not
%       (default: false)
%  'txtPosOffset': 1x3 vector - additional label display offset in x,y,z
%       with respect to the value in center{r}+posOffset (used to make sure
%       labels can be placed above nodes and edges in certain projections)
%       (default: 0 0 0)
%  'flipLabelsRL': replace 'L' with 'R' at end of region name and
%       vice-versa (default: false)
%
%
% * OPTIONS TO CONTROL HOW OBJECTS ARE SCALED
%
%  'linearWeight': scalar - linear scaling for edge weights display
%       (default: 2)
%  'normMax': scalar or empty - if scalar, maximum edge weights with
%        respect to which all others are scaled. If empty, maximum is
%        computed from data in CM (default: [])
%  'gamma': scalar between 0 and 1 - nonlinear scaling for display of edge
%       weights (towards 0: small edge weights are inflated, towards 1: small
%       edge weights are left small) (default: 0.5)
%  'f_matD_to_sphS':  conversion factor between node degrees and
%       sphere size for degree plot (default: 1/15)
%   'e_matD_to_sphS': exponent for conversion between node degrees and
%   sphere sizes (default: 2);
%  'FctForSizSph': scalar or nVerticesx1, uses same value for all nodes if scalar, 
%       (default: empty)
%  'computeDegreesOnFullMat': boolean - compute the weights and degrees
%   on full (not sparsified through e.g. nodesXonly) adjacency matrix.
%   Default true
%
%
% * OPTIONS TO CONTROL COLOURS OF VARIOUS ELEMENTS
%
%  'doAlpha': logical - if true, edges with smaller weights are displayed
%       more transparent (default: true)
%  'doColourCode': logical - if true, edges weight is represented with
%       colour coding in addition to line thickness (default: true)
%  'cmapchoiceE': string - colormap to use for mapping edge weights to
%       colours. Use standard colormap weights (default: 'autumn')
%  'cmapchoiceV': string - colormap to use for mapping Vertex weights to
%       colours. Use standard colormap weights (default: 'autumn')
%  'FctForColSph': nVerticesx1 double, colors spheres on vertices based on
%       input instead of degree of sphere (default: [])
%  'CdataVal': scalar above 0 - value to put in faces of graph edges to
%       map into a colour map (used if doColourCode==false, useful not to clash
%       with brain mesh values or other objects) (default: 2)
%       if = 0, connections colored in gray shade (default: black)
%       (works only with useSplitcMap = false)
%  'GrayShade': scalar [0 1] : colors connections in gray shade (0: black,
%        1: white) if 'CdataVal'=0
%  'SphBlack': boolean, if true plots spheres below threshold as gray shade
%       (c.f. 'GrayShade') instead of omitting them (default: false)
%  'BGColour': ColorSpec - background colour for figure (default: 1 1 1)
%  'ColAxis': [clim cmax], manually set limits of colorbar (default: [])
%  'useSplitCMap': logical - rescale values with respect to splitID so
%       that different colourmaps can be applied to different subgraphs
%       ! Experimental code - use freezecolors from Matlab FEX instead
%       (default: false)
%  'splitID': scalar in {1,2} - subgraph identifier, used with
%       useSplitCMap (default: [])
%       ! Experimental code
%  'alphagamma': scalar between 0 and 1 - further non-linear scaling
%       applied to edge transparency (allows finer control) (default: .1)
% 
% 
% * OPTIONS REGARDING 3D CORTICAL MESH
%
%  'bmesh_show': logical - show brain mesh or not (default: false)
%  'bmesh_hemi': string in 'l','r','lr' - select which hemisphere to
%       display (default: 'lr')
%  'bmesh_faceColour': ColorSpec - colour for the faces of the brain mesh
%       triangles (default: .8 .8 .8)
%  'bmesh_faceAlpha': scalar between 0 and 1 - alpha transparency value to
%       use for brain mesh triangles (0: fully transparent, 1: fully
%       opaque) (default: .4)
%  'bmesh_edgeColour': ColorSpec - colour for the edges of the brain mesh
%       triangles (default: .5 .5 .5)
%  'bmesh_edgeAlpha': scalar between 0 and 1 - alpha transparency value to
%       use for brain mesh triangle edges (0: fully transparent, 1: fully
%       opaque) (default: .05)
%  'bmesh_useConvexHull': boolean - whether to compute and display the
%        convex hull rather than the full mesh (default: false)
%  'bmesh_convexHullOutlinePlanes': string ('', 'a','s','c'): if not empty,
%        draw an outline of the convex hull on the corresponding plane (axial,
%        sagittal, or coronal). Needs the bmesh_useConvexHull option to work.
%        (default: '')
%  'bmesh_convexHullOutlinePlanesLineWidth':  width of convex hull outline
%        plot (default: 1)
%   'bmesh_reducePatchVal': scalar - fraction of patches to retain for
%   drawing brain meshes - increasing this from default of 0.015 may crash
%   Matlab.
%
%
% * GENERAL AND MISC OPTIONS
%
% 'createFigure: logical - whether a new figure should be created or not
%       (default: true)
% 'posOffset': 3x1 vector - display offset in x,y,z with respect to the
%   value in center{r} (default: 0 0 0)
% 'Title': add a title to the plot (default: [])
%
%
% OUTPUT:
%   myH: a cell array of graphics handles.
%       myH{1}: scalar, contains handle to the figure if createFigure is
%           true, empty otherwise
%       myH{2}: vector, contains 1 or 2 handles for the patch objects
%           representing the hemispheric meshes or convex hull thereof
%           if bmesh_hemi is set, otherwise is empty
%       myH{3}: vector, contains 1 to 3 handles for the plot3 line objects
%           representing the a,s,c projections of the convex hull of the
%           hemispheric meshes, if bmesh_convexHullOutlinePlanes is set,
%           otherwise is empty
%       myH{4}: vector, contains handles to all surface objects (degree
%           spheres and connections)
%       myH{5}: vector, contains handles to all text objects
%   varargout(1):
%       wCM: a reweighted and thresholded CM corresponding to
%           what show_CM displays. the diagonal corresponds to region
%           sphere sizes.
%
% VERSION HISTORY:
% v1.0 (date unknown) Dimitri Van de Ville
% - initial release
% v2.0 oct-nov 2009 Jonas Richiardi
% - allow for adjacency matrix plotting (zero values = no edge)
% - nicer color mapping
% - 3D view with patches and transparency control
% - all weights rescaled to max 1
% - optional node label and weights display
% - gamma exponent for non-linearity selectable
% - alpha transparency on/off
% - single-colour or colour mapped to link strength
% - returns figure handles
% - many more...
% v2.0.1 2009 Nov 10 Jonas Richiardi
% - code cleanup for combined mesh + connection view
% v2.1 2010 April 21 Jonas Richiardi
% - convex hull computation and plotting, and outline of convex hull.
% v3.0 jul 2010 Nora Leonardi
% - print only connections of 1 node
% v3.1 jul 2010 Jonas Richiardi
% - cleaner handle return structure. Separate handles for figure, hemispheric
%   meshes, hull outline projectsion, surface-plotted objects, and text, to
%   allow the use of separate renderers (painters / opengl) when exporting
%   to eps. This allows the generate of high-quality text and surface
%   objects for publications, and keeps editorial offices and compositors
%   happy (see jHighQRenderingToEPS.m).
% - other minor cleanups and improvements
% v3.1.1 jul 2010 Jonas Richiardi
% - bugfix on 3.1.0: handles of only half coronal and axial outlines were
%   available
% v3.1.2 Jul 2010 Jonas Richiardi
% - don't floor and display very small weight connections, instead don't
% display them at all
% v 3.2 Mar 2011 JR
% - minor bugfixes
% v3.2.1 Nov 2011 JR
% - partial merge with Nora's branch - FctForSizSph
% v3.2.2 Nov 2011 NL
% - FctForSizSph can now be used to control sphere size in parallel
% with FctForColSph controlling the colour
% v3.2.3 Mar 2013 JR
% - fix logical condition code for nodesXonly
% v3.3 April 2013 JR
% - reorganised options, now sorted alphabetically
% - harmonised capitalisation between incoming options and internal names
% - improved documentation, but still work in progress
% v3.4 Aug 2013 NL
% - degrees calculated using abs(CM)
% - removed Cdataval and Grayshade warning as called with default inputs
% - changed pp.FctForSizSph: input is scalar or nVerticesx1

% showWeights off by default, text() is a slow function
[pp.alphagamma,pp.BGColour,pp.bmesh_convexHullOutlinePlanes,...
    pp.bmesh_convexHullOutlinePlanesLineWidth,pp.bmesh_edgeAlpha,...
    pp.bmesh_edgeColour,pp.bmesh_faceAlpha,pp.bmesh_faceColour,...
    pp.bmesh_hemi,pp.bmesh_reducePatchVal, pp.bmesh_show,...
    pp.bmesh_useConvexHull,pp.CdataVal,pp.cmapchoiceE,pp.cmapchoiceV,...
    pp.ColAxis,pp.computeDegreesOnFullMat,pp.createFigure,pp.doAlpha, ...
    pp.doColourCode,pp.doPlotConnections,pp.doPlotDegreeSpheres,...
    pp.e_matD_to_sphS, pp.f_matD_to_sphS,pp.FctForColSph,...
    pp.FctForSizSph,pp.flipLabelsRL,pp.gamma,pp.GrayShade,...
    pp.linearWeight,pp.myFontSize,pp.nodesXonly,pp.normMax,...
    pp.plotDegreeSpheresOnlyAboveThreshold,pp.posOffset,...
    pp.pruneTinyConnections,pp.pruningWeight,pp.showLabels,...
    pp.showTopNlabels,pp.showTopNlabelsBy,pp.showWeights,pp.SphBlack,...
    pp.splitID,pp.Title,pp.txtPosOffset,pp.useShortNames,...
    pp.useSplitCMap]=...
    process_options(varargin,'alphagamma',0.1,'BGColour',[1 1 1],...
    'bmesh_convexHullOutlinePlanes','',...
    'bmesh_convexHullOutlinePlanesLineWidth',1,'bmesh_edgeAlpha',0.05,...
    'bmesh_edgeColour',[0.5 0.5 0.5],'bmesh_faceAlpha',0.05,...
    'bmesh_faceColour',[0.8 0.8 0.8],'bmesh_hemi','lr',...
    'bmesh_reducePatchVal',0.015,'bmesh_show',false,...
    'bmesh_useConvexHull',false,'CdataVal',2,'cmapchoiceE','autumn',...
    'cmapchoiceV','autumn','ColAxis',[],'computeDegreesOnFullMat',true,...
    'createFigure',true,'doAlpha',true,'doColourCode',true,...
    'doPlotConnections',true,'doPlotDegreeSpheres',false,...
    'e_matD_to_sphS',2,'f_matD_to_sphS', 1/15,'FctForColSph',[],...
    'FctForSizSph',[],'flipLabelsRL',false,'gamma',0.5,...
    'GrayShade',0,'linearWeight',2,'myFontSize',12,'nodesXonly',[],...
    'normMax',[],'plotDegreeSpheresOnlyAboveThreshold',0,...
    'posOffset',[0 0 0]','pruneTinyConnections',false,...
    'pruningWeight',0.01,'showLabels',true,...
    'showTopNlabels',[],'showTopNlabelsBy','edgeW','showWeights',false,...
    'SphBlack',false,'splitID',[],'Title',[],'txtPosOffset',[0 0 0]',...
    'useShortNames',false,'useSplitCMap',false);

%% CONSTANTS
% handle cell array identification constants in myH
H_F=1;  % handle to the figure
H_H=2;  % handle to hemispheric meshes
H_O=3;  % handle to brain outlines
H_S=4;  % handle to surface objects
H_T=5;  % handle to text objects


%% sanity check
% check codebook
allCBfields=fieldnames(cod);
needCBfields={'id','num','name','sname','center'};
for f=1:numel(needCBfields)
    if all(cellfun('isempty',strfind(allCBfields,needCBfields{f})))
        error(['Codebook is lacking field: ' needCBfields{f}]);
    end
end

% check params are within permissible ranges
if (pp.gamma <0) || (pp.gamma > 1)
    error('gamma must be between 0 and 1');
end
if (pp.linearWeight<=0)
    error('linearWeight must be above 0');
end
if (pp.bmesh_faceAlpha<0) || (pp.bmesh_faceAlpha>1)
    error('bmesh_faceAlpha must be between 0 and 1');
end
if (pp.myFontSize<7) || (pp.myFontSize>50)
    error('bmesh_faceAlpha must be between 7 and 50');
end
if ~isempty(pp.showTopNlabels) && pp.showLabels && ...
        ((pp.showTopNlabels <1) || (pp.showTopNlabels > cod.num))
    error('showTopNlabels must be between 1 and cod.num');
end
if pp.CdataVal<0
    error('CdataVal must be greater than 0');
end

if ~isempty(pp.bmesh_convexHullOutlinePlanes)
    if pp.bmesh_useConvexHull==false || pp.bmesh_show==false
        error('bmesh_useConvexHull and bmesh_show must be true to use bmesh_convexHullOutlinePlanes');
    end
    for p=1:numel(pp.bmesh_convexHullOutlinePlanes)
        if all(cellfun('isempty',strfind({'a','s','c'},pp.bmesh_convexHullOutlinePlanes(p))))
            error('bmesh_convexHullOutlinePlanes must be empty, a, s, or c');
        end
    end
end

% if pp.FctForSizSph && isempty(pp.FctForColSph)
%     error('FctForColSph cannot be empty if FctForSizSph is true');
% end

% check nodes to plot exist
if ~isempty(pp.nodesXonly)
    if any(pp.nodesXonly > cod.num | pp.nodesXonly==0)
        error('elements of nodesXonly must be valid indices of CM (1 to %d)', cod.num);
    end
end

% check splitmap
if (pp.useSplitCMap==true)
    if isempty(pp.splitID)
        error('Must specify splitID when using useSplitCMap');
    elseif ~ismember(pp.splitID,[1 2]);
        error('splitID must be 1 or 2');
    end
end
%if (pp.alphagamma <0) || (pp.alphagamma > 1)
%    error('alphagamma must be between 0 and 1');
%end
% check offsets
if ~all(size(pp.posOffset)==[3 1])
    if ~all(size(pp.posOffset)==[1 3])
        error('Please supply column vectors for posOffset');
    else
        pp.posOffset=pp.posOffset';
    end
end
if ~all(size(pp.txtPosOffset)==[3 1])
    if ~all(size(pp.txtPosOffset)==[1 3])
        error('Please supply column vectors for txtPosOffset');
    else
        pp.txtPosOffset=pp.txtPosOffset';
    end
end

% check CDataVal=0 if GrayShade not empty
% if pp.CdataVal~=0 && ~isempty(pp.GrayShade)
%     warning('GrayShade not used since CdataVal > 0');
% end

if pp.SphBlack==true && isempty(pp.plotDegreeSpheresOnlyAboveThreshold)
    error('Please specifiy a threshold below which spheres are plotted in gray shade: plotDegreeSpheresOnlyAboveThreshold');
end

% compute complementary colour for text
pp.compColour=abs([1 1 1]-pp.BGColour);

% replace long names by short ones if needed
if (pp.useShortNames==true)
    cod.name=cod.sname;
end

if ~isempty(pp.FctForSizSph) && length(pp.FctForSizSph)==1
    pp.FctForSizSph=pp.FctForSizSph*ones(size(CM,1),1);
end

%% cut adjacency matrix if only connetions of 1 node to be shown

if ~isempty(pp.nodesXonly)
    A=zeros(cod.num);
    Nnodes=length(pp.nodesXonly);
    for k=1:Nnodes
        colx=CM(:,pp.nodesXonly(k));
        rowx=CM(pp.nodesXonly(k),:); % could also reuse colx as CM sym
        A(:,pp.nodesXonly(k))=colx;
        A(pp.nodesXonly(k),:)=rowx;
    end
    CMold = CM; % keep old matrix to calc node degrees
    CM = A;
else
    CMold = CM;
end

%% main selector
% display adjacency matrix
if mode==1,
    
    if (pp.createFigure==true)
        figure;
        caxis([-1 1]);
    end
    
    imagesc(CM);
    h=get(gcf,'CurrentAxes');
    
    if (pp.showLabels==true)
        set(h,'Ytick',[1:length(CM)]);
        
        ml=0;
        for iter=1:cod.num;
            ml=max(ml,length(cod.name{iter}));
        end;
        
        dispOffset=5;
        ml=ml+dispOffset; % save space for staggering
        % pre-fill str with blanks
        str=repmat([' '],[cod.num ml]);
        
        for iter=1:cod.num;
            if mod(iter,2)==0
                offset=dispOffset;
            else
                offset=1;
            end
            str(iter,offset:offset+length(cod.name{iter})-1)=cod.name{iter};
        end;
        
        set(h,'YTickLabel',str,'FontSize',8);
    end
    
    caxis([-1 1]);
    axis square; axis tight;
    
    % plot connections in  2D or 3D
else
    % normalise by max
    if isempty(pp.normMax)
        maxVal=max(max(CM));
        maxValOld=max(max(CMold));
    else
        maxVal=pp.normMax;
        maxValOld=pp.normMax;
    end
    CM=CM./maxVal;
    CMold=CMold./maxValOld;
    % get column vectors
    for r=1:cod.num
        cod.center{r}=cod.center{r}(:);
    end
    
    % 2D brain view
    if mode==2,
        warning('returned handle vector is not compatible with show_cm v3.1');
        myH=zeros(3,1);
        % display
        myH(1)=brain_plot(CM,cod,[1 3],pp);
        title('coronal'); xlabel('left-right'); ylabel('ventral-dorsal');
        myH(2)=brain_plot(CM,cod,[1 2],pp);
        title('axial'); xlabel('left-right'); ylabel('posterior-anterior');
        myH(3)=brain_plot(CM,cod,[2 3],pp);
        title('sagittal'); xlabel('posterior-anterior'); ylabel('ventral-dorsal');
        
        % 3D brain view
        
    elseif mode==3
        if (pp.bmesh_show==true)
            nHemis=numel(pp.bmesh_hemi);
        else
            nHemis=0;
        end
        myH=cell(5,1);
        if nHemis==0
            myH{H_H}=[];
        else
            myH{H_H}=zeros(nHemis,1); % hemispheric meshes or hulls
        end
        if (pp.createFigure==true)
            myH{H_F}=figure;
        else
            myH{H_F}=[];
        end
        hold on;
        caxis([0 1]);
        
        % show brain mesh
        if (pp.bmesh_show==true)
            myH{H_H}=show_FSmesh(pp.bmesh_hemi,...
                'white','createFigure',false,'FaceColor',pp.bmesh_faceColour,...
                'FaceAlpha',pp.bmesh_faceAlpha,'EdgeColor',pp.bmesh_edgeColour,...
                'EdgeAlpha',pp.bmesh_edgeAlpha,'useConvexHull',pp.bmesh_useConvexHull, ...
                'reducePatchVal',pp.bmesh_reducePatchVal);
            % setup mesh colours
            for hi=1:nHemis
                set(myH{H_H}(hi),'FaceColor',pp.bmesh_faceColour,...
                    'FaceAlpha',pp.bmesh_faceAlpha,'EdgeColor',pp.bmesh_edgeColour,...
                    'EdgeAlpha',pp.bmesh_edgeAlpha);
            end
        end
        
        % plot graph
        [myH{H_S},myH{H_T},wCM]=brain_plot3d(CM,cod,pp,CMold);
        varargout(1)={wCM};
        
        myH{H_T}(end+1)=xlabel('left-right','Color',pp.compColour);
        myH{H_T}(end+1)=ylabel('posterior-anterior','Color',pp.compColour);
        myH{H_T}(end+1)=zlabel('ventral-dorsal','Color',pp.compColour);
        if ~isempty(pp.Title), myH{H_T}(end+1)=title(pp.Title); end
        
        % plot hull outline planes if needed
        if (pp.bmesh_useConvexHull==true) && ~isempty(pp.bmesh_convexHullOutlinePlanes)
            myH{H_O}=[]; % zeros(numel(pp.bmesh_convexHullOutlinePlanes)*nHemis,1);
            Xs=cell(nHemis,1); Ys=cell(nHemis,1); Zs=cell(nHemis,1);
            kas=cell(nHemis,1);kss=cell(nHemis,1);kcs=cell(nHemis,1);
            for hi=1:nHemis
                % get all 3D vertices locations in convex hull
                xd = get(myH{H_H}(hi),'XData');
                yd = get(myH{H_H}(hi),'YData');
                zd = get(myH{H_H}(hi),'ZData');
                Xs{hi}=[xd(1,:) xd(2,:) xd(3,:)];
                Ys{hi}=[yd(1,:) yd(2,:) yd(3,:)];
                Zs{hi}=[zd(1,:) zd(2,:) zd(3,:)];
                
                % now we can compute the convex hull of the projections
                kas{hi}=convhull(Xs{hi},Ys{hi}); % convex hull of the axial projection
                kss{hi}=convhull(Ys{hi},Zs{hi}); % sagittal
                kcs{hi}=convhull(Xs{hi},Zs{hi}); % coronal
                
                % we could also create slice plane with values - see "exploring volumes with slice
                %planes"
            end
            %figure; scatter3(Xs{1},Ys{1},Zs{1},'k'); hold on;
            %scatter3(Xs{2},Ys{2},Zs{2},'k'); axis equal;
            sagHasBeenPlotted=false;
            for hi=1:nHemis
                % delete old 3D convex hull graphics object and set its
                % handle to nana
                delete(myH{H_H}(hi));
                myH{H_H}(hi)=nan;
                hold on;
                
                % and now plot this stuff to check
                for p=1:numel(pp.bmesh_convexHullOutlinePlanes)
                    switch pp.bmesh_convexHullOutlinePlanes(p)
                        case 'a'
                            myH{H_O}(end+1)=plot3(Xs{hi}(kas{hi}),Ys{hi}(kas{hi}),zeros(numel(kas{hi}),1),...
                                'k--','LineWidth',pp.bmesh_convexHullOutlinePlanesLineWidth); % plot axial
                        case 's'
                            if ~sagHasBeenPlotted % avoid plotting twice - hemispheres one behind the other
                                myH{H_O}(end+1)=plot3(zeros(numel(kss{hi}),1),Ys{hi}(kss{hi}),Zs{hi}(kss{hi}),'k--',...
                                    'LineWidth',pp.bmesh_convexHullOutlinePlanesLineWidth); % plot sagittal
                                sagHasBeenPlotted=true;
                            end
                        case 'c'
                            myH{H_O}(end+1)=plot3(Xs{hi}(kcs{hi}),zeros(numel(kcs{hi}),1),Zs{hi}(kcs{hi}),'k--',...
                                'LineWidth',pp.bmesh_convexHullOutlinePlanesLineWidth); % plot coronal
                        otherwise
                            error('bmesh_convexHullOutlinePlanes must be a,s,c. There is no way you should be here!');
                    end
                end % all convex hull outline planes
            end % all hemispheres
        end % plot convex hull outline plane
    else
        error('Unknown visualisation mode');
    end
end




function myH=brain_plot(CM,cod,coor,pp) % 2D

warning('show_cm:ancientCode','Ancient code. Check myH return format and plotting details.');

myH=figure;hold on;
for iter=1:cod.num,
    scatter(cod.center{iter}(coor(1)),cod.center{iter}(coor(2)),'o','MarkerEdgeColor',pp.compColour);
    if (pp.showLabels==true)
        text(cod.center{iter}(coor(1)),cod.center{iter}(coor(2)),...
            cod.name{iter},'Interpreter','none','HorizontalAlignment',...
            'center','Color',pp.compColour);
    end
    
end;
for iter=1:cod.num,
    for iter2=iter+1:cod.num,
        if ~isnan(CM(iter,iter2)) && CM(iter,iter2)~=0
            % compute display params
            myW=abs(CM(iter,iter2));
            edgeWeight=pp.linearWeight*myW^pp.gamma;
            %colVal=exp(2*(1-myEdges(e).Weight));
            colVal=1-myW;
            if colVal < 0.2, colVal=0.2; end
            if colVal > 0.9, colVal=0.9; end
            h=plot([cod.center{iter}(coor(1)) cod.center{iter2}(coor(1))],[cod.center{iter}(coor(2)) cod.center{iter2}(coor(2))],'b-');
            set(h,'LineWidth',edgeWeight);
            set(h,'Color',[colVal colVal colVal]);
            if (pp.showWeights==true)
                midPt=mean([cod.center{iter}([coor(1) coor(2)]) ...
                    cod.center{iter2}([coor(1) coor(2)])],2);
                text(midPt(1),midPt(2),num2str(edgeWeight,'%1.2f'));
            end
        end;
    end;
end;
axis equal;
hold off;





function [myHs,myHt,wCM]=brain_plot3d(CM,cod,pp,CMold)
% IN:
%   CM: square matrix, a (potentially sparsified) adjacency matrix for a graph
%   cod: struct, a codebook annotating regions
%   pp: struct, plotting parameters
%   CMold: square matrix, original full adjacency matrix
%
% OUT:
%   myHs: vector, handles to all surface objects
%   myHt: vector, handles to all text objects

myHs=[];    % surface objects
myHt=[];    % text objects

nVertices=size(CM,1);

wCM=zeros(nVertices)*nan;

%posOffset=0; %100;
nCylPts=20;

% if (pp.createFigure==true)
%     mH=figure;
%     hold on;
%     caxis([0 1]);
% else
%     mH=[];
%     hold on;
%     caxis([0 1]);
% end
set(gca,'Color',pp.BGColour,'XColor',pp.compColour,'YColor',pp.compColour,...
    'ZColor',pp.compColour);
set(gcf,'Color',pp.BGColour);

%% plot region centroids and vertex degrees

if pp.computeDegreesOnFullMat==true
    % check symmetry
    if (all(all(CMold==CMold')))
        CMsym=CMold;
    else
        CMsym=(CMold+CMold')/2;
    end
else
    if (all(all(CM==CM')))
        CMsym=CM;
    else
        CMsym=(CM+CM')/2;
    end
end

% compute soft degrees
%myDegrees=sum(triu(CMsym,1),2);
myDegrees=sum(abs(CMsym)); % use abs weight
[sXbase,sYbase,sZbase] = sphere(30); % coord for sphere -> draw sphere with surf (filled) or mesh (grid)
for iter=1:cod.num,
    cx=cod.center{iter}(1)+pp.posOffset(1);
    cy=cod.center{iter}(2)+pp.posOffset(2);
    cz=cod.center{iter}(3)+pp.posOffset(3);
    if isempty(pp.FctForSizSph)
        thisDeg=myDegrees(iter)+1;
    else % use value i fct to size
        thisDeg=pp.FctForSizSph(iter)+1;
    end
    %plot3(cx,cy,cz,'o','MarkerEdgeColor',pp.compColour./3,...
    %    'MarkerSize',10); % thisDeg
    
    if pp.doPlotDegreeSpheres
        % XXX the else/if below seem unnecessarily complicated
        if abs(thisDeg-1)>=pp.plotDegreeSpheresOnlyAboveThreshold*(1/nVertices)
            foo_plotDS=true;
        else
            foo_plotDS=false;
        end
        
        if ~isempty(pp.FctForColSph) && pp.SphBlack==true % OVERRIDE ABOVE, plots all
                foo_plotDS=true;
        end
        
%         if isempty(pp.FctForColSph) % use degree to decide if to plot
%             if abs(thisDeg-1)>=pp.plotDegreeSpheresOnlyAboveThreshold*(1/nVertices)
%                 foo_plotDS=true;
%             else
%                 foo_plotDS=false;
%             end
%         elseif ~isempty(pp.FctForColSph) 
%             if pp.SphBlack==true % plots all
%                 foo_plotDS=true;
%             else  % uses degree as above
%                 if ((thisDeg-1)>=pp.plotDegreeSpheresOnlyAboveThreshold*(1/nVertices))
%                     foo_plotDS=true;
%                 else
%                     foo_plotDS=false;
%                 end
%                 %if (abs(pp.FctForColSph(iter))>=pp.plotDegreeSpheresOnlyAboveThreshold)
%                 %    foo_plotDS=true;
%                 %else
%                 %    foo_plotDS=false;
%                 %end
%             end
%         end
        %         else
        %             foo_plotDS=true;
        %         end
        if foo_plotDS==true
            thisDegS=(thisDeg*pp.f_matD_to_sphS)^pp.e_matD_to_sphS;%^4;
            wCM(iter,iter)=thisDegS;
            sX=(thisDegS*sXbase)+cx;
            sY=(thisDegS*sYbase)+cy;
            sZ=(thisDegS*sZbase)+cz;
            if ~isempty(pp.FctForColSph) % plot function on vertex instead of degree
                myCData=pp.FctForColSph(iter);
                if ~isempty(pp.ColAxis), caxis(pp.ColAxis);
                end
            else
                myCData=((thisDegS/pp.f_matD_to_sphS)/(nVertices+1)); % normalise between 0 and 1 : max degree is nVertices
            end
            if pp.SphBlack==true && ~(abs(pp.FctForColSph(iter))>=pp.plotDegreeSpheresOnlyAboveThreshold)
                myCoData=ones([size(sZ,1) size(sZ,2) 3]);
                myCoData=myCoData*pp.GrayShade; % same color as connections
            else
                myCoData=repmat(myCData,size(sZ));
            end
            myHs(end+1)=surf(sX,sY,sZ,'FaceAlpha',1.0,'CData',myCoData,'edgecolor','none');
        end
    end
end;

%% colormap
% some colormaps need to be inverted so that low edge weights are bright
% and high vals are dark
invCMs={'autumn','bone','hot','gray','copper','pink'};

if all(cellfun('isempty',strfind(invCMs,pp.cmapchoiceV)))
    colormap(pp.cmapchoiceV);
else
    colormap(flipud(colormap(pp.cmapchoiceV)));
end

if strfind(pp.cmapchoiceV,'hsv')
    colormap([hsv;hsv]);
end

freezeColors;

%% plot connections and weights
edgeWeights=zeros(size(CM),'single');
plotThisConnection=false;
disp('Processing nodes...');
for iter=1:cod.num,
    if (mod(iter,10)==0)
        fprintf('%s ',num2str(iter));
        %disp(['Processing node ' num2str(iter) '...']);
    end
    for iter2=iter+1:cod.num,
        if ~isnan(CM(iter,iter2)) && CM(iter,iter2)~=0
            % compute display params
            myWraw=abs(CM(iter,iter2));
            myW=myWraw^pp.gamma;
            edgeWeight=pp.linearWeight*myW; %^2;
            edgeWeights(iter,iter2)=edgeWeight;
            %colVal=myW;
            colVal = CM(iter,iter2) ^ pp.gamma;
            
            if pp.pruneTinyConnections==true && myW<pp.pruningWeight
                plotThisConnection=false;
                %disp(['Pruning ' num2str(myW)]);
            else
                plotThisConnection=true;
            end
          %  disp([iter, iter2, colVal])
            if plotThisConnection==true
                wCM(iter,iter2)=edgeWeight; % save for returning outside of show_cm
                if abs(colVal) < pp.pruningWeight, colVal=pp.pruningWeight; end
                if abs(colVal) > 0.99, colVal=sign(colVal)*0.99; end
                if pp.doAlpha==true
                    alphaVal=myWraw^pp.alphagamma; % ^3 for more rapid decrease
                    if alphaVal<pp.pruningWeight, alphaVal=pp.pruningWeight; end
                    if alphaVal>0.99, alphaVal=0.99; end
                else
                    alphaVal=0.99;
                end
                
                
                %h=plot3([cod.center{iter}(1) cod.center{iter2}(1)],[cod.center{iter}(2) cod.center{iter2}(2)],[cod.center{iter}(3) cod.center{iter2}(3)],'b-');
                %set(h,'LineWidth',edgeWeight*3);
                %set(h,'Color',[colVal colVal colVal]);
                [fooX, fooY, fooZ] = cylinder2P(edgeWeight, ...
                    nCylPts,[cod.center{iter}]'+pp.posOffset',...
                    [cod.center{iter2}]'+pp.posOffset');
                if (pp.doColourCode==true)
                    if pp.useSplitCMap==true
                        % map split 1 to 0-0.5 and split 2 to 0.5-1
                        colVal=(colVal/2)+0.5*(pp.splitID-1);
                        % what falls on 0.5 is mapped down for split 1, up for
                        % split 2
                        if (colVal==0.5)
                            if (pp.splitID==1)
                                colVal=colVal-eps;
                            else
                                colVal=colVal+eps;
                            end
                        end
                    end
                    myCdata=repmat(colVal,2,nCylPts);
                else
                    % assign a fixed scalar CDATA value to this vertex
                    if pp.useSplitCMap==true
                        % assign different CDATA values to different splits
                        myCdata=repmat(pp.CdataVal+1*(pp.splitID-1),2,nCylPts);
                    else
                        % all vertices use same CDATAvalue (pp.CdataVal)
                        if pp.CdataVal==0 % show connections as black
                            myCdata=ones([2 nCylPts 3]);
                            myCdata=myCdata*pp.GrayShade;
                        else
                            myCdata=repmat(pp.CdataVal,2,nCylPts);
                        end
                    end
                end
                if (pp.doPlotConnections==true)
                    
                    % each side is defined by XData YData ZData (all 2xnCylPts), where
                    % 3D points are defined by X(i,j) Y(i,j) Z(i,j)
                    % and lines joining them are X(:,j) Y(:,j) Z(:,j)
                    % sides can be drawn by surface(X(:,1:2),Y(:,1:2),Z(:,1:2))
                    % force double for cData to avoid complaint - maybe
                    % input data is single...
                    myHs(end+1)=surf(fooX,fooY,fooZ,'FaceAlpha',alphaVal,'EdgeAlpha',alphaVal','CData',double(myCdata));
                else
                    %myHs(end+1)=[];
                end
                if (pp.showWeights==true)
                    midPt=mean([cod.center{iter}+pp.posOffset cod.center{iter2}+pp.posOffset],2);
                    myHt(end+1)=text(midPt(1),midPt(2),midPt(3),[num2str(colVal,'%1.2f')]); % ...
                    %    '=>' num2str(edgeWeight,'%1.2f')]);
                end
            end % plotThisConnection
        end;
    end;
end;




%% plot labels
% compute subset to show (long-winded way)
if ~isempty(pp.showTopNlabels) && pp.showLabels
    labelsToShow=zeros(cod.num,1);
    if strcmp(pp.showTopNlabelsBy,'edgeW')
        edgeWeights_v=jUpperTriMatToVec(edgeWeights,1);
        edgeWeights_v_srt=sort(edgeWeights_v,'descend');
        cutoffW=edgeWeights_v_srt(pp.showTopNlabels);
        if cutoffW==0 % there might be fewer non-zero edges than showTopNlabels
            warning('Not enough non-zero weights for pp.showTopNlabels parameter, selecting fewer');
            nonZeroWeights=find(edgeWeights_v_srt);
            if ~isempty(nonZeroWeights)
                cutoffW=edgeWeights_v_srt(nonZeroWeights(end));
            else
                cutoffW=Inf;
            end
        end
        % point out edges whose nodes are worth labelling
        connsToShow=edgeWeights>=cutoffW;
        % ugly
        for i=1:cod.num
            for j=i:cod.num
                if (connsToShow(i,j)==true)
                    labelsToShow(i)=1;
                    labelsToShow(j)=1;
                end
            end
        end
    elseif strcmp(pp.showTopNlabelsBy,'vertexW')
        % XXX hack
        myDegrees_srt=sort(myDegrees,'descend'); % myDegrees
        cutoffW=myDegrees_srt(pp.showTopNlabels);
        labelsToShow=myDegrees >= cutoffW;
    else
        error('Unknown option for plotting labels by');
    end
else
    labelsToShow=ones(cod.num,1);
end
% do the display
if (pp.showLabels==true)
    for iter=1:cod.num
        if (labelsToShow(iter)==1)
            processedLabel=cod.name{iter};
            if pp.flipLabelsRL
                if strcmp(processedLabel(end),'L')
                    processedLabel(end)='R';
                else
                    processedLabel(end)='L';
                end
            end
            myHt(end+1)=text(double(cod.center{iter}(1)+pp.posOffset(1)+pp.txtPosOffset(1)),...
                double(cod.center{iter}(2)+pp.posOffset(2)+pp.txtPosOffset(2)),...
                double(cod.center{iter}(3)+pp.posOffset(3)+pp.txtPosOffset(3)),...
                processedLabel,...
                'FontWeight','normal','FontSize',pp.myFontSize,...
                'Interpreter','none','HorizontalAlignment','center', ...
                'Color',pp.compColour);
        end
    end
end

fprintf('%s\n','Done');
axis equal;
hold off;
%view(66,10);
grid off;

lighting flat;
shading flat;

% colormap edge
if all(cellfun('isempty',strfind(invCMs,pp.cmapchoiceE)))
    colormap(pp.cmapchoiceE);
else
    colormap(flipud(colormap(pp.cmapchoiceE)));
end

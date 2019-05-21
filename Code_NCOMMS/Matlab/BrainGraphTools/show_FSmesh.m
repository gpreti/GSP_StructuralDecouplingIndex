function myH=show_FSmesh(hemi,meshType,varargin)
% SHOW_FSMESH show a freesurfer 3-dimensional hemisphere mesh
% 
% IN:
%   hemi - a char array containing 'l','r', or 'lr', indicating which hemi-
%       sphere to display
%   meshtype - a char array containing mesh type ('white')
%   'faceColor': ColorSpec - colour for the faces of the brain mesh
%       triangles
%   'faceAlpha': scalar between 0 and 1 - alpha transparency value to
%       use for brain mesh triangles (0: fully transparent, 1: fully
%       opaque)
%   'edgeColor': ColorSpec - colour for the edges of the brain mesh
%       triangles
%   'edgeAlpha': scalar between 0 and 1 - alpha transparency value to
%       use for brain mesh triangle edges (0: fully transparent, 1: fully
%       opaque)
%   'reducePatchVal': scalar, fraction of patches to keep - increasing
%   from the default value may crash Matlab due to a bug in Matlab.
% OUT:
%   myH: a vector of handles to patch objects
%
% v1.0 Nov 2009 Jonas Richiardi
% - inital release based on read_surf/trisurf snippet by Markus Gschwind
% v1.1 2009 Dec 10 Jonas Richiardi
% - code clean up
% - fix for alpha values crashing in 2-hemispheres display (known bug
% 484445): use subdivision algorithm to "downsample" Delaunay triangulation 
% This is a better solution in this case than those suggested in
% http://www.mathworks.com/support/bugreports/484445
% v1.2 2010 Jul Jonas Richiardi
% - fix a bug due to reducepatch does not setting CData and FaceVertexCData
% after reducing the vertex and face count (despite the doc saying it does
% not set the Faces and Vertices properties). Now CData and FaceVertexCData
% are in sync.

[createFigure,FaceColor,EdgeColor,FaceAlpha,EdgeAlpha,useConvexHull,...
    reducePatchVal]=process_options(varargin,...
    'createFigure',true,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.5 0.5 0.5],...
    'FaceAlpha',0.2,'EdgeAlpha',0.05,'useConvexHull',false,'reducePatchVal',...
    0.015);

%% sanity check
hemiSupported={'l','r','lr'};
meshTypeSupported={'white'};
if all(cellfun('isempty',strfind(hemiSupported,hemi)))
    error(['Unsupported hemi type: ' hemi]);
end
if all(cellfun('isempty',strfind(meshTypeSupported,meshType)))
    error(['Unsupported mesh type: ' meshType]);
end


%% setup host-specific options and paths
location=whereAmIRunning;
if strcmp(location,'jmac')
    % set root path to FreeSurfer application
    FSROOTPATH='/Applications/freesurfer/';
else
    error('Please modify whereAmIRunning.m for your own settings');
end
% generate path where to find the surface mesh files
FSSURFPATH=fullfile(FSROOTPATH,'subjects','fsaverage','surf');

if exist('read_surf')~=2
    addpath(fullfile(FSROOTPATH,'matlab'));
end

if (createFigure==true)
    figure;
end
hold on;


%% load surfaces
% surfaces are encoded by Delaunay triangulation (region-adjacency graph
% of Voronoi tesselation - Delaunay graph is the dual of the Voronoi
% diagram (Okabe et al 200))
% faces is the m x 3 triangulation matrix - the set of simplices (triangle)
% that make the triangulation. Each row (in simplex-vertex format) encodes
% the three vertices a triangular face via indexing into a matrix containing
% the 3D space location of vertices
% e.g. 1 4 2 says that this triangles is made up of vertices 1, 4, and 2,
% whose 3D spatial locations are found in the vertices matrix.
myH=zeros(numel(hemi),1);
for h=1:numel(hemi)
    if strcmp(hemi(h),'r')
        f_surf=['rh.' meshType];
    elseif strcmp(hemi(h),'l')
        f_surf=['lh.' meshType];
    else
        error('Wrong hemi specification. How on earth did you trigger this code path?');
    end
    % read surface using freesurfer function
    % note we could access matlab functionality for the TriRep class via
    % trirep(faces+,vertices(:,1)... etc)
    [vertices,faces]=read_surf(fullfile(FSSURFPATH,f_surf));
    vertices=single(vertices); % save memory
    nVertices = size(vertices,1);
    if (useConvexHull==true)
        disp('Computing convex hull...');
        myConvHull = convhulln(double(vertices));
        myH(h)=trisurf(myConvHull+1,vertices(:,1),vertices(:,2),vertices(:,3),...
            'FaceColor',FaceColor,'EdgeColor',EdgeColor);
        disp('Doing mesh subdivision...');
        %reducepatch(myH(h),0.05);
    else
        % plot Delaunay triangulated mesh
        myH(h)=trisurf(faces+1,vertices(:,1),vertices(:,2),vertices(:,3),...
            ones(nVertices,1),'FaceColor',FaceColor,'EdgeColor',EdgeColor);
        
        % do subdivision to reduce number of faces
        % another solution would be to set(gca,'DrawMode','Fast'); as per
        % 484445 bug report
        disp('Doing mesh subdivision...');
        reducepatch(myH(h),reducePatchVal); % higher may cause crash, still too many faces

        % CData and FaceVertexCData are now desynchronised, so fix them
        % manually
        nf=get(myH(h),'Faces');
        nv=get(myH(h),'Vertices');
        % 1) fix CData to be 3xsize(Faces,1), all ones or actually FaceColor
        set(myH(h),'CData',repmat(FaceColor,3,size(nf,1)));
        % 2) fix FaceVertexCData to be nVertices x 1
        set(myH(h),'FaceVertexCData',repmat(FaceColor,size(nv,1),1));
    end
    
    xlabel('left-right'); ylabel('posterior-anterior'); zlabel('ventral-dorsal');
    clear vertices faces;
    if (h==1)
        axis equal;
        view(0,90);
    end
end
%set(gca,'DrawMode','Fast');

% set transparency now, with reduced number of faces (does not crash)
for h=1:numel(hemi)
    set(myH(h),'EdgeAlpha',EdgeAlpha,'FaceAlpha',FaceAlpha,...
        'DiffuseStrength',1,'SpecularStrength',0);
end



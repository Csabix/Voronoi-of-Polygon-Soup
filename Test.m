classdef Test
properties (Constant)
    vertexRegionColor    = rgb(209, 171, 125);
    edgeRegionColor      = rgb(107, 153, 195);
    lineColor            = rgb(2, 46, 102);
    parabolaColor        = rgb(99, 79, 60);
    tripointColor        = rgb(194, 195, 197);
    polygonColor         = rgb(201, 25, 89);
    delaunayColor        = Test.edgeRegionColor;
    graphEdgeColor       = Test.vertexRegionColor;

    lineWidth            = 3;
    pointSize            = 3;
    polyPtSize           = 4;
    polyFaceAlpha        = 0;

    useColorbar          = false;
    useAxis              = false;
    bringTripoints2Front = false;
    colorTripoints       = false;

    drawRegions          = true;
    drawTripoints        = true;
    drawDelaunay         = false;
    drawnGraphEdges      = 1;     % 0: none, 1: full graph, 2: medial axis, 3: inner MA, 4: outer MA
    truncateMedialAxis   = false; % if true no MA edges will touch the polygon
end

methods (Static)

function p = regularNGon(n)
    assert(n>=3);
    t = linspace(0,2*pi,n+1); t = t(1:end-1);
    p = [cos(t') sin(t')];
end

function p = randRadPoly(n)
    assert(n>=3);
    radii = 1 + 2*rand(n,1);
    th = rand(3,1)/3 + [0;1/3;2/3];
    thetas = 2*pi * sort([th; rand(n-3,1)]);
    p = [cos(thetas).*radii, sin(thetas).*radii];
end

function input = rotated_square_grid(n,m)
    arguments
        n (1,1) double {mustBePositive}
        m (1,1) double {mustBePositive} = n 
    end
    sq = [-.5 -.5; .5 -.5; .5 .5; -.5 .5]; % unit square
    input = cell(n*m,1);
    ths = 2*pi*rand(n*m,1);
    for i = 0:m-1
        for j = 0:n-1
            id = n*i+j+1;
            a = ths(id);
            input{id} = sq * [cos(a) sin(a); -sin(a) cos(a)] + 2*[i j];
        end
    end
end

function pts = concat_polys(args)
    arguments
        args (:,1) cell
    end
    n = length(args) - 1 + sum(cellfun(@(a) size(a,1),args));
    pts = NaN(n,2);
    startId = 1;
    for arg = args'
        num = size(arg{1},1);
        pts(startId:startId+num-1,:) = arg{1};
        startId = startId + num + 1;
    end
end

function [pts,Ps] = drawRandPoly(varargin)
    rng(0); Ps = {}; Test.concatCallback; % reset n
    if length(varargin) == 1 && any(isnan(varargin{1}),"all")
        data = varargin{1};
        rows = size(data,1);
        nanIndices = [find(isnan(data(:,1))); rows + 1];
        last = 0;
        for i = nanIndices'
            startIdx = last + 1;
            endIdx = i - 1;
            last = i;
            if endIdx - startIdx < 2; continue; end
            Ps{end+1} = Polygon(data(startIdx:endIdx,:),Test.polygonColor,LineWidth=Test.lineWidth,FaceAlpha=Test.polyFaceAlpha,MarkerSize=Test.polyPtSize);
        end
    else
        for i = 1:length(varargin)
            n = varargin{i};
            if length(n) > 1
                P = n;
            else
                P = Test.randRadPoly(n) + [6*i-6,0];
            end
            Ps{end+1} = Polygon(P,Test.polygonColor,LineWidth=Test.lineWidth,FaceAlpha=Test.polyFaceAlpha,MarkerSize=Test.polyPtSize);
        end
    end
    pts = CustomValue('polyArr',Ps{:}, @Test.concatCallback);

    % Create a dependent delaunay triangulation structure
end

function poly = readLandmark(k)
    persistent data;
    if isempty(data)
        data = readstruct('data/Historic_Landmarks.geojson','FileType','json');
        data = data.features;
    end
    poly = cell2mat(data(k).geometry.coordinates{1}');
end

function pts = concatCallback(varargin)
    pts = [];
    for arg = varargin
        pts = [pts; arg{1}(1:end-1,:); NaN NaN];
    end
    persistent n;
    if nargin == 0; n = []; return; end
    len = size(pts, 1) - numel(varargin);
    if isempty(n) || n ~= len
        n = len;
        clim([1 2*n]);
        %col = net(haltonset(3),2*n);
        col = rand(2*n,3);
        col(1:n,:)     = ones(n,1)*Test.vertexRegionColor;  % col(1:n,:)    .*[.3 .3 .3] + [.7 .4 .0];
        col(n+1:end,:) = ones(n,1)*Test.edgeRegionColor;    % col(n+1:end,:).*[.3 .3 .3] + [.0 .4 .7];
        %colormap(col*.5+.5);
        colormap(col);
        if Test.useColorbar
            colorbar('Ticks',linspace(1.5,2*n-0.5,2*n),'TickLabels',1:2*n);
        end
        if Test.useAxis
            axis on;
        else
            axis off;
        end
    end
end

function id = selectTriangle(CD)
    % triangle selector point
    S = Point('selectorPt',[0 0],'r','LabelVis','off');
    % selected triangle index
    id = Scalar('triangleId',CD,S,@(CD,S) pointLocation(CD.dt,S));
    % selected circle
    Circle(Point(CD,id,@(cd,id) cd.dCenters(id,:),'Vis','off'),...
           Scalar(CD,id,@(cd,id) cd.dRadii(id)));
end

function id = selectVertex(CD)
    % vertex selector point
    S = Point('selectorPt',[0 0],'r','LabelVis','off');
    id = Scalar('pointId',CD,S,@(CD,S) nearestNeighbor(CD.dt,S));
    % selected vertex
    Point('selectedPt',CD,id, @(CD, id) CD.pVertices(id,:),'g',16,'LabelVis','off');
end

function id = selectEdge(CD)
    % edge selector point
    S = Point('selectorPt',[0 0],'r','LabelVis','off');
    id = Scalar('edgeId',CD,S,@Test.nearestEdge);
    % selected edge
    SegmentSequence('selectedEdge',CD,id, @(CD, id) [CD.pVertices(CD.pEdges(id,1),:); CD.pVertices(CD.pEdges(id,2),:)],'g',15);
end

function inds = regioncolor(x,y,cd)
    dist = inf(size(x));
    inds = zeros(size(x));
    a = cd.pVertices(cd.pEdges(:,1),:);
    b = cd.pVertices(cd.pEdges(:,2),:);
    N = size(a,1);
    abx = b(:,1)-a(:,1); aby = b(:,2)-a(:,2);
    abilen = 1./sqrt(abx.*abx + aby.*aby);
    nx =  aby.*abilen;   ny = -abx.*abilen;
    
    for i = 1:N
        ax = a(i,1)-x; ay = a(i,2)-y;
        dvert = ax.*ax + ay.*ay;
        inds(dvert < dist) = i;
        dist = min(dist,dvert);

        bx = b(i,1)-x; by = b(i,2)-y;
        mask = ax.*abx(i) + ay.*aby(i) < 0 & bx.*abx(i) + by.*aby(i) > 0;
        dedge = (nx(i).*ax+ny(i).*ay).^2;
        mask = mask & dist > dedge;
        inds(mask) = i+N;
        dist(mask) = dedge(mask);
    end
end

function id = nearestEdge(cd, p)
    a = cd.pVertices(cd.pEdges(:,1),:);
    b = cd.pVertices(cd.pEdges(:,2),:);
    d2 = geometry.distSquared2seg(p,a,b);
    [~,id] = min(d2);
end



end % methods (Static)
end % classdef Test


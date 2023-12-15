addpath Geomatplot\
%% test tripoint graph creation
clf; disp 'TripointGraph test'

testcase = 'glyph';
glyph = 'A';
n = 50; % 2^n !!!
m = 4;
switch testcase
    case 'glyph'
        % read glyph file
        glyphs = readGlyphs(['glyphs/' glyph '.glyph']);
        poly = glyphsToPolygons(glyphs, 2);
        input = {poly};
    case 'rand'
        % random polygons
        input = {n};
    case 'ngon'
        % regular polygons
        input = {Test.regularNGon(n)};
    case 'given'
        %input = {[0 0; 1 0; .3 1; .4 .3] * 4 - 2};
        %input = {[0 0; 2 0; 2 1; 0 1]};
        %input = {[0 0; 0 1; 1 1; 1 0]}; % square / 4-point
        %input = {[0 0; 0 1; 2 1; 2 0]}; % rectangle
        %input = {[0 0; -.1 0.5; 0 1; 1 1; 1 0]}; % 4-point
        %input = {[0 0; -.1 0.5; 0 1; 1 1; 1.1 0.5; 1 0]}; % 2 neighbouring 4-points
        %input = {[0 0; 0 1; 1 1; 1 .75; .25 .75; .25 .25; 1 .25; 1 0]}; % 4-point -- ideal edge
        %input = {[-5 -1; -5 1; 5 1; 5 2; 6 2; 6 -2; 5 -2; 5 -1],[-.2 -.2; 0 .2; .2 -.2]}; % double shared region boundary connected to 5-point
        input = {[1 -1; 0.3 -2.1; -1 -1; -3 -1; -3 1; -1 1; 0.1 2.1; 1 1; -0.4 0]}; % VEE parallel test
    case '344'
        input = {...
            [-1.4755    1.7781;   -2.2012   -2.1051;    1.7330   -2.1731], ...
            [4.3250    2.1733;    5.5894   -0.6875;    8.8136   -2.3786;    8.5924    1.9362], ...
            [15.5726    2.0911;   15.0181   -0.5615;   11.7534   -2.2047;   17.2211   -2.1500]};
    case 'cross'
        input = {[1 -1; 3 -1; 3 1; 1 1; 1 3; -1 3; -1 1; -3 1; -3 -1; -1 -1; -1 -3; 1 -3]}; % cross, 9 4-points
    case 'fivepoint'
        input = {[0 -1; 2 -1; 2 1; 0 1; 0 2; -1 2; -1 -2; 0 -2]}; % 5-point
    case 'hexa'
        input = {[-1.2473    2.5024;   -1.0243    0.0407;   -2.0827   -2.6493;    0.0332   -0.9764;    3.6047   -1.0391;    1.0858   -0.3389]}; % concave hexagon
    case 'rect'
        input = {[-3 -1; -3 1; 3 1; 3 -1],[-.3 -.3; 0 .3; .3 -.3]}; % small triangle in long rectangle
    case 'squaregrid'
        sq = [-.5 -.5; .5 -.5; .5 .5; -.5 .5]; % unit square
        input = cell(n^2,1);
        for i = 0:n-1
            for j = 0:n-1
                input{n*i+j+1} = sq + 2*[i j];
            end
        end
    case 'rotgird_n_m'
        input = Test.rotated_square_grid(n,m);
    case 'rotgrid'
        sq = [-.5 -.5; .5 -.5; .5 .5; -.5 .5]; % unit square
        input = cell(n^2,1);
        ths = 2*pi*rand(n^2,1);
        for i = 0:n-1
            for j = 0:n-1
                id = n*i+j+1;
                a = ths(id);
                input{id} = sq * [cos(a) sin(a); -sin(a) cos(a)] + 2*[i j];
            end
        end
    case 'squares'
        sq = [-.5 -.5; .5 -.5; .5 .5; -.5 .5]; % unit square
        input = cell(n,1);
        for i = 1:n
            input{i} = sq * (2*i-1);
        end
    otherwise
        error('Unknown test case');
end

[Poly,Ps] = Test.drawRandPoly(input{:});

xyRange = [min(Poly.value,[],1); max(Poly.value,[],1);];
center = mean(xyRange,1);
margin = min(0.3*(xyRange(2,:) - center));
xyRange = xyRange + [-margin -margin; margin margin];
xlim(xyRange(:,1)'); ylim(xyRange(:,2)');

% Create a dependent delaunay triangulation structure
CD = CustomValue('connectionData',Poly,@ConnectionData);

% Draw Image
%c1 = Point('c1',xyRange(1,:),'b',5); c2 = Point('c2',xyRange(2,:),'b',5);
%Image(CD,@Test.regioncolor,c1,c2,'Res',512);

% Create Triangulation graph
TG = CustomValue('TipointGraph',CD, @TripointGraph);
% profile trigraph = TripointGraph(ConnectionData(Poly.value));

% create regions
g = TG.parent;
if Test.drawRegions
    for i = 1:2*size(CD.val.pVertices,1)
        lab = g.getNextLabel('region');
        cb = getRegionCB(i,100);
        args = struct;
        args.Color = i;
        args.Linestyle = 'none';
        args.Facealpha = 1;
        dpolygon(g,lab,{TG},cb,args);
    end
end

% Draw tripoint graph
if Test.drawnGraphEdges == 1 % full graph
    %   linear edges
    SegmentSequence('lines',TG,@callb_graphLines,'-',Test.lineWidth,Color=Test.lineColor);
    %   parabolic edges
    bezierSubdiv = 100;
    SegmentSequence('parabolas',TG,@(tg) callb_graphParabolas(tg,bezierSubdiv),bezierSubdiv,'-',Test.lineWidth,Color=Test.parabolaColor);
elseif Test.drawnGraphEdges >= 2 % medial axis
    if Test.drawnGraphEdges == 2 % full
        if Test.truncateMedialAxis
            edgeFilter = @(tg) isMAEdge(tg) & ~isTouchingEdge(tg);
        else
            edgeFilter = @(tg) isMAEdge(tg);
        end
    elseif Test.drawnGraphEdges == 3 % inner
        if Test.truncateMedialAxis
            edgeFilter = @(tg) isInEdge(tg) & isMAEdge(tg) & ~isTouchingEdge(tg);
        else
            edgeFilter = @(tg) isInEdge(tg) & isMAEdge(tg);
        end
    elseif Test.drawnGraphEdges == 4 % outer
        if Test.truncateMedialAxis
            edgeFilter = @(tg) isOutEdge(tg) & isMAEdge(tg) & ~isTouchingEdge(tg);
        else
            edgeFilter = @(tg) isOutEdge(tg) & isMAEdge(tg);
        end
    end
    %   linear edges
    SegmentSequence('lines',TG,@(tg) callb_graphLines_filtered(tg,edgeFilter(tg)),'-',Test.lineWidth,Color=Test.lineColor);
    %   parabolic edges
    bezierSubdiv = 100;
    SegmentSequence('parabolas',TG,@(tg) callb_graphParabolas_filtered(tg,edgeFilter(tg),bezierSubdiv),bezierSubdiv,'-',Test.lineWidth,Color=Test.lineColor);
end

if Test.drawDelaunay
    % Draw the triangulation edges
    % the edges(dt) function returns a n x 2 index matrix we need to transpose for the right ordering
    delseq = SegmentSequence('delaunay',CD,@(CD) CD.pVertices(CD.dEdges',:),'c--',Test.lineWidth,'Color',Test.delaunayColor);
    topDelaunay = @() uistack(delseq.fig,'top');
    bottomPolygon = @() cellfun(@(ps) uistack(ps.fig,'bottom'),Ps);
    topDelaunay();
    PointSequence(CD,@(cd) cd.pVertices,Test.polygonColor,Test.pointSize);
    bottomPolygon();
end
% draw tripoints
if Test.drawTripoints
    if ~Test.colorTripoints
        tpts = PointSequence('tripoints',TG, @(TG) TG.tPos,Test.tripointColor,Test.lineWidth);
        tripointFigs = tpts.fig;
    else
        tpts1 = PointSequence('tripoints1',TG, @(TG) TG.tPos(TG.tSignedDist<0,:),[1 0 0],Test.lineWidth);
        tpts2 = PointSequence('tripoints2',TG, @(TG) TG.tPos(TG.tSignedDist>0,:),[0 1 0],Test.lineWidth);
        tpts3 = PointSequence('tripoints3',TG, @(TG) TG.tPos(TG.tSignedDist==0,:),[0 0 1],Test.lineWidth);
        tripointFigs = [tpts1.fig; tpts2.fig; tpts3.fig];
    end
    polyFigs = cellfun(@(p)p.fig,Ps);
    topTris = @() uistack(tripointFigs, 'top');
    topPolys = @() uistack(polyFigs, 'top');
    
    if Test.bringTripoints2Front
        topTris();
    end
end
xyRange2 = xyRange';
%a = annotation('rectangle',[.23 .1 .575 .83],'Color',0*[.99 .99 .99]);
%uistack(a,'bottom');
xlim(xyRange(:,1)'); ylim(xyRange(:,2)');

%%

function edges = callb_graphLines(tg)
    edges = reshape(tg.gpEdge(tg.gEdgeIsLinear,1:4)',2,[])';
end
function edges = callb_graphLines_filtered(tg, mask)
    edges = reshape(tg.gpEdge(tg.gEdgeIsLinear & mask,1:4)',2,[])';
end
% mask for medial axis edges
function mask = isMAEdge(tg)
    mask = ((tg.gEdgeFoots(:,1)>0) == (tg.gEdgeFoots(:,2)>0)) | ~tg.gEdgeIsLinear;
end
% mask for inside graph edges
function mask = isInEdge(tg)
    mask = tg.tSignedDist(tg.gEdges(:,1)) <= 0 & tg.tSignedDist(tg.gEdges(:,2)) <= 0;
end
% mask for outside graph edges
function mask = isOutEdge(tg)
    mask = tg.tSignedDist(tg.gEdges(:,1)) >= 0 & tg.tSignedDist(tg.gEdges(:,2)) >= 0;
end
% mask for graph edges touching the polygon
function mask = isTouchingEdge(tg)
    mask = tg.tSignedDist(tg.gEdges(:,1)) == 0 | tg.tSignedDist(tg.gEdges(:,2)) == 0;
end

function edges = callb_graphParabolas(tg,N)
    B = tg.gpEdge(~tg.gEdgeIsLinear,1:6);
    t = linspace(0,1,N);
    Px = (B(:,1) * (1-t).*(1-t) + B(:,5) * 2*t.*(1-t) + B(:,3) * t.*t)';
    Py = (B(:,2) * (1-t).*(1-t) + B(:,6) * 2*t.*(1-t) + B(:,4) * t.*t)';
    edges = [Px(:) Py(:)];
end
function edges = callb_graphParabolas_filtered(tg,mask,N)
    B = tg.gpEdge(~tg.gEdgeIsLinear & mask,1:6);
    t = linspace(0,1,N);
    Px = (B(:,1) * (1-t).*(1-t) + B(:,5) * 2*t.*(1-t) + B(:,3) * t.*t)';
    Py = (B(:,2) * (1-t).*(1-t) + B(:,6) * 2*t.*(1-t) + B(:,4) * t.*t)';
    edges = [Px(:) Py(:)];
end

function cb = getRegionCB(i,sub)
    function [X,Y] = CB(TG)
        [X,Y] = TripointGraph.getRegionBoundary(TG.value.gRegions(i),sub);
    end
    cb = @CB;
end
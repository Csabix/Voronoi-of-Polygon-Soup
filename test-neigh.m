addpath Geomatplot\
%% Circle-edges intersection test
clf; disp 'Circle-edges intersection test'

% create random polygon
Poly = Test.drawRandPoly(30);
% Create a dependent delaunay triangulation structure
CD = CustomValue('connectionData',Poly,@ConnectionData);
% selected triangle index
id = Test.selectTriangle(CD);

% the edges that intersect the circle
SegmentSequence(CD,id, @(CD,id)CD.getEdgeSequence(CD.edgesInCircle{id}),'r-',6);

%% Vertex neightbour vertex test
clf; disp 'Vertex neightbour vertex test'

% create random polygon
Poly = Test.drawRandPoly(30);
% Create a dependent delaunay triangulation structure
CD = CustomValue('connectionData',Poly,@ConnectionData);
% selected vertex index
id = Test.selectVertex(CD);

% neightbouring vertices
seq1 = PointSequence(CD,id, @(CD,id) CD.pVertices(CD.dVertNeigbourVert{id},:),'m',10);

%% nearbyEdgesV
clf; disp 'nearbyEdgesV test'

% create random polygon
Poly = Test.drawRandPoly(30);
% Create a dependent delaunay triangulation structure
CD = CustomValue('connectionData',Poly,@ConnectionData);
% selected vertex index
id = Test.selectVertex(CD);

% nearbyEdgesV
SegmentSequence(CD,id, @(CD,id)CD.getEdgeSequence(CD.nearbyEdgesV{id}),'r',7);

%% nearbyEdgesE
clf; disp 'nearbyEdgesE test'

% create random polygon
Poly = Test.drawRandPoly(30);
% Create a dependent delaunay triangulation structure
CD = CustomValue('connectionData',Poly,@ConnectionData);
% selected edge index
id = Test.selectEdge(CD);

% nearbyEdgesE
SegmentSequence(CD,id, @(CD,id)CD.getEdgeSequence(CD.nearbyEdgesE{id}),'r',7);






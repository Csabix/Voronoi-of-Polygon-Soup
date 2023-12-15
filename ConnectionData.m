classdef ConnectionData

properties
    pVertices           (:,2) double % pVertId  -> [x,y]
    pEdges              (:,2) int32  % pEdgeId  -> [start pVertId, end pVertId]
    pVertEdges          (:,2) int32  % pVertId  -> [pEdgeId, pEdgeId]  the two connecting edges at the vertex

    dTriangles          (:,3) int32  % dTriId   -> [vertId,vertId,vertId]  tringle vertices
    dEdges              (:,2) int32  % dEdgeId  -> [vertId,vertId]
    dTriangleNeighbours (:,3) int32  % dTriId   -> [dTriId,dTriId,dTriId]  tringle indices of neightbour
    dVertexNeighbours   (:,1) cell   % vertId   -> {[dTriId;...]}
    dCenters            (:,2) double % dTriId   -> [x,y]
    dRadii              (:,1) double % dTriId   -> r
    edgesInCircle       (:,1) cell   % dTriId   -> {[pEdgeId;...]}   polygon edges in the Delaunay circles
    dVertNeigbourVert   (:,1) cell   % vertId   -> {[vertId;...]}
    dCirclesOnEdge      (:,1) cell   % pEdgeId  -> {[dTriId;...]}
    convHull            (:,1) int32  % id       -> pVertId

    nearbyEdgesV        (:,1) cell   % vertId   -> {[pEdgeId;...]}   polygon edges near this vertex
    nearbyEdgesE        (:,1) cell   % pEdgeId  -> {[pEdgeId;...]}   polygon edges near this edge
    trianglesNearEdge   (:,1) cell   % pEdgeId  -> {[dTriId;...]}    triangles near this edge

    dt % the Delaunay Triangulation
    ps % the polyshape object
end

methods
    function o= ConnectionData(poly)
        arguments
            poly (:,2) double
        end
        if nargin == 0; return; end
        o.ps = polyshape(poly);
        o.pVertices = o.ps.Vertices(~any(isnan(o.ps.Vertices),2),:);

        [o.pEdges,o.pVertEdges] = calc_pEdges(o);
        DT = delaunayTriangulation(o.pVertices);
        o.dt = DT;
        o.dTriangles = DT.ConnectivityList();
        o.dEdges = DT.edges();
        o.dTriangleNeighbours = DT.neighbors();
        o.dVertexNeighbours = DT.vertexAttachments();
        [o.dCenters, o.dRadii] = DT.circumcenter();
        [o.edgesInCircle,o.dCirclesOnEdge] = calc_edgesInCircle_and_dCirclesOnEdge(o);
        o.dVertNeigbourVert = calc_dVertNeigbourVert(o);
        o.convHull = convexHull(DT);

        o.nearbyEdgesV = calcNearbyEdgesV(o);
        [o.nearbyEdgesE,o.trianglesNearEdge] = calcNearbyEdgesE(o);
    end
    
    % for drawing a set of edges
    function edges = getEdgeSequence(o, pEdgeIds)
        arguments
            o (1,1) ConnectionData
            pEdgeIds (:,1) int32
        end
        es = o.pEdges(pEdgeIds,:)';
        edges = o.pVertices(es(:),:);
    end
end

methods (Access=private)

    function [pEdges,pVertEdges] = calc_pEdges(o)
        psv = [o.ps.Vertices;NaN NaN];
        N = numsides(o.ps);
        pEdges = zeros(N,2);
        pVertEdges = zeros(N,2);
        currFirstV = 1; psVertId = 1;        
        for pEdgeId = 1:N
            pVertEdges(pEdgeId,2) = pEdgeId;
            if all(isnan(psv(psVertId+1,:)))
                pEdges(pEdgeId,:) = [pEdgeId currFirstV];
                pVertEdges(currFirstV,1) = pEdgeId;
                psVertId    = psVertId + 1;
                currFirstV = pEdgeId + 1;
            else
                pEdges(pEdgeId,:) = [pEdgeId pEdgeId+1];
                pVertEdges(pEdgeId+1,1) = pEdgeId;
            end
            psVertId = psVertId + 1;
        end
    end

    function [edgesInCircle,dCirclesOnEdge] = calc_edgesInCircle_and_dCirclesOnEdge(o)
        edgesInCircle = cell(size(o.dTriangles,1),1);
        dCirclesOnEdge = cell(size(o.pEdges,1),1);
        for triId = 1:size(o.dTriangles,1)
            edgesInCircle{triId} = ConnectionData.intersect_circle2segments(o.dCenters(triId,:),o.dRadii(triId),...
                o.pVertices(o.pEdges(:,1),:),o.pVertices(o.pEdges(:,2),:));
            for pEdgeId = edgesInCircle{triId}'
                dCirclesOnEdge{pEdgeId}(end+1) = triId;
            end
        end
    end

    function dVertNeigbourVert = calc_dVertNeigbourVert(o)
        dVertNeigbourVert = cell(size(o.dVertexNeighbours,1),1);
        for vertId = 1:size(o.dVertexNeighbours,1)
             ids = o.dTriangles(o.dVertexNeighbours{vertId},:);
             ids = unique(ids(:));
             dVertNeigbourVert{vertId} = ids(ids~=vertId);
        end
    end

    function nearbyEdgesV = calcNearbyEdgesV(o)
        arguments
            o (1,1) ConnectionData
        end
        numVerts = size(o.pVertices,1);
        nearbyEdgesV = cell(numVerts,1);
        for vertId = 1:numVerts
            nearbyEdgesV{vertId} = unique(vertcat(o.edgesInCircle{o.dVertexNeighbours{vertId}}));
        end
    end

    function [nearbyEdgesE,trianglesNearEdge] = calcNearbyEdgesE(o)
        arguments
            o (1,1) ConnectionData
        end
        numEdges = size(o.pEdges,1);
        nearbyEdgesE = cell(numEdges,1);
        trianglesNearEdge = cell(numEdges,1);
        for pEdgeId = 1:numEdges
            tris = o.dCirclesOnEdge{pEdgeId};
         %   neighTris = o.dTriangleNeighbours(tris,:);
         %   tris = union(neighTris(neighTris~=0)',tris);
            es = vertcat(o.edgesInCircle{tris});
            es = setdiff(es, pEdgeId);
            nearbyEdgesE{pEdgeId} = es;
            trianglesNearEdge{pEdgeId} = tris;
        end
    end
end

methods (Static)
    % returns the indices of segments (between `a` and `b` Nx2) which intersect the
    % (center,radius) cirle
    % note that if a segment is contained in the circle it is not returned
    function indices = intersect_circle2segments(center, radius, a, b)
        EPS = 1e-8;
        ab = b-a;
        ac = a - center;
        A = dot(ab,ab,2);
        B = 2*dot(ab,ac,2);
        d = B.*B-4.*A.*(dot(ac,ac,2)-radius.^2);        
        L = d >= 0;
        d  = L .* sqrt(d);
        t = [L;L] .* ( -0.5*([B;B]+[d;-d])./[A;A] );
        L = [L;L]  &  -EPS <= t  &  t <= 1 + EPS;
        L = L(1:end/2) | L(end/2+1:end);
        indices = int32(find(L));
    end
end

end


classdef TripointGraph
properties
    % tripoints
    tPos   (:,2) double % tId -> [x,y]
    tType  (:,3) int32  % tId -> [type,type,type]   type = 0: infinity, 1: vertex, 2: edge, monoton incr.
    tFoots (:,3) int32  % tId -> [id,id,id]         id = iId/pVertId/pEdgeId (generating countour elements)
    tSignedDist (:,1) double  % tId -> sDist        sDist = signed distance to the polygon (negative inside)
    tTypeNums (1,6) int32  % VVV,VVE,VEE,EEE,IVV,IVE -> num of tripoints
    tIdealDir (:,2) double % iId -> [vx,vy]  ray (to the ideal tripoint)
    % graph
    gEdges        (:,2) int32   % gEdgeId -> [tId,tId]
    gEdgeFoots    (:,2) int32   % gEdgeId -> [id,id]   id = -pVertId / pEdgeId (generating countour elements)
    % graph edge drawing helpers
    gEdgeIsLinear (:,1) logical % gEdgeId -> bool      true: linear, false: parbolic
    gpEdge        (:,6) double  % gEdgeId -> [ax ay bx by cx cy]   a: starting point, b: endpoint, c: 2nd bezier control point (or NaN if linear)
    % regions
    gRegions      (:,1) struct  % gRegionId -> region struct  (gRegionId: 1..N -> vertices, N+1..2*N -> edges)

    % input connection data
    CD     (1,1) %ConnectionData 
end

methods
function o = TripointGraph(CD)
    arguments
        CD (1,1) ConnectionData
    end
    o.CD = CD;
    
    % create tripoints
    [o.tPos, o.tType, o.tFoots, o.tTypeNums, o.tIdealDir, o.tSignedDist] = calcTripoints(o);
    
    checkTripointCount(o);
    
    % create graph edges
    [o.gEdges,o.gEdgeFoots,o.tPos,o.tSignedDist] = calcGraph(o);

    [o.gEdgeIsLinear, o.gpEdge] = calcGrapgEdgeDrawingHelpers(o);

    o.gRegions = calcRegions(o);
end

function [tPos, tType, tFoots, tTypeNums, tIdealDir, tSDist] = calcTripoints(o)
    tTypeNums = [0 0 0 0 0 0];
    % VVV
    [VVVPos,VVVd]   = o.CD.dt.circumcenter();
    mask = true(size(VVVPos,1),1);
    for dTriId = 1:size(o.CD.dTriangles,1)
        newTPos = VVVPos(dTriId,:);
        dist2 = VVVd(dTriId).^2;
        edges2check = o.CD.pEdges(o.CD.edgesInCircle{dTriId},:);
        edgeV1Pos = o.CD.pVertices(edges2check(:,1),:);
        edgeV2Pos = o.CD.pVertices(edges2check(:,2),:);
        edgeVExclude = o.CD.dTriangles(dTriId,:);
        
        mask(dTriId) = TripointGraph.checkTripointDistFromEdges(newTPos,dist2,edges2check,edgeV1Pos,edgeV2Pos,edgeVExclude);
    end
    VVVPos = VVVPos(mask,:);
    tTypeNums(1) = size(VVVPos,1);
    VVVType  = zeros(tTypeNums(1),3) + [1 1 1];
    VVVFoots = o.CD.dTriangles(mask,:);
    VVVSDist = VVVd(mask);
    % tripoint is outside <=> foot vertex is convex
    for i = 1:tTypeNums(1)
        vId = VVVFoots(i,1);
        es = o.CD.pVertEdges(vId,:);
        vIds = [o.CD.pEdges(es(1),1), vId, o.CD.pEdges(es(2),2)];
        Vs = o.CD.pVertices(vIds,:);
        norm = glm.normal(Vs(1,:), Vs(2,:));
        isConvex = norm*(Vs(1,:) - Vs(3,:))' > 0;
        VVVSDist(i) = VVVSDist(i)*(2*isConvex - 1);
    end
    
    % VVE
    VVEPos   = [];
    VVEType  = [1 1 2];
    VVEFoots = [];
    VVESDist = [];
    for vert1Id = 1:size(o.CD.pVertices,1)
        vert1Pos = o.CD.pVertices(vert1Id,:);
        neighVs = o.CD.dVertNeigbourVert{vert1Id};
        for vert2Id = neighVs'
            if vert1Id >= vert2Id; continue; end
            vert2Pos = o.CD.pVertices(vert2Id,:);
            neighTriIds = o.CD.dt.edgeAttachments(double(vert1Id), double(vert2Id)); % 1 or 2 triangles
            neighEs = unique(vertcat(o.CD.edgesInCircle{neighTriIds{1}}));
            checkVertices = setdiff(reshape(o.CD.dTriangles(neighTriIds{1}',:),1,[]),[vert1Id,vert2Id]);
            for edge3Id = neighEs'
                edge3VIds = o.CD.pEdges(edge3Id,:);
                vert31Pos = o.CD.pVertices(edge3VIds(1),:);
                vert32Pos = o.CD.pVertices(edge3VIds(2),:);
                
                p1_online = vert1Id == edge3VIds(1) || vert1Id == edge3VIds(2);
                p2_online = vert2Id == edge3VIds(1) || vert2Id == edge3VIds(2);
               
                if p1_online && p2_online
                    % the points are the endpoints of the segment
                    continue; % no new point
                elseif p1_online && ~p2_online
                    % the first point is an endpoint of the segment
                    n = glm.normal(vert31Pos,vert32Pos);
                    TPos = geometry.equidistRP(vert1Pos,n,vert2Pos);
                elseif ~p1_online && p2_online
                    % the second point is an endpoint of the segment
                    n = glm.normal(vert31Pos,vert32Pos);
                    TPos = geometry.equidistRP(vert2Pos,n,vert1Pos);
                elseif ~p1_online && ~p2_online
                    % general case: both points are off the segment
                    [newTPos1, newTPos2] = geometry.equidistPPL(vert1Pos,vert2Pos,vert31Pos,vert32Pos);
                    if ~geometry.pointIsOverEdge(newTPos1, vert31Pos, vert32Pos)
                        newTPos1 = [NaN NaN];
                    end
                    if ~geometry.pointIsOverEdge(newTPos2, vert31Pos, vert32Pos)
                        newTPos2 = [NaN NaN];
                    end
                    TPos = [newTPos1; newTPos2];
                else
                    warning('shouldn''t get here');
                end
                
                for temp = TPos' % 1 or 2 new points
                    newTPos = temp';
                    if any(isnan(newTPos)); continue; end
                    dist2 = glm.dot2(newTPos-vert1Pos);

                    % check distance of all close vertices
                    
                    vertices2check = checkVertices(checkVertices ~= vert1Id & checkVertices ~= vert2Id);
                    vertPos = o.CD.pVertices(vertices2check,:);
                    
                    keep = TripointGraph.checkTripointDistFromVertices(newTPos,dist2,vertPos);
                    if ~keep; continue; end
                    
                    % check distance of all close edges

                    edges2check = o.CD.pEdges(neighEs(neighEs~=edge3Id),:);
                    edgeV1Pos = o.CD.pVertices(edges2check(:,1),:);
                    edgeV2Pos = o.CD.pVertices(edges2check(:,2),:);
                    edgeVExclude = [vert1Id, vert2Id];

                    keep = TripointGraph.checkTripointDistFromEdges(newTPos,dist2,edges2check,edgeV1Pos,edgeV2Pos,edgeVExclude);
                    if ~keep; continue; end

                    norm = glm.normal(vert31Pos, vert32Pos);
                    isOutside = norm*(newTPos - vert31Pos)' > 0;

                    % save newTPos
                    VVEPos = vertcat(VVEPos, newTPos);
                    VVEFoots = vertcat(VVEFoots, [vert1Id, vert2Id, edge3Id]);
                    VVESDist = vertcat(VVESDist, (2*isOutside-1)*sqrt(dist2));
                end
            end
        end
    end
    tTypeNums(2) = size(VVEPos,1);
    VVEType  = zeros(tTypeNums(2),3) + VVEType;
    
    % VEE
    VEEPos =   [];
    VEEType =  [1 2 2];
    VEEFoots = [];
    VEESDist = [];
    for vert1Id = 1:size(o.CD.pVertices,1)
        vert1Pos = o.CD.pVertices(vert1Id,:);
        neigh1Es = o.CD.nearbyEdgesV{vert1Id};
        vert1NeighVs = o.CD.dVertNeigbourVert{vert1Id};
        for edge2Id_id = 1:size(neigh1Es,1)-1
            edge2Id = neigh1Es(edge2Id_id);
            edge2VIds = o.CD.pEdges(edge2Id,:);
            vert21Pos = o.CD.pVertices(edge2VIds(1),:);
            vert22Pos = o.CD.pVertices(edge2VIds(2),:);

            for edge3Id_id = edge2Id_id+1 : size(neigh1Es,1)
                edge3Id = neigh1Es(edge3Id_id);
                edge3VIds = o.CD.pEdges(edge3Id,:);
                vert31Pos = o.CD.pVertices(edge3VIds(1),:);
                vert32Pos = o.CD.pVertices(edge3VIds(2),:);
            
                % 1 point, 2 lines
                P1isP21 = vert1Id == edge2VIds(1);
                P1isP22 = vert1Id == edge2VIds(2);
                P1isP31 = vert1Id == edge3VIds(1);
                P1isP32 = vert1Id == edge3VIds(2);
                summa = P1isP21 + P1isP22 + P1isP31 + P1isP32;
                if summa == 0
                    % general case: the point is off the segments
                    [n1, n2] = geometry.equidistPLL(vert1Pos, vert21Pos, vert22Pos, vert31Pos, vert32Pos);
                    edges2Check = [vert21Pos, vert22Pos; vert31Pos, vert32Pos];
                elseif summa == 1
                    % the point is the endpoint of one of the segments
                    if P1isP21
                        n = glm.normal(vert1Pos,vert22Pos);
                        [n1,n2] = geometry.equidistRL(vert1Pos,n,vert31Pos,vert32Pos);
                        edges2Check = [vert31Pos, vert32Pos];
                    elseif P1isP22
                        n = glm.normal(vert1Pos,vert21Pos);
                        [n1,n2] = geometry.equidistRL(vert1Pos,n,vert31Pos,vert32Pos);
                        edges2Check = [vert31Pos, vert32Pos];
                    elseif P1isP31
                        n = glm.normal(vert1Pos,vert32Pos);
                        [n1,n2] = geometry.equidistRL(vert1Pos,n,vert21Pos,vert22Pos);
                        edges2Check = [vert21Pos, vert22Pos];
                    elseif P1isP32
                        n = glm.normal(vert1Pos,vert31Pos);
                        [n1,n2] = geometry.equidistRL(vert1Pos,n,vert21Pos,vert22Pos);
                        edges2Check = [vert21Pos, vert22Pos];
                    else
                        warning('shouldn''t get here');
                    end
                else
                    % the segments share an endpoint (and it is the first parameter)
                    % save trivial tripoint newTPos
                    VEEPos = vertcat(VEEPos, vert1Pos);
                    VEEFoots = vertcat(VEEFoots, [vert1Id, edge2Id, edge3Id]);
                    VEESDist = vertcat(VEESDist, 0);
                    continue;
                end
                
                ns = [n1;n2];
                for temp = ns' % 1 or 2 new points
                    newTPos = temp';
                    if any(isnan(newTPos)); continue; end
                    
                    keep = true; 
                    for checkEdge = edges2Check'
                        if ~geometry.pointIsOverEdge(newTPos, checkEdge(1:2)', checkEdge(3:4)')
                            keep = false;
                            break;
                        end
                    end 
                    if ~keep; continue; end
                    
                    dist2 = glm.dot2(newTPos-vert1Pos);

                    % check distance of close vertices
                    
                    vertices2check = vert1NeighVs(vert1NeighVs ~= vert1Id);
                    vertPos = o.CD.pVertices(vertices2check,:);
                    
                    keep = TripointGraph.checkTripointDistFromVertices(newTPos,dist2,vertPos);
                    if ~keep; continue; end

                    % check distance of close edges

                    edges2check = o.CD.pEdges(neigh1Es(neigh1Es~=edge2Id & neigh1Es~=edge3Id),:);
                    edgeV1Pos = o.CD.pVertices(edges2check(:,1),:);
                    edgeV2Pos = o.CD.pVertices(edges2check(:,2),:);
                    edgeVExclude = [vert1Id];

                    keep = TripointGraph.checkTripointDistFromEdges(newTPos,dist2,edges2check,edgeV1Pos,edgeV2Pos,edgeVExclude);
                    if ~keep; continue; end

                    norm = glm.normal(vert31Pos, vert32Pos);
                    isOutside = norm*(newTPos - vert31Pos)' > 0;

                    % save newTPos
                    VEEPos = vertcat(VEEPos, newTPos);
                    VEEFoots = vertcat(VEEFoots, [vert1Id, edge2Id, edge3Id]);
                    VEESDist = vertcat(VEESDist, (2*isOutside-1)*sqrt(dist2));
                end
            end
        end
    end
    tTypeNums(3) = size(VEEPos,1);
    VEEType  = zeros(tTypeNums(3),3) + VEEType;
    
    % EEE
    EEEPos =   [];
    EEEType =  [2 2 2];
    EEEFoots = [];
    EEESDist = [];
    for edge1Id = 1:size(o.CD.pEdges,1)
    	edge1VIds = o.CD.pEdges(edge1Id,:);
        vert11Pos = o.CD.pVertices(edge1VIds(1),:);
        vert12Pos = o.CD.pVertices(edge1VIds(2),:);
        neigh1Es = o.CD.nearbyEdgesE{edge1Id};
        e1Tris = o.CD.trianglesNearEdge{edge1Id};
        neigh1Vs = unique(reshape(o.CD.dTriangles(e1Tris,:),[],1));
        neigh1VPoss = o.CD.pVertices(neigh1Vs,:);
        for edge2Id_id = 1:size(neigh1Es,1)-1
            edge2Id = neigh1Es(edge2Id_id);
            if edge2Id <= edge1Id; continue; end
            edge2VIds = o.CD.pEdges(edge2Id,:);
            vert21Pos = o.CD.pVertices(edge2VIds(1),:);
            vert22Pos = o.CD.pVertices(edge2VIds(2),:);
            for edge3Id_id = edge2Id_id+1 : size(neigh1Es,1)
                edge3Id = neigh1Es(edge3Id_id);
                edge3VIds = o.CD.pEdges(edge3Id,:);
                vert31Pos = o.CD.pVertices(edge3VIds(1),:);
                vert32Pos = o.CD.pVertices(edge3VIds(2),:);
                
                % disp([edge1Id, edge2Id, edge3Id]);
                
                newTPos = geometry.equidistLLL(vert11Pos,vert12Pos,vert21Pos,vert22Pos,vert31Pos,vert32Pos);
                
                edges2Check = [vert11Pos,vert12Pos; vert21Pos,vert22Pos; vert31Pos,vert32Pos];
                
                keep = true; 
                for checkEdge = edges2Check'
                    if ~geometry.pointIsOverEdge(newTPos, checkEdge(1:2)', checkEdge(3:4)')
                        keep = false;
                        break;
                    end
                end 
                if ~keep; continue; end
                
                dist2 = geometry.distSquared2seg(newTPos, vert11Pos, vert12Pos);

                % check distance of close vertices
                
                keep = TripointGraph.checkTripointDistFromVertices(newTPos,dist2,neigh1VPoss);
                if ~keep; continue; end
                
                % check distance of close edges

                edges2check = o.CD.pEdges(neigh1Es(neigh1Es~=edge1Id & neigh1Es~=edge2Id & neigh1Es~=edge3Id),:);
                edgeV1Pos = o.CD.pVertices(edges2check(:,1),:);
                edgeV2Pos = o.CD.pVertices(edges2check(:,2),:);
                edgeVExclude = [];

                keep = TripointGraph.checkTripointDistFromEdges(newTPos,dist2,edges2check,edgeV1Pos,edgeV2Pos,edgeVExclude);
                if ~keep; continue; end

                norm = glm.normal(vert31Pos, vert32Pos);
                isOutside = norm*(newTPos - vert31Pos)' > 0;
                
                % save newTPos
                EEEPos = vertcat(EEEPos, newTPos);
                EEEFoots = vertcat(EEEFoots, [edge1Id, edge2Id, edge3Id]);
                EEESDist = vertcat(EEESDist, (2*isOutside-1)*sqrt(dist2));
            end
        end
    end
    tTypeNums(4) = size(EEEPos,1);
    EEEType = zeros(tTypeNums(4),3) + EEEType;
    
    % ideal tripoints
    hull = o.CD.convHull;
    hullCentroid = mean(o.CD.pVertices(hull(1:end-1),:),1);
    hullEdges = [hull(1:end-1), circshift(hull(1:end-1),1)];
    nonOriginalHullEdges = setdiff(hullEdges, o.CD.pEdges,'rows');
    originalHullEdges = setdiff(hullEdges,nonOriginalHullEdges,'rows');
    
    numIdealsIVV = size(nonOriginalHullEdges,1);
    numIdealsIVE = 2*size(originalHullEdges,1);
    numIdeals = numIdealsIVV + numIdealsIVE;
    tTypeNums(5) = numIdealsIVV;
    tTypeNums(6) = numIdealsIVE;

    tIdealDir = zeros(numIdeals, 2); % [vx vy]
    IdePos = zeros(numIdeals, 2);
    IdeType = [zeros(numIdealsIVV,3) + [0 1 1]; zeros(numIdealsIVE,3) + [0 1 2]];
    IdeFoots = zeros(numIdeals, 3);
    IdeSDist = repmat(+Inf, numIdeals, 1); % ideal tripoints are outside

    idealIndex = 1;
    % IVV
    for edge = nonOriginalHullEdges'
        ab = o.CD.pVertices(edge,:);
        n = - glm.normalPosSide(ab(1,:),ab(2,:),hullCentroid);
        tIdealDir(idealIndex,:) = n;
        IdePos(idealIndex,:) = (ab(1,:)+ab(2,:))/2 + 100 * n;
        IdeFoots(idealIndex,:) = [idealIndex edge'];
        idealIndex = idealIndex+1;
    end

    % IVE
    for e = originalHullEdges'
        edge = sort(e);
        ab = o.CD.pVertices(edge,:);
        n = - glm.normalPosSide(ab(1,:),ab(2,:),hullCentroid);
        % calc edge index: if v1 + 1 == v2 -> e = v1; else e = v2
      % if edge(1) + 1 == edge(2); edgeId = edge(1); else; edgeId = edge(2); end
        edgeId = edge(1 + (edge(1) + 1 ~= edge(2)));

        tIdealDir(idealIndex,:) = n;
        IdePos(idealIndex,:) = ab(1,:) + 100*n;
        IdeFoots(idealIndex,:) = [idealIndex edge(1) edgeId];
        idealIndex = idealIndex+1;

        tIdealDir(idealIndex,:) = n;
        IdePos(idealIndex,:) = ab(2,:) + 100*n;
        IdeFoots(idealIndex,:) = [idealIndex edge(2) edgeId];
        idealIndex = idealIndex+1;
    end


    % concatenate all tripoints
    tPos   = [VVVPos  ; VVEPos  ; VEEPos  ; EEEPos  ; IdePos  ];
    tType  = [VVVType ; VVEType ; VEEType ; EEEType ; IdeType ];
    tFoots = [VVVFoots; VVEFoots; VEEFoots; EEEFoots; IdeFoots];
    tSDist = [VVVSDist; VVESDist; VEESDist; EEESDist; IdeSDist];
end

function p = getFootEntityMidpoint(o,footIds)
    p = zeros(length(footIds),2);
    for i=1:length(footIds)
        footId = footIds(i);
        if footId < 0 % vertex
            p(i,:) = o.CD.pVertices(-footId,:);
        elseif footId > 0 % edge
            p(i,:) = mean(o.CD.pVertices(o.CD.pEdges(footId,:),:));
        else
            assert(false);
        end
    end
end

function checkTripointCount(o)
    n = size(o.CD.pVertices,1);
    calcAll = sum(o.tTypeNums);
    shouldBe = 4*n-2;
    if calcAll < shouldBe
        warning('Missing tripoints!');
        disp (o.tTypeNums);
        disp (['sum = ', num2str(calcAll)]);
        disp (['should be = ', num2str(shouldBe)]);
    end
end

function [gEdges,gEdgeFoots,tPos,tSignedDist] = calcGraph(o)
    % N = size(o.CD.pVertices,1);
    % [0 1 2] -> [0 -1 1]
    tPos = o.tPos;
    tSignedDist = o.tSignedDist;
    mappedfoots = o.tFoots.*(2*o.tType-3).*int32(o.tType ~= 0);
    [tIdPairs,gEdgeFoots,tdIdRest,gRestFoots] = TripointGraph.findPairsWithTwoCommonElements(mappedfoots);

    threeRegulars = [];
    extrafoots = [];
    mergeTripoints = zeros(length(mappedfoots),1); % tId -> setId (setId == 0: not in any set)
    numMergeSets = 0; % largest setId
    for i = 1:length(tdIdRest)
        tIds = tdIdRest{i};         % length>=4 (if 3 regular)
        if length(tIds) < 2; continue; end
        trispos = o.tPos(tIds,:);
        hasMergedPoints = false;
        for j = 1:length(trispos)-1
            for k = j+1:length(trispos)
                if glm.dot2(trispos(j,:)-trispos(k,:)) < 1e-10
                    ids = tIds([j k]);
                    SIds = mergeTripoints(ids);
                    if SIds(1) == 0 && SIds(2) == 0 % new set
                        numMergeSets = numMergeSets+1;
                        mergeTripoints(ids) = numMergeSets;
                    elseif SIds(1) == 0 % insert tripoint1 to set
                        mergeTripoints(ids(1)) = SIds(2);
                    elseif SIds(2) == 0 % insert tripoint2 to set
                        mergeTripoints(ids(2)) = SIds(1);
                    else
                        if SIds(1) ~= SIds(2) % merge the two sets
                            mergeTripoints(mergeTripoints == SIds(2)) = SIds(1);
                        end % else: they are already in the same set
                    end
                    
                    hasMergedPoints = true;
                end
            end
        end
        if ~hasMergedPoints && length(tIds) > 2
            threeRegulars = [threeRegulars i];
            extrafoots = [extrafoots; gRestFoots(i,:)];
        else
            
        end
    end

    if numMergeSets ~= 0
        mergedTripointFoots = {};
        % merge tripoints
        for setId = 1:numMergeSets
            bb = mergeTripoints==setId;
            foots = unique(mappedfoots(bb,:));
            if isempty(foots); continue; end
            triPos = mean(o.tPos(bb,:));
            tPos(end+1,:) = triPos;
            tSignedDist(end+1,:) = mean(o.tSignedDist(bb));
            assert(all(foots~=0));
            mappedfoots(bb,:) = NaN;

            footEntityMidPoints = o.getFootEntityMidpoint(foots);
            footDirs = footEntityMidPoints - triPos;
            mappedDirectionsToLine = atan2(footDirs(:,2),footDirs(:,1));
            [~,sortedIds] = sort(mappedDirectionsToLine);

            mergedTripointFoots{end+1} = foots(sortedIds);
        end

        [tIdPairs,gEdgeFoots,tdIdRest,gRestFoots] = TripointGraph.findPairsWithTwoCommonElements(mappedfoots,mergedTripointFoots,mergeTripoints~=0);
        threeRegulars = [];
        extrafoots = [];
        for i = 1:length(tdIdRest)
            if length(tdIdRest{i}) > 2
                threeRegulars = [threeRegulars i];
                extrafoots = [extrafoots; gRestFoots(i,:)];
            end
        end
    end
    
    % Sometimes two regions have a multiple common boundaries (fairly rare).
    % In this case tdIdRest cell array contains indices for at least 4 tripoints with two common
    % regions. At least one of these is an edge region due to vertex region's convexity, so we can
    % project the 4 or more tripoints onto this edge and connect every pair in sorted order.
    for i = threeRegulars
        tIds = tdIdRest{i};         % length>=4 and mod(length,2) == 0
        idNum = length(tIds); 
        assert(idNum >= 4 && mod(idNum,2) == 0);
        trispos = tPos(tIds,:);
        boots = extrafoots(i,:);
        boots = boots(boots>0);     % select edge. only need one if there was two
        ab = o.CD.pVertices(o.CD.pEdges(boots(1),:)',:);
        tvals = geometry.paramOnLine(trispos,ab(1,:),ab(2,:));
        [~,idx] = sort(tvals);
        tIdPairs = [tIdPairs; reshape(tIds(idx),2,[])']; %#ok<AGROW>
        gEdgeFoots = [gEdgeFoots; repmat(extrafoots(i,:),idNum/2,1)]; %#ok<AGROW>
    end

    % filter ideal-ideal connections
    nonIdealMask = all(gEdgeFoots ~= 0,2);
    tIdPairs = tIdPairs(nonIdealMask,:);
    gEdgeFoots = gEdgeFoots(nonIdealMask,:);

    gEdges = tIdPairs;
end

function [gEdgeIsLinear, gpEdge] = calcGrapgEdgeDrawingHelpers(o)
    N = size(o.gEdges,1);
    gEdgeIsLinear = zeros(N,1,'int32');
    gpEdge = zeros(N,6);
    for i = 1:N
        tIds = o.gEdges(i,:);
        B0 = o.tPos(tIds(1),:);
        B2 = o.tPos(tIds(2),:);
        gpEdge(i,1:4) = [B0 B2];
        B1 = 0.5 * B0 + 0.5 * B2; % midpoint
        footIds = o.gEdgeFoots(i,:);
        isLinear = true;
        if footIds(1) < 0 && footIds(2) > 0
            % possibly parabolic
            edgeVertices = o.CD.pEdges(footIds(2),:);
            if all(edgeVertices ~= -footIds(1))
                isLinear = false;
                P = o.CD.pVertices(edgeVertices,:);
                F = o.CD.pVertices(-footIds(1),:);
                tpB0 = geometry.paramOnLine(B0, P(1,:), P(2,:));
                pB0 = (1-tpB0) * P(1,:) + tpB0 * P(2,:);
                tpB2 = geometry.paramOnLine(B2, P(1,:), P(2,:));
                pB2 = (1-tpB2) * P(1,:) + tpB2 * P(2,:);

                B1 = geometry.lineintersect(0.5*(pB0+F),F-pB0,0.5*(pB2+F),F-pB2);
            end
        end
        gEdgeIsLinear(i) = isLinear;
        gpEdge(i,5:6) = B1;
    end
end

% regionId: 1..N: pEdgeId, -N..-1: -pVertId
function region = getRegion(o, regionId)
    region = struct;
    edgesMask = any(o.gEdgeFoots == regionId,2);
    edges = o.gEdges(edgesMask,:);                  % localEdgeId -> [tId, tId]
    edgeMapping = find(edgesMask);                  % localEdgeId -> gEdgeId
    dic = configureDictionary('int32','cell');      % tId -> [localEdgeId, localEdgeId]
    for i = 1:size(edges,1)
        for j = 1:2
            key = edges(i,j);
            if isKey(dic,key)
                dic{key} =  [dic{key} i];
            else
                dic{key} = i;
            end
        end
    end
    vals = values(dic);                             % dicId -> {[localEdgeId,...]}
    dicKeys = keys(dic);                            % dicId -> tId
    lens = cellfun(@numel,vals);
    assert(~any(lens>2));
    localEdgeIdPairs = cell2mat(vals(lens==2));     % dic2Id -> [localEdgeId,localEdgeId]
    pairKeys = dicKeys(lens==2);                    % dic2Id -> tId
    localEdgeIdSingles = cell2mat(vals(lens==1));   % dic1Id -> localEdgeId
    singleKeys = dicKeys(lens==1);                  % dic1Id -> tId
    numSingles = size(localEdgeIdSingles,1);
    assert(numSingles == 0 || numSingles == 2);
    loopEndsAtInfinity = (numSingles == 2);
    % first edge and vertex on the loop
    vertLoop = zeros(size(dicKeys,1),1);
    if loopEndsAtInfinity
        currentEdge = localEdgeIdSingles(1);          % localEdgeId
        currentVert = singleKeys(1);                  % tId
        edgeLoop = zeros(size(dicKeys,1)-1,1);
    else
        currentEdge = localEdgeIdPairs(1,1);          % localEdgeId
        currentVert = pairKeys(1);                    % tId
        edgeLoop = zeros(size(dicKeys,1),1);
    end
    edgeLoop(1) = currentEdge;
    vertLoop(1) = currentVert;
    % calculate loop
    for loopId = 1:size(dicKeys)-1
        edgeVerts = edges(currentEdge,:);                    % [tId, tId]
        nextVert = edgeVerts(edgeVerts ~= currentVert);      % tId
        vertLoop(loopId+1) = nextVert;
        currentVert = nextVert;
        connectedEdges = dic{currentVert};                      % [localEdgeId, localEdgeId] or [localEdgeId]
        nextEdge = connectedEdges(connectedEdges~=currentEdge); % localEdgeId or []
        if isempty(nextEdge)
            break;
        end
        edgeLoop(loopId+1) = nextEdge;
        currentEdge = nextEdge;
    end
    edgeIds = edgeMapping(edgeLoop);
    region.gEdgeIsLinear = o.gEdgeIsLinear(edgeIds);
    region.gpEdge = o.gpEdge(edgeIds,:);
    % set the edges in the correct orientation
    for i = 1:size(edgeLoop,1)
        if edges(edgeLoop(i),1) ~= vertLoop(i)
            region.gpEdge(i,1:4) = region.gpEdge(i,[3 4 1 2]);
        end
    end
    % close the loop
    if loopEndsAtInfinity
        a = region.gpEdge(end,3:4);
        b = region.gpEdge(1,1:2);
        region.gpEdge(end+1,:) = [a, b, (a+b)/2];
        region.gEdgeIsLinear(end+1) = true;
    end
end

function regions = calcRegions(o)
    numVerts = size(o.CD.pVertices,1);
    regions = repmat(struct('gEdgeIsLinear',[],'gpEdge',[]),2*numVerts,1);
    % vertex regions
    for i = 1:numVerts
        regions(i) = o.getRegion(-i);
    end
    % edge regions
    for i = 1:numVerts
        regions(numVerts + i) = o.getRegion(i);
    end
end

end

methods (Static)
    
function keep = checkTripointDistFromEdges(newTPos,dist2,edges2check,edgeV1Pos,edgeV2Pos,edgeVExclude)
    keep = true;
    for i = 1:size(edges2check,1)
        edgeCheckVIds = edges2check(i,:);
        [dist2_to_check,closestIndex] = geometry.distSquared2seg_meta(newTPos, edgeV1Pos(i,:), edgeV2Pos(i,:));
        if closestIndex ~= 0
            closestPointID = edgeCheckVIds(closestIndex);
            if any(edgeVExclude == closestPointID)
                continue;
            end
        end
        if dist2_to_check < dist2 - 1e-9
            keep = false;
            break;
        end
    end
end
function keep = checkTripointDistFromVertices(newTPos,dist2,vertPos)
    keep = all(sum((vertPos-newTPos).^2,2) >= dist2 - 1e-9);
end
% nontriplets are assumed to be ordered (so the common elements should be adjacent)
% if skipTriplets is true for an index I, triplets(I,:) is skipped in the algorithm
function [tIdPairs,gEdgeFoots,tdIdRest,gRestFoots] = findPairsWithTwoCommonElements(triplets,nontriplets,skipTriplets)
    arguments
        triplets     (:,3) int32
        nontriplets  (:,1) cell = {}
        skipTriplets (:,1) logical = []
    end
    doSkip = false;
    if ~isempty(skipTriplets)
        doSkip = true;
        assert(size(skipTriplets,1) == size(triplets,1));
    end

    triplets = sort(triplets,2);
    idmax = max(triplets(:,3));
    idmin = min(triplets(:,1));
    if ~isempty(nontriplets)
        idmax = max(idmax, max(cellfun(@max,nontriplets)));
        idmin = min(idmin, min(cellfun(@min,nontriplets)));
    end
    N = (idmax-idmin+1);
    dic = configureDictionary('int32','cell');
    dicKeys = triplets(:,[1 2 1])-idmin + N*(triplets(:,[2 3 3])-idmin); % key = a-m + N*(b-m)
    for tId = 1:size(triplets,1)
        if doSkip && skipTriplets(tId); continue; end
        for i = 1:3
            key = dicKeys(tId,i);
            if isKey(dic,key)
                dic{key} =  [dic{key} tId];
            else
                dic{key} = tId;
            end
        end
    end
    numTripoints = size(triplets,1);
    for tId2 = 1:size(nontriplets,1)
        nont = nontriplets{tId2};
        n = length(nont);
        for i = 1:n
            j = i+1; if j > n; j = 1; end
            ab = sort(nont([i j]));
            key = ab(1)-idmin + N*(ab(2)-idmin); % key = a-m + N*(b-m)
            if isKey(dic,key)
                dic{key} =  [dic{key} tId2+numTripoints];
            else
                dic{key} = tId2+numTripoints;
            end
        end
    end
    vals = values(dic); lens = cellfun(@numel,vals); 
    tIdPairs = cell2mat(vals(lens==2));
    tdIdRest = vals(lens~=2);
    dicKeys = keys(dic);
    pairKeys = dicKeys(lens==2);
    gEdgeFoots = [mod(pairKeys,N) idivide(pairKeys,N)] + idmin;
    restKeys = dicKeys(lens~=2);
    gRestFoots = [mod(restKeys,N) idivide(restKeys,N)] + idmin;
end

function [X,Y] = getRegionBoundary(region, bezierSubdiv)
    numLinSegments = sum(region.gEdgeIsLinear);
    numBezSegments = sum(~region.gEdgeIsLinear);
    X = zeros(numLinSegments + bezierSubdiv * numBezSegments,1);
    Y = zeros(numLinSegments + bezierSubdiv * numBezSegments,1);
    B = region.gpEdge;
    writeId = 1;
    t = linspace(0,1,bezierSubdiv+1); t(end)=[];
    for i = 1:size(region.gEdgeIsLinear,1)
        if region.gEdgeIsLinear(i)
            X(writeId) = region.gpEdge(i,1);
            Y(writeId) = region.gpEdge(i,2);
            writeId = writeId + 1;
        else
            Px = (B(i,1) * (1-t).*(1-t) + B(i,5) * 2*t.*(1-t) + B(i,3) * t.*t)';
            Py = (B(i,2) * (1-t).*(1-t) + B(i,6) * 2*t.*(1-t) + B(i,4) * t.*t)';
            X(writeId + (0:bezierSubdiv-1)) = Px;
            Y(writeId + (0:bezierSubdiv-1)) = Py;
            writeId = writeId + bezierSubdiv;
        end
    end
end

end
end


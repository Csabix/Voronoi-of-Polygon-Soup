function ret = BruteForce(P)
arguments
    P (:,2) double = [0 0; 1 0; .3 1; .4 .3];
end

N = size(P, 1);
if all(P(1,:) == P(end,:))
    N = N-1;
end
indices = cell(1,2*N);
entities = cell(1,2*N);
chars = char('A' + (0:N-1));
names = num2cell(chars);

for i = 1:N
    indices{i} = i;
    entities{i} = {P(i,:)};
end

for i = 1:N
    j = i+1;
    if j == N+1; j = 1; end
    indices{N+i} = {i,j};
    entities{N+i} = {P(i,:), P(j,:)};
    names{N+i} = [chars(i) chars(j)];
end

ret = struct;
ret.tripoints = [];
ret.trinames = {};
triIndices = [];

for i = 1 : 2*N-2
    for j = i+1 : 2*N-1
        for k = j+1 : 2*N
            args = entities([i j k]);
            args = [args{1}(:); args{2}(:); args{3}(:)];
            inds = indices([i j k]);
            name = [names{i} '_' names{j} '_' names{k}];
            skipSegmentCheck = -1; % index of the segment
            % set fun and rearrange args (if necessary)
            if k <= N
                % all points
                fun = @geometry.equidistPPP;
            elseif i > N
                % all lines
                fun = @geometry.equidistLLL;
            elseif j <= N
                % 2 points, 1 line
                % is the first point an endpoint of the segment?
                p1_online = inds{1} == inds{3}{1} || inds{1} == inds{3}{2};
                % is the second point an endpoint of the segment?
                p2_online = inds{2} == inds{3}{1} || inds{2} == inds{3}{2};
                if p1_online && p2_online
                    % the points are the endpoints of the segment
                    continue;
                elseif p1_online && ~p2_online
                    % the first point is an endpoint of the segment
                    fun = @equidistPPL_RP;
                    skipSegmentCheck = k;
                elseif ~p1_online && p2_online
                    % the second point is an endpoint of the segment
                    args = args([2 1 3 4]);
                    fun = @equidistPPL_RP;
                    skipSegmentCheck = k;
                elseif ~p1_online && ~p2_online
                    % general case: both points are off the segment
                    fun = @geometry.equidistPPL;
                else
                    warning('shouldn''t get here');
                end
            else
                % 1 point, 2 lines
                P1isP2 = inds{1} == inds{2}{1};
                P1isP3 = inds{1} == inds{2}{2};
                P1isP4 = inds{1} == inds{3}{1};
                P1isP5 = inds{1} == inds{3}{2};
                summa = P1isP2 + P1isP3 + P1isP4 + P1isP5;
                if summa == 0
                    % general case: the point is off the segments
                    fun = @geometry.equidistPLL;
                elseif summa == 1
                    % the point is the endpoint of one of the segments
                    fun = @equidistPLL_RL;
                    if P1isP2
                        args = args([1 3 4 5]);
                        skipSegmentCheck = j;
                    elseif P1isP3
                        args = args([1 2 4 5]);
                        skipSegmentCheck = j;
                    elseif P1isP4
                        args = args([1 5 2 3]);
                        skipSegmentCheck = k;
                    elseif P1isP5
                        args = args([1 4 2 3]);
                        skipSegmentCheck = k;
                    else
                        warning('shouldn''t get here');
                    end
                else
                    % the segments share an endpoint (and it is the first
                    % parameter)
                    continue;
                end
            end
            outs = cell(1,abs(nargout(fun)));
            % calculate tripoints
            [outs{:}] = fun(args{:});

            for ii = 1:length(outs)
                tripoint = outs{ii};
                % drop point if NaN
                if any(isnan(tripoint)); continue; end
                   
                keep = true;
                % check if the triple point is 'over' every segment
                for partindex = [i j k]
                    partIsPoint = partindex <= N;
                    if partIsPoint || skipSegmentCheck == partindex; continue; end
                    es = entities{partindex};
                    val1 = glm.dot(tripoint - es{1}, es{2} - es{1});
                    val2 = glm.dot(tripoint - es{2}, es{1} - es{2});
                    if val1 < 0 || val2 < 0
                        keep = false;
                        break;
                    end
                end
                if ~keep
                    continue;
                end
                
                % distance of tripoint from the 3 generating parts (it
                % should be the same from all)
                if i <= N
                    dist = glm.distance(tripoint, entities{i}{1});
                else
                    dist = geometry.dist2seg(tripoint, entities{i}{1}, entities{i}{2});
                end
                
                % check if any vertex is closer to the tripoint than the
                % generating parts
                for jj = 1 : N
                    if jj == i || jj == j || jj == k; continue; end
                    dist_to_jj = glm.distance(tripoint, entities{jj}{1});
                    if dist_to_jj < dist
                        keep = false;
                        break;
                    end
                end
                % check if any segment is closer to the tripoint than the
                % generating parts
                for jj = N+1 : 2*N
                    if jj == i || jj == j || jj == k; continue; end
                    [dist2_to_jj,closestIndex] = geometry.distSquared2seg_meta(tripoint, entities{jj}{1}, entities{jj}{2});
                    if closestIndex ~= 0
                        closestPointID = indices{jj}{closestIndex};
                        if closestPointID == i || closestPointID == j || closestPointID == k
                            continue;
                        end
                    end
                    if sqrt(dist2_to_jj) < dist
                        keep = false;
                        break;
                    end
                end
                if keep
                    ret.tripoints = vertcat(ret.tripoints,tripoint);
                    ret.trinames{end+1} = name;
                    triIndices = vertcat(triIndices,[i j k]);
                end
            end
        end
    end
end


% calc vertex tripoints
vertTripoints = zeros(N,3);
vertTripoints(1,:) = [1, N+1, 2*N];
for i = 2:N
    v = [i, N+i-1, N+i];
    vertTripoints(i,:) = v;
end

% calc ideal tripoints
hull = convhull(P(1:N,:));
idealTripoints = [];
for i = 1:length(hull)-1
    ind1 = hull(i);
    ind2 = hull(i+1);
    third = NaN;
    if ind1 + 1 == ind2
        third = N+ind1;
    end
    ideal = [ind1, ind2, third];
    if ind1 == N && ind2 == 1
        ideal = [N, 1, 2*N];
    end
    idealTripoints = vertcat(idealTripoints, ideal);
end

% all tripoints
triPoints = [triIndices; vertTripoints; idealTripoints];
positions = [ret.tripoints; P(1:N,:)]; % ideal points are added later
firstVertex = size(triIndices,1) + 1;
firstIdeal = firstVertex + size(vertTripoints,1);
NT = size(triPoints,1);

connections = [];
for i = 1:NT-1
    a = triPoints(i,:);
    for j = i+1:NT
        b = triPoints(j,:);
        % check if a and b contain 2 common indices
        xx = sort([a b]);
        yy = circshift(xx,1);
        if sum(xx==yy) == 2
            j2 = j;
            if j >= firstIdeal % if b is an ideal tripoint
                vert = positions(i,:);
                normal = -glm.normal(P(b(1),:),P(b(2),:));
                idealPos = vert + 2 * normal;
                positions = vertcat(positions, idealPos);
                j2 = size(positions,1);
            end
            connections = vertcat(connections,[i j2]);
        end
    end
end

ret.linesX = positions(connections',1); % alternates between start and endpoint
ret.linesY = positions(connections',2);

end

% p: common point (on the q0q1 line), b: other point, q0q1: line
function o = equidistPPL_RP(p, b, q0, q1)
    n = glm.normal(q0,q1);
    o = geometry.equidistRP(p,n,b);
end

% p0 is the common point
function [o1,o2] = equidistPLL_RL(p0,p1,q0,q1)
    n = glm.normal(p0,p1);
    [o1,o2] = geometry.equidistRL(p0,n,q0,q1);
end

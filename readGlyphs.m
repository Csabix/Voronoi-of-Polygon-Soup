function glyphs = readGlyphs( file )

fileID = fopen(file);
numGlyps = fread(fileID, 1, 'int32');
if numGlyps == 0; return; end
glyphs(numGlyps) = struct('contours',struct('segments',struct('points',[]))); % prealloc
for g = 1:numGlyps
    numContours = fread(fileID, 1, 'int32');
    if numContours == 0; continue; end
    glyphs(g).contours(numContours) = struct('segments',struct('points',[])); % prealloc
    for c = 1:numContours
        numSegments = fread(fileID, 1, 'int32');
        if numSegments == 0; continue; end
        glyphs(g).contours(c).segments(numSegments) = struct('points',[]); % prealloc
        for s = 1:numSegments
            numPoints = fread(fileID, 1, 'int32');
            points = fread(fileID, numPoints*2, 'float'); % [x0,y0,x1,y1..]
            points = reshape(points, [2 numPoints]) .* [1; -1];
            glyphs(g).contours(c).segments(s).points = points;
        end
    end
end
fclose(fileID);

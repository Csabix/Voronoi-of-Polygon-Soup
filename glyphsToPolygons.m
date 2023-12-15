
function polygons = glyphsToPolygons( glyphs, segments )
arguments
    glyphs,
    segments = 4
end
    bezSettings.segments = segments;
    
    polygons = zeros([0 2]);
    numGlyphs = size(glyphs, 2);
    for i = 1:numGlyphs
        polygons = [polygons; glyphToPolygons(glyphs(i), bezSettings)];
    end
end

function polygons = glyphToPolygons( glyph, bezSettings )
    polygons = zeros([0 2]);
    numContours = size(glyph.contours, 2);
    for i = 1:numContours
        polygons = [polygons; contourToPolygon(glyph.contours(i), bezSettings)];
    end
end

function polygon = contourToPolygon( contour, bezSettings )
    polygon = zeros([0 2]);
    numSegments = size(contour.segments, 2);
    for i = 1:numSegments
        polygon = [polygon; segmentToPoints(contour.segments(i), bezSettings)];
    end
    polygon(end+1,:) = NaN([1 2]);
end

function points = segmentToPoints( segment, bezSettings )
    s = size(segment.points);
    if s(2) == 2
        points = segment.points(:,1)';
    elseif s(2) == 3
        points = bezier2Points(segment.points, bezSettings)';
    end
end

% interpolate bezier2, the last controll point is 
function points = bezier2Points( controllPoints, bezSettings )
    n = bezSettings.segments;
    t = (0:(n-1))/n;
    s = 1 - t;
    points  = controllPoints(:,1) * s.*s + controllPoints(:,2) * 2*s.*t + controllPoints(:,3) * t.*t;
end


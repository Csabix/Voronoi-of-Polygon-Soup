classdef glm        % A namespace for geometric utility functions
methods (Static)    % ___________________________________________
    
function v = length(x)
    % arguments
    %     x (:,2) double
    % end
    %v = sqrt(dot(x,x,2));
    v = sqrt(sum(x.^2,2));
end
function v = distance(x,y)
    % arguments
    %     x (:,2) double; y (:,2) double
    % end
    v = sqrt(sum((x-y).^2,2));
end
function x = normalize(x)
    % arguments
    %     x (:,2) double
    % end
    x = x ./ sqrt(sum(x.^2,2));
end
function v = cosangle(x,y)
    arguments
        x (:,2) double; y (:,2) double
    end
    v = dot(glm.normalize(x),glm.normalize(y),2);
end
function v = dot(x,y)
    % arguments
    %     x (:,2) double; y (:,2) double
    % end
    v = sum(x.*y,2);
end
function v = dot2(x)
    % arguments
    %     x (:,2) double;
    % end
    v = sum(x.^2,2);
end
function v = mix(a,b,x)
    arguments
        a (:,:) double; b (:,:) double; x (:,:) double
    end
    v = a.*(1-x)+b.*x;
end
function I = reflect(I,N)    % same as the GLSL reflect
    arguments
        I (:,2) double; N (:,2) double;
    end
    I = I - 2 * dot(N, I, 2) .* N;
end

function xr = mirror(x,p1,p2)
    % mirrors x onto point p1 or onto line p1p2
    arguments
        x  (:,2) double; p1 (:,2) double; p2 (:,2) double = [];
    end
    if nargin ==3
        n = (p2-p1)*[0 1;-1 0];
        xr = x + 2*n*dot(p1-x,n,2) ./ dot(n,n,2);
    else
        xr = 2*p1-x;
    end
end

function n = normal(p1,p2)
    % normal of a line trough p1,p2. optional third position so that the
    % normal points towards that, eg.: dot(n,x_pos-p1) > 0
    % arguments
    %     p1 (:,2) double; p2 (:,2) double;
    % end
    n = (p2-p1) ./ sqrt(sum((p2-p1).^2,2))*[0 1;-1 0];
end

function n = normalPosSide(p1,p2,pos_side)
    % normal of a line trough p1,p2. optional third position so that the
    % normal points towards that, eg.: dot(n,x_pos-p1) > 0
    % arguments
    %     p1 (:,2) double; p2 (:,2) double; pos_side (:,2) double = [];
    % end
    n = (p2-p1) ./ sqrt(sum((p2-p1).^2,2))*[0 1;-1 0];
    if sum(n.*(pos_side-p1),2) < 0
        n = -n;
    end
end

function b = isparallel(a1,a2,b1,b2)
    % are the two segments paralel a1a2, b1b2?
    % Or are a1,a2,b1 points colinear?
    arguments
        a1 (:,2) double; a2 (:,2) double
        b1 (:,2) double; b2 (:,2) double = a2;
    end
    b = abs(dot(a1-a2,(b1-b2)*[0 1;-1 0],2)) < 1e-10;
end

function [g,w1,w2] = angularbisect(a1,a2,b1,b2)
    % g is the itnersection, w1, w2 are the angular bisector directions
    arguments
        a1 (:,2) double; a2 (:,2) double
        b1 (:,2) double; b2 (:,2) double % = a2 %default value not needed
    end
    if glm.isparallel(a1,a2,b1,b2); error 'bisector of parallel lines'; end
    na = glm.normal(a1,a2); nb = glm.normal(b1,b2);
    g = glm.lineintersect(a1,na,b1,nb);
    w1 = na+nb;    w2 = w1*[0 1;-1 0];
end

function [g,w] = midline(a1,a2,b1,b2)
    % midline of two parallell lines in g+t*w form
    arguments
        a1 (:,2) double; a2 (:,2) double
        b1 (:,2) double; b2 (:,2) double
    end
    if ~glm.isparallel(a1,a2,b1,b2); error 'midline of intersecting lines'; end
    g = 0.25*(a1+a2+b1+b2);
    w = a2-a1; % possible improvement
end


function x = lineintersect(p,n,q,m)
    % intersect two implicit (point, normal) lines
    % assumes that n and m are not parallel
    % arguments
    %     p (:,2) double; n (:,2) double
    %     q (:,2) double; m (:,2) double
    % end
    x = ([n;m]\[n*p'; m*q'])';
end

function [x1,x2] = vquadroots(a,b,c)
    % arguments
    %     a (:,:) double; b (:,:) double; c (:,:) double
    % end
    % Returns roots to the ax^2+bx+c=0 equation, such that x1 <= x2. x2 or both may be NaNs.
    x1 = NaN(size(a));
    x2 = x1;

    d = b.*b-4.*a.*c;
    mask = d>=0;

    sisqrt = sign(a(mask)).*sqrt(d(mask));
    ainv = 0.5./a(mask);
    bd0 = -b(mask);

    x1(mask) = (bd0-sisqrt).*ainv;
    x2(mask) = (bd0+sisqrt).*ainv;
    
    mask = abs(d) < 1e-13;
    x1(mask) = -0.5*b(mask)./a(mask);
    x2(mask) = NaN;

    mask = abs(a) < 1e-13 & abs(b) > 1e-13;
    x1(mask) = -c(mask)./b(mask);
    x2(mask) = NaN;
end

function [x1,x2] = quadroots(a,b,c)
    % arguments
    %     a (:,:) double; b (:,:) double; c (:,:) double
    % end
    % Returns roots to the ax^2+bx+c=0 equation, such that x1 <= x2. x2 or both may be NaNs.
    x1 = NaN; x2 = NaN;
    if abs(a)<1e-13 && abs(b) > 1e-13
        x1 = -c/b;
    else
        d = b*b-4*a*c;
        if abs(d)<1e-13
            x1 = 0.5*b/a;
        elseif d>=0
            x1 = 0.5*(-b-sign(a)*sqrt(d))/a;
            x2 = 0.5*(-b+sign(a)*sqrt(d))/a;
        end
    end
end


end
end
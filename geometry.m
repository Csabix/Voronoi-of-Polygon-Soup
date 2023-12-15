classdef geometry < glm % A namespace for geometric computations
methods (Static)        % ______________________________________    

function bt = bezier2(t,b0,b1,b2)
    % b0,b1,b2 may be vec2-s and t may be a column vector
     arguments
         t  (:,1) double
         b0 (1,2) double; b1 (1,2) double; b2 (1,2) double
     end
    bt = b0.*(1-t).^2 + 2*b1.*t.*(1-t) + b2.*t.^2;
end

function pt = parabola(t,p1,p2,f)
    % Eval parabola with directrix p1p2 and focus f at t column vector
    arguments
        t  (:,1) double
        p1 (1,2) double; p2 (1,2) double; f  (1,2) double
    end
    n = glm.normalPosSide(p2,p1,f);
    lin = glm.mix(p1,p2,t);
    pt =  lin + 0.5/dot(f-p1,n)*dot(lin-f,lin-f,2).*n;
end


function [g,w1,w2] = bisectlines(a1,a2,b1,b2)
    % if parallel same as glm.midline, otherwise same as glm.angularbisect
    % if parallel, w2 contains NaNs
    % arguments (Input)
    %     a1 (:,2) double; a2 (:,2) double
    %     b1 (:,2) double; b2 (:,2) double % = a2 %default value not needed
    % end
    % arguments (Output)
    %     g (:,2) double; w1 (:,2) double; w2 (:,2) double;
    % end
    %va = glm.normalize(a2-a1);
    va = (a2-a1)./sqrt((a2-a1)*(a2-a1)');
    %nb = glm.normalize(b2-b1)*[0 1;-1 0];
    nb = (b2-b1)*[0 1;-1 0]./sqrt((b2-b1)*(b2-b1)');
    if abs(sum(va.*nb,2)) < 1e-7 % dot
        g = 0.25*(a1+a2+b1+b2);
        w1 = va; w2 = [NaN NaN];
    else
        na = va*[0 1; -1 0];
        %g = glm.lineintersect(a1,na,b1,nb);
        g = ([na;nb]\[na*a1'; nb*b1'])';
        w1 = na+nb; w2 = w1*[0 1;-1 0];
    end
end

function o = equidistPPP(a,b,c)
    % arguments
    %     a (1,2) double; % point
    %     b (1,2) double; % point
    %     c (1,2) double; % point
    % end
    n = a-b; m = b-c;
    o = 0.5*[(a+b)*n' (b+c)*m']/[n;m]';
end


function [x1,x2] = equidistPPL(a,b,p,q)
    % point a - point b - line pq
    % arguments
    %     a (1,2) double;                 % point
    %     b (1,2) double;                 % point
    %     p (1,2) double; q (1,2) double; % line
    % end
    %n = glm.normal(p,q);
    n = (q-p)*[0 1;-1 0]./sqrt((q-p)*(q-p)');
    v = (b-a)*[0 1;-1 0];
    p = 0.5*(a+b); % overwrites p, but ok
    vn =  v*n';
    pqn = (p-q)*n';
    A = vn.^2-v*v';
    B = 2*(vn*pqn - v*(p-a)');
    C = pqn.^2 - (p-a)*(p-a)';
    [t1,t2] = glm.quadroots(A,B,C);
    x1 = p + t1*v;
    x2 = p + t2*v;
end

function [x1,x2] = equidistPLL(a,p0,p1,q0,q1)
    % point a - line p0p1 - line q0q1
    % arguments
    %     a (1,2) double;                   % point
    %     p0 (1,2) double; p1 (1,2) double; % line
    %     q0 (1,2) double; q1 (1,2) double; % line
    % end
    n = (p1-p0)*[0 1;-1 0]./sqrt((p1-p0)*(p1-p0)');
    [p,v,w] = geometry.bisectlines(p0,p1,q0,q1);
    pqn = (p-p0)*n';
    vn = v*n';
    A = vn.^2-v*v';
    B = 2*(vn*pqn - v*(p-a)');
    C = pqn.^2 - (p-a)*(p-a)';
    [t1,t2] = glm.quadroots(A,B,C);
    if isnan(t1)
        v = w;
        vn = v*n';
        A = vn.^2-v*v';
        B = 2*(vn*pqn - v*(p-a)');  
        [t1,t2] = glm.quadroots(A,B,C);
    end
    x1 = p + t1*v;
    x2 = p + t2*v;
end

function [g,w1] = angularbisectforLLL(a1,na,b1,nb)
    g = ([na;nb]\[na*a1'; nb*b1'])';
    w1 = na+nb;
end
function o = equidistLLL(p0,p1,q0,q1,r0,r1)
    vp = (p1-p0) ./ sqrt(sum((p1-p0).^2,2)); np = vp*[0 1;-1 0];
    vq = (q1-q0) ./ sqrt(sum((q1-q0).^2,2)); nq = vq*[0 1;-1 0];
    vr = (r1-r0) ./ sqrt(sum((r1-r0).^2,2)); nr = vr*[0 1;-1 0];
    rpq = abs([vp*nq' vq*nr' vr*np']);
    mn = min(rpq); mx = max(rpq);

    if mx < 1e-7
        o = [NaN, NaN];
        return;
    elseif rpq(1) == mn
        [a,va] = geometry.angularbisectforLLL(r0,nr,p0,np);
        [b,vb] = geometry.angularbisectforLLL(r0,nr,q0,nq);
    elseif rpq(2) == mn
        [a,va] = geometry.angularbisectforLLL(p0,np,q0,nq);
        [b,vb] = geometry.angularbisectforLLL(p0,np,r0,nr);
    else %if rpq(3) == mn
        [a,va] = geometry.angularbisectforLLL(q0,nq,r0,nr);
        [b,vb] = geometry.angularbisectforLLL(q0,nq,p0,np);
    end
    na = va*[0 1; -1 0];
    if abs(na*vb') < 1e-10
        o = [NaN NaN];
    else
        nb = vb*[0 1; -1 0];
        o = ([na;nb]\[na*a'; nb*b'])';
    end
end

function b = pointIsOverEdge(newTPos, ev1, ev2)
    b = (newTPos - ev1)*(ev2 - ev1)' >= 0 && ...
        (newTPos - ev2)*(ev1 - ev2)' >= 0 && ...
        ~any(isnan(newTPos));
end

function o = equidistRP(p,v,a)
    % arguments
    %     p (1,2) double; v (1,2) double;  % ray (edge normal at its end)
    %     a (1,2) double;                  % point
    % end
    s = -2*glm.dot(p-a,v);
    if abs(s) < 1e-10
        o = [NaN NaN];
    else
        o = p + v * glm.dot2(p-a)/(s);
    end
end

function [x1,x2] = equidistRL(p,v,a,b)
    % arguments
    %     p (1,2) double; v (1,2) double; % ray (edge normal at its end)
    %     a (1,2) double; b (1,2) double; % line
    % end
    %n = glm.normal(a,b);
    n = (b-a)*[0 1;-1 0]./sqrt((b-a)*(b-a)');
    vn  = v*n';
    pan = (p-a)*n';
    A = vn.^2- v*v';
    B = 2*vn.*pan;
    C = pan.^2;
    [t1,t2] = glm.quadroots(A,B,C);
    x1 = p + v*t1;
    x2 = p + v*t2;
end

function d = dist2seg(p,a,b)
    d = glm.distance(p,a+(b-a).*min(max(glm.dot(p-a,b-a)./glm.dot2(b-a),0),1));
end

function d2 = distSquared2seg(p,a,b)
    d2 = sum((-p + a+(b-a).*min(max(sum((p-a).*(b-a),2)./sum((b-a).^2,2),0),1)).^2,2);
end

% distance squared from point p to segment ab
% also returns whether p is closest to a or b or the middle of the segment
% (e = 1, 2 and 0 respectively)
function [d2,e] = distSquared2seg_meta(p,a,b)
    t = min(max(sum((p-a).*(b-a),2) ./ sum((b-a).^2,2) ,0),1);
    e = 1*(t==0) + 2*(t==1);
    d2 = sum((a+(b-a).*t-p) .^2,2);
end

function t = paramOnLine(p,a,b)
    arguments
        p (:,2) double;
        a (1,2) double; b (1,2) double;
    end
    t = glm.dot(p-a,b-a)./glm.dot2(b-a);
end

end
end
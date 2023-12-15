addpath Geomatplot\

%%
clf; disp 'test02/01 polygon triple points'
getPolyVertex; % clears persistent vars

vals = [0.4990 0.5811; -0.0358 0.5847; 0.1852 0.8773; -0.2525 1.0110; -0.0915 -0.1292; 1.2522 0.4222; 0.4092 0.9141; 0.8473 0.5839; 0.3721 0.3073];
%vals = [0 0; 1 0; .3 1; .4 .3];
poly = Polygon(vals,'g','FaceAlpha',0.1,'LineWidth',1);
for i = 1:10
    Point(poly,@(p) getPolyVertex(p, i),[0 .8 0]);
end

% BruteForce returns a struct with multiple fields
cval = CustomValue(poly,@BruteForce);

seq = PointSequence(cval,@(v) v.tripoints,'y',5);

% for i = 1:25
%     center = Point(seq,@(s) s(i,:),'LabelV','off','Visible','off');
%     Circle(center,poly);
%     Text(center,poly,@(p) callbfun2(p,i),'Interpreter','none','Col','g');
% end

c1 = Point('c1',[-.5 -.25],'b',5); c2 = Point('c2',[1.5 1.25],'b',5);
Image(poly,@regioncolor,c1,c2,'Res',1024);
xlim([-.5 1.5]); ylim([-.25 1.25]);

SegmentSequence(cval,@(v) [v.linesX(:) v.linesY(:)],'c');

%SegmentSequence(poly,@delaunaycallb,'m',3);

% g = poly.parent;
% s = struct('Visible','on','Color','c','LineWidth',2);
% lns = dcurve(g,'circs',{poly},@delaunaycircs,s);
% global index;

%%

function [x,y] = delaunaycallb(poly)
    es = edges(delaunayTriangulation(poly(1:end-1,:)));
    x = [poly(es(:,1),1) poly(es(:,2),1)];
    y = [poly(es(:,1),2) poly(es(:,2),2)];
end

function [x,y] = delaunaycircs(t,poly)
    global index;
    if isempty(index); index = 1; end
    poly = poly.value;
    dT = delaunayTriangulation(poly(1:end-1,:));
    [C,r] = circumcenter(dT,index);
    t = 2*pi*t';
    x = C(:,1) + r.*cos(t);
    y = C(:,2) + r.*sin(t);
    x = [x NaN(length(r),1)]';
    y = [y NaN(length(r),1)]';
    x = x(:);    y = y(:);
end

% drop the last duplicated vertex 
function o = getPolyVertex(poly, i)
    persistent n;
    if nargin == 0; n = []; o=[]; return; end
    if isempty(n) || n~=size(poly,1)-1
        n = size(poly, 1)-1;
        clim([1 2*n]);
        col = net(haltonset(3),2*n);
        col(1:n,:)     = col(1:n,:)    .*[.5 .5 .5] + [.5 .1 .0];
        col(n+1:end,:) = col(n+1:end,:).*[.5 .5 .5] + [.0 .1 .5];
        colormap(col*.5+.5);
        colorbar('Ticks',linspace(1.5,2*n-0.5,2*n),'TickLabels',1:2*n);
    end
    if n+1 == i; throw('intentional');end
    o = poly(i,:);
end

function inds = regioncolor(x,y,poly)
    dist = inf(size(x));
    inds = zeros(size(x));
    N = size(poly,1)-1;

    abx = diff(poly(:,1))'; aby = diff(poly(:,2))';
    abilen = 1./sqrt(abx.*abx + aby.*aby);
    nx =  aby.*abilen;      ny = -abx.*abilen;
    
    for i = 1:N
        ax = poly(i+0,1)-x; ay = poly(i+0,2)-y;
        dvert = ax.*ax + ay.*ay;
        inds(dvert < dist) = i;
        dist = min(dist,dvert);

        bx = poly(i+1,1)-x; by = poly(i+1,2)-y;
        mask = ax.*abx(i) + ay.*aby(i) < 0 & bx.*abx(i) + by.*aby(i) > 0;
        dedge = (nx(i).*ax+ny(i).*ay).^2;
        mask = mask & dist > dedge;
        inds(mask) = i+N;
        dist(mask) = dedge(mask);
    end
end
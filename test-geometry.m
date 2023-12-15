addpath Geomatplot\
%%
clf; disp 'test01/01 PPL'
a = Point('a',[0.6 .55]);
b = Point('b',[1 .5]);
p = Point('p',[0.0 0.0]);
q = Point('q',[1.0 0.0]);
Segment(p,q,'b',2)
PerpendicularBisector(a,b)
Curve(p,q,a,@geometry.parabola)
Curve(p,q,b,@geometry.parabola)
%E = Point;
%F = Point;
%poly = Polygon(A,B,C,'b--');
%vals = poly.value;
%cntr = Contour(vals)
x1 = Point('x1',a,b,p,q,@geometry.equidistPPL);
Circle(x1, a);

%%
clf; disp 'test01/02 PLL'
a  = Point('f' ,[0.6 .55]);
p0 = Point('p0',[1.0 .5]);
p1 = Point('p1',[0.0 0.0]);
q0 = Point('q0',[0.0 1.0]);
q1 = Point('q1',[1.0 1.0]);
Segment(p0,p1,'b',2)
Segment(q0,q1,'b',2)
Curve(p0,p1,a,@geometry.parabola)
Curve(q0,q1,a,@geometry.parabola)
x1 = Point('x1',a,p0,p1,q0,q1,@geometry.equidistPLL);
Circle(x1,a);

%%
clf; disp 'test01/03 PPP'
a = Point('a',[0 0]);
b = Point('b',[1 0]);
c = Point('c',[.5 1]);
x = Point('x', a,b,c, @geometry.equidistPPP);
Circle(x,a);

%%
clf; disp 'test01/04 LLL'
a1 = Point('a1',[0 0]);
a2 = Point('a2',[.5 0]);
b1 = Point('b1',[1 .5]);
b2 = Point('b2',[1 1]);
c1 = Point('c1',[.5 1]);
c2 = Point('c2',[0 .5]);
l1 = Line(a1,a2);
Line(b1,b2);
Line(c1,c2);
x = Point('x', a1,a2,b1,b2,c1,c2, @geometry.equidistLLL);
Circle(x,l1);
%%
clf; disp 'test01/04 RP'
a = Point('a',[0.6 .55]);
p = Point('p',[0.0 0.0]);
q = Point('q',[1.0 0.0]);
%Line(a,b);
Line(p,q);
Curve(p,q,a,@geometry.parabola)
pq = Point(p,q,@(a,b) (b-a)*[0 1;-1 0],[1,1,1],0.1,'LabelV','off');
nulla = Point(@() [0 0],[1,1,1],0.1,'LabelV','off');
v = Eval((pq-nulla)/Distance(p,q));
o = Point('o',p,v,a,@geometry.equidistRP);
Circle(o,a);
%%
clf; disp 'test01/05 RL'
a = Point('a',[0 1]);
b = Point('b',[0 0]);
p = Point('p',[0.5 0.5]);
q = Point('q',[1.0 0.5]);
Line(a,b); Line(p,q);
Curve(a,b,p,@geometry.parabola)
pq = Point(p,q,@(a,b) (b-a)*[0 1;-1 0],[1,1,1],0.1,'LabelV','off');
nulla = Point(@() [0 0],[1,1,1],0.1,'LabelV','off');
v = Eval((pq-nulla)/Distance(p,q));
o = Point('o',p,v,a,b,glm.selectout(@geometry.equidistRL,2));
Circle(o,p);

%%
clf; disp 'test dist2seg'
A = Point;
B = Point;
C = Point;
Segment(B,C);
Circle(A,Scalar(A,B,C,@geometry.dist2seg));

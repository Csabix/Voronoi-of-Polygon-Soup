%%
disp 'Performance measurement test - bruteforce - rand radial polygons (timeit)'
% bruteforce-times.mat

M = 2;     % M different random
ns = 3:1:30; % n-gons

rng(0);
test0_times = zeros(size(ns,2),M+1);
for i = 1:size(ns,2)
    n = ns(i);
    fprintf('%3d-gons:   ',n);
    test0_times(i,1) = n;
    for j = 1:M
        for k = 1:5
            poly = Test.randRadPoly(n);
            
            doAll = @() BruteForce(poly);
            try
                test0_times(i,j+1) = timeit(doAll);
                break;
            catch
            end
        end
    end
    fprintf('%f s\n',mean(test0_times(i,2:M+1)));
end

%%
disp 'Performance measurement test - regular N-gons (timeit)'
% regular-n-gon-times.mat

ns = 3:30; % n-gons
test1_times = zeros(size(ns,2),2);
for i = 1:size(ns,2)
    n = ns(i);
    fprintf('%3d-gon:   ',n);
    test1_times(i,1) = n;

    poly = Test.regularNGon(n);
    
    doAll = @() TripointGraph(ConnectionData(poly));
    test1_times(i,2) = timeit(doAll);
    fprintf('%f s\n',test1_times(i,2));
end

%%
disp 'Performance measurement test - rand radial polygons (timeit)'
% rand-radial-poly-times.mat

M = 7;     % M different random
ns = 5:5:50; % n-gons

rng(0);
test2_times = zeros(size(ns,2),M+1);
for i = 1:size(ns,2)
    n = ns(i);
    fprintf('%3d-gons:   ',n);
    test2_times(i,1) = n;
    k = 5;
    for j = 1:M
        poly = Test.randRadPoly(n);
        
        doAll = @() TripointGraph(ConnectionData(poly));
       try
          test2_times(i,j+1) = timeit(doAll);
       catch
          if k > 0
            j = j-1;
            k = k-1;
          end
       end
    end
    fprintf('%f s\n',mean(test2_times(i,2:M+1)));
end

%%
disp 'Performance measurement test - glyphs (timeit)'
% glyph-times.mat

ns = 1:5; % Bezier-subdivisions

a = dir('glyphs');
a = a(arrayfun(@(b)length(b.name),a)==7);
test4_times = zeros(size(a,1), 2*size(ns,2));
test4_names = arrayfun(@(b) string(b.name),a);

for i = 1:length(a)
    b = a(i);
    fprintf('%s\n', b.name);

    glyph = readGlyphs("glyphs/"+b.name);
    for j = 1:size(ns,2)
        n = ns(j);
        poly = glyphsToPolygons(glyph, n);
        try
            cd = ConnectionData(poly);
            test4_times(i,2*j-1) = size(cd.pVertices,1);
            fprintf('%3d   ', test4_times(i,2*j-1));
            doAll = @() TripointGraph(ConnectionData(poly));
            test4_times(i,2*j) = timeit(doAll);
            fprintf('%f s\n', test4_times(i,2*j));
        catch
        end
    end
end

%%
disp 'Performance measurement test - rotated square grid (timeit)'
% rotated-grid-times.mat

M = 3;     % M different random orientations
ns = 1:7; % (1,1),(1,2)..(1,7),(2,2),(2,3)..(2,7),(3,3)..

num = size(ns,2)*(size(ns,2)+1)/2;

test5_times = zeros(num,M+1);
rng(0);
row = 1;
for i=1:size(ns,2)
    n = ns(i);
    for j=i:size(ns,2)
        m = ns(j);
        fprintf('(%2d,%2d):   ',n,m);
        test5_times(row,1) = n*m*4;
        for k = 1:M
            poly = Test.concat_polys(Test.rotated_square_grid(n,m));
            
            doAll = @() TripointGraph(ConnectionData(poly));
            test5_times(row,k+1) = timeit(doAll);
        end
        fprintf('%f s\n',mean(test5_times(row,2:M+1)));
        row = row+1;
    end
end

%%

% create plot with all test results

labels = {"Regular polygons","Random radial polygons","Glyphs","Rotated squares","Brute force"};
styles = {'c+', 'gx', 'ob', '*r','ms'};
data = cell(1,5);
data{1} = test1_times;

N2 = size(test2_times,2) - 1;
data{2} = [repmat(test2_times(:,1),[N2 1]), reshape(test2_times(:,2:N2+1),[],1)];
data{2} = data{2}(data{2}(:,2) ~= 0, :);

data{3} = [reshape(test4_times(:,1:2:5),[],1), reshape(test4_times(:,2:2:6),[],1)];

data{4} = [repmat(test5_times(:,1),[3 1]), reshape(test5_times(:,2:4),[],1)];

N5 = size(test0_times,2) - 1;
data{5} = [repmat(test0_times(:,1),[N5 1]), reshape(test0_times(:,2:N5+1),[],1)];
data{5} = data{5}(data{5}(:,2) ~= 0, :);

clf;
indices = [5,4,3,2,1];
for ii = 1:length(indices)
    i = indices(ii);
    plot(data{i}(:,1), data{i}(:,2), styles{i});
    if ii == 1
         hold on;
    end
end

legend(labels(indices),'Location','northeast');
xlabel('Number of vertices');
ylabel('Time (s)');
hold off;


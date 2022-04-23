close all;
s = load('4_2.mat');
tris = s.triangles + 1;
figure()
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunction','edgecolor','none','facecolor','interp');
hold on;

all_mins = s.all_mins' + 1;
min_vtxs = s.vertices(all_mins, :);

select_mins = s.select_mins' + 1;
smin_vtxs = s.vertices(select_mins, :);

all_maxes = s.all_maxes' + 1;
max_vtxs = s.vertices(all_mins, :);

select_maxes = s.select_maxes' + 1;
smax_vtxs = s.vertices(select_maxes, :);

all_saddles = s.all_saddles' + 1;
saddle_vtxs = s.vertices(all_mins, :);

select_saddles = s.select_saddles' + 1;
ssaddle_vtxs = s.vertices(select_saddles, :);

scatter3(min_vtxs(:,1),min_vtxs(:,2),min_vtxs(:,3), 'bo');
scatter3(max_vtxs(:,1),max_vtxs(:,2),max_vtxs(:,3), 'ro');
scatter3(saddle_vtxs(:,1),saddle_vtxs(:,2),saddle_vtxs(:,3), 'go');

scatter3(smin_vtxs(:,1),smin_vtxs(:,2),smin_vtxs(:,3), 'bo','filled');
scatter3(smax_vtxs(:,1),smax_vtxs(:,2),smax_vtxs(:,3), 'ro','filled');
scatter3(ssaddle_vtxs(:,1),ssaddle_vtxs(:,2),ssaddle_vtxs(:,3), 'go','filled');

figure()
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.colors,'edgecolor','none','facecolor','flat');
hold on;




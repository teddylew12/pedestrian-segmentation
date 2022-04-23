s = load('4_2.mat');
tris = s.triangles + 1;
figure()
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.colors,'edgecolor','none','facecolor','flat');
hold on;

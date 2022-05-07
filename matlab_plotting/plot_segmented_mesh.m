s = load('../off/base_mesh.mat');
tris = s.triangles + 1;
figure()
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.colors,'edgecolor','none','facecolor','flat');
hold on;


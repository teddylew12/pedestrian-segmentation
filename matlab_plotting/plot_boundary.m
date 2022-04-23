s = load('eigenfuncs.mat');
tris = s.triangles + 1;
figure()
boundary = s.colors;
boundary(boundary ~=1)=0;
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',boundary,'edgecolor','none','facecolor','flat');
title("Level Set Triangles")
view(90,0);
close all;
s = load('eigenfuncs.mat');
tris = s.triangles + 1;
figure()
subplot(2,2,1)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
title("First Eigenfunction")
subplot(2,2,2)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
title("Second Eigenfunction")
subplot(2,2,3)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
title("Third Eigenfunction")
subplot(2,2,4)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
title("Fourth Eigenfunction")
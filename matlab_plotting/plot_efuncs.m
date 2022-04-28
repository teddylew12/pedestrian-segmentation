
s = load('../runs/4_2.mat');
tris = s.triangles + 1;
figure()
subplot(2,2,1)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
title("First Eigenfunction")
subplot(2,2,2)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
title("Second Eigenfunction")
subplot(2,2,3)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
title("Third Eigenfunction")
subplot(2,2,4)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
title("Fourth Eigenfunction")


s1 = load('../runs/matching2.mat');
tris1 = s1.triangles + 1;
figure()
subplot(2,2,1)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
title("First Eigenfunction")
subplot(2,2,2)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
title("Second Eigenfunction")
subplot(2,2,3)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
title("Third Eigenfunction")
subplot(2,2,4)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata',s1.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
title("Fourth Eigenfunction")
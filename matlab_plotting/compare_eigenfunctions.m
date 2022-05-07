s = load('../off/base_mesh.mat');
tris = s.triangles + 1;
figure()
subplot(2,4,1)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
view(-180,-90)
title("First Eigenfunction")
subplot(2,4,2)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
view(-180,-90)
title("Second Eigenfunction")
subplot(2,4,3)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
view(-180,-90)
title("Third Eigenfunction")
subplot(2,4,4)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
view(-180,-90)
title("Fourth Eigenfunction")


s1 = load('../runs/first_run5/initial_eigenfunctions.mat');
tris1 = s1.triangles + 1;
subplot(2,4,5)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
view(-180,-90)
subplot(2,4,6)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
view(-180,-90)
subplot(2,4,7)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
view(-180,-90)
subplot(2,4,8)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata',s1.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
view(-180,-90)




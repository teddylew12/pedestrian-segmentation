
s = load('../runs/standing1.mat');
tris = s.triangles + 1;
figure()
subplot(2,2,1)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
view(0,225)
title("First Eigenfunction")
subplot(2,2,2)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
view(0,225)
title("Second Eigenfunction")
subplot(2,2,3)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata', s.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
view(0,225)
title("Third Eigenfunction")
subplot(2,2,4)
patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
view(0,225)
title("Fourth Eigenfunction")


s1 = load('../runs/running1_initial.mat');
tris1 = s1.triangles + 1;
figure()
subplot(2,2,1)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,2),'edgecolor','none','facecolor','interp');
view(0,225)
title("First Eigenfunction")
subplot(2,2,2)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,3),'edgecolor','none','facecolor','interp');
view(0,225)
title("Second Eigenfunction")
subplot(2,2,3)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata', s1.eigenfunctions(:,4),'edgecolor','none','facecolor','interp');
view(0,225)
title("Third Eigenfunction")
subplot(2,2,4)
patch('faces',tris1,'vertices',s1.vertices,'facevertexcdata',s1.eigenfunctions(:,5),'edgecolor','none','facecolor','interp');
view(0,225)
title("Fourth Eigenfunction")



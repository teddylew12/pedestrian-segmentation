s = load('allcolors2.mat');
tris = s.triangles + 1;
for i = 1:size(s.allcolors,2)
p=patch('faces',tris,'vertices',s.vertices,'facevertexcdata',s.allcolors(:,i),'edgecolor','none','facecolor','flat');
saveas(p,sprintf("animate/%d.png",i))
end

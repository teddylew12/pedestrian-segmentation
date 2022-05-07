s = load('../off/base_mesh.mat');
tris = s.triangles + 1;
figure()
view(-180,90)
subplot(1,2,1)
patch('faces',tris,'vertices',s.vertices,'edgecolor','red','facecolor','none');
hold on;

xyz = s.landmark_matrix(:,5:7);
for i=1:5
    txt = sprintf("%d",i);
    text(xyz(i,1),xyz(i,2),xyz(i,3),txt)
end

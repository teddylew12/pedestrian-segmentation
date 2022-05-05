s = load('../runs/standing1.mat');
tris = s.triangles + 1;
figure()
patch('faces',tris,'vertices',s.vertices,'edgecolor','red','facecolor','none');
hold on;

xyz = double(s.xyz);
for i=1:5
    txt = sprintf("%d",i);
    text(xyz(i,1),xyz(i,2),xyz(i,3),txt)
end

s1 = load('../runs/sitting1.mat');
tris1 = s1.triangles + 1;
figure()
patch('faces',tris1,'vertices',s1.vertices,'edgecolor','red','facecolor','none');
hold on;

xyz1 = double(s1.xyz);
for i=1:5
    txt = sprintf("%d",i);
    text(xyz1(i,1),xyz1(i,2),xyz1(i,3),txt)
end

s2 = load('../runs/running1.mat');
tris = s2.triangles + 1;
figure()
patch('faces',tris,'vertices',s2.vertices,'edgecolor','red','facecolor','none');
hold on;

xyz2 = double(s2.xyz);
for i=1:5
    txt = sprintf("%d",i);
    text(xyz2(i,1),xyz2(i,2),xyz2(i,3),txt)
end

s2 = load('../runs/lying1.mat');
tris = s2.triangles + 1;
figure()
patch('faces',tris,'vertices',s2.vertices,'edgecolor','red','facecolor','none');
hold on;

xyz2 = double(s2.xyz);
for i=1:5
    txt = sprintf("%d",i);
    text(xyz2(i,1),xyz2(i,2),xyz2(i,3),txt)
end
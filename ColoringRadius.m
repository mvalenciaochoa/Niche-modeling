function classifiers=ColoringRadius(in,show,outlier,outlier2,alpha)
if nargin <3
    outlier=false;
    outlier2=false;
end
out1=[];
out2=[];
T2 = in{1}.T2;
Z = in{1}.Z;
R = in{1}.R;

reps = size(Z);
caps = reps(3);
template = Z(:, :, 1);
data = NaN(length(template(:)), caps);

for i = 1 : caps
    template = Z(:, :, i);
    data(:, i) = template(:);
end

nanDetector = sum(data, 2);
pointer = ~isnan(nanDetector);
idx = find(pointer==1);


points = T2{:,4:end};
normalizers=[max(points(:,:));min(points(:,:))];
points(:,:)=(points(:,:)-normalizers(2,:))./(normalizers(1,:)-normalizers(2,:));
data(idx,:) =  (data(idx,:)-normalizers(2,:))./(normalizers(1,:)-normalizers(2,:));

if outlier
    [~,~,RD,chi_crt]=DetectMultVarOutliers(points(:,:));
    id_out=RD>chi_crt(4);
    out1=points(id_out,:);
    points=points(~id_out,:);
end

[coeff,~,~,~,explained]=pca(points(:,:));
pin=points(:,:)*coeff(:,1:3);

if ~isempty(out1)
    out1=out1*coeff(:,1:3);
end


if outlier2
    %siz=round(size(pin,1)*0.3);
    [~,~,RD,chi_crt]=DetectMultVarOutliers(pin);
    id_out=RD>chi_crt(4);
    out2=pin(id_out,:);
    pin=pin(~id_out,:);
end

nodes = boundary(pin(:,1),pin(:,2),pin(:,3),alpha);

boundPointsIndex = unique(nodes)';
boundPoints = points(boundPointsIndex,:);
pointsSize = length(boundPointsIndex);
samples = length(points);
radius = zeros(pointsSize,samples);
map = ones(reps(1), reps(2));

for i=1:pointsSize
    for j=1:samples
        radius(i,j)=norm(points(boundPointsIndex(i),:)-points(j,:));
    end
end

radiusClass = min(radius)
radiusIndex = find(radiusClass>0)
radiusClass = radiusClass(:,setdiff(1:end,boundPointsIndex))

response = NaN(1,length(radiusClass));
intensity = NaN(1,length(idx));

for i=1:length(idx)
    for j=1:length(radiusClass)
        response(j) = norm(points(radiusIndex(j),:)-data(idx(i),:));
    end
    intensity(i) = sum(response<=radiusClass);
end

intensity = (intensity - min(intensity))./(max(intensity)-min(intensity));

final = NaN(length(template(:)),1);

final(idx)=intensity;

map(:) = final(:);

classifiers.nodes = boundPoints;
classifiers.index = radiusIndex;
classifiers.radius = radiusClass;
classifiers.normalizers = normalizers;
classifiers.T2 = T2;
classifiers.map = map;

outT=[out1;out2];

if show
    trisurf(nodes,pin(:,1),pin(:,2),pin(:,3), 'Facecolor','cyan','FaceAlpha',0.8); axis equal;
    hold on
    plot3(pin(:,1),pin(:,2),pin(:,3),'.r')
    hold off
    

    figure
    clf
    geoshow(map, R, 'DisplayType','surface');
    contourcmap('jet',0:0.05:1, 'colorbar', 'on', 'location', 'vertical')
end
    

if isempty(outT)
    grap=false;
else
    grap=true;
end
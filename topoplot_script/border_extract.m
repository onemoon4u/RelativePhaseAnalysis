function borderCoords = border_extract(filename, resolutionfactor)
%border_extract(filename)
%   Description:    Extracts coordinates for a border.  Function works best
%                   if the image a simple black and white bitmap of the
%                   border.  Resolution factor can be used to reduce the
%                   amount of coordinates generated, the higher the factor,
%                   the fewer the points generated. (ie resolutionfactor of
%                   4 grabs every other pixel as a point) 


%import the figure
I = imread(filename);
%iterate through and find points that are adjacent to a white space
borderCoords = [];
borderTemp = zeros(2,1);
%flag for reducing resolution
flag = 0;
size(I,1)
size(I,2)
for i = 2:size(I,1) - 1
    for j = 4:size(I,2) - 1
        if(I(i,j,1) == 255)
            continue;
        end
        
        if((I(i, j+1) == 255 ) && (I(i,j-1) == 0 || I(i,j-2) == 0 || I(i,j-3) == 0 || I(i+1,j) == 0 || I(i-1,j) == 0))
            if(flag < resolutionfactor)
                flag = flag + 1;
            else
                borderTemp(1,end + 1) = i;
                borderTemp(2, end) = j;
                flag = 0;
            end
        end
        
    end
end

% reorder the coodinates by nearest distance
borderCoords(:,end+1) = borderTemp(:,2);
borderTemp(:,1) = [];

minDistance = 100000;
index = 0;
ending = size(borderTemp,2);
for i = 1:ending
    x1 = borderCoords(1,i);
    y1 = borderCoords(2,i);
    for j = 1:size(borderTemp,2);
        dist = sqrt((x1 - borderTemp(1,j))^2 + (y1 - borderTemp(2,j))^2);
        if(dist < minDistance)
            minDistance = dist;
            index = j;
        end
    end
    
    if(minDistance < 1000)
        borderCoords(:,end+1) = borderTemp(:,index);
        borderTemp(:,index) = [];
    else
        break;
    end
    
    minDistance = 100000;
end

%close the loop
borderCoords(:,end+1) = borderCoords(:,1);

%flip on the x axis
borderCoords(1,:)  = -borderCoords(1,:);
scatter(borderCoords(2,:), borderCoords(1,:), '.');
plot3(borderCoords(1,:), borderCoords(2,:), ones(size(borderCoords,2)));

end
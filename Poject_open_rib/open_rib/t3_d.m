function [im3_t] = t3_d(input)
 current_data = input;
 data = smooth3(current_data,'box',1);
 p1 = patch(isosurface(data,.5),...
      'Facecolor', 'yellow','Edgecolor','green');
 p2 = patch(isocaps(data,.5),...
      'Facecolor','interp','Edgecolor','red');
 isonormals(data,p1)
 view(3);camlight; lighting phong
 im3_t = data;
end


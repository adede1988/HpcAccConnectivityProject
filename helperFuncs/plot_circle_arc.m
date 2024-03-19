function [center_x, center_y, start_angle, end_angle] = plot_circle_arc(x1, y1, x2, y2, radius)
    % Step 1: Calculate midpoint between the two given points
    mid_x = (x1 + x2) / 2;
    mid_y = (y1 + y2) / 2;
    
    % Step 2: Calculate the slope of the line passing through the two given points
    slope = (y2 - y1) / (x2 - x1);
    
    % Step 3: Calculate the negative reciprocal of the slope to find the slope of the perpendicular bisector
    perpendicular_slope = -1 / slope;
    
    % Step 4: Calculate the y-intercept of the perpendicular bisector passing through the midpoint
    b = mid_y - perpendicular_slope * mid_x;
    
    % Step 5: Calculate the x-coordinate of the center of the circle
    center_x = b / (1 + perpendicular_slope^2);
    
    % Step 6: Calculate the y-coordinate of the center of the circle
    center_y = perpendicular_slope * center_x + b;
    
    % Calculate the distance from the center to one of the given points
    distance = sqrt((center_x - x1)^2 + (center_y - y1)^2);
    
    % Scale the coordinates of the center using the radius
    scale_factor = radius / distance;
    center_x = center_x * scale_factor;
    center_y = center_y * scale_factor;
    
    % Calculate the start and end angles of the arc
    start_angle = atan2(y1 - center_y, x1 - center_x);
    end_angle = atan2(y2 - center_y, x2 - center_x);
    
    % Adjust angles to ensure the arc passes through the given points
    if start_angle > end_angle
        end_angle = end_angle + 2*pi;
    end
    
    % Display the center coordinates and angles
%     disp(['Center coordinates: (', num2str(center_x), ',', num2str(center_y), ')']);
%     disp(['Start angle: ', num2str(start_angle), ' radians']);
%     disp(['End angle: ', num2str(end_angle), ' radians']);
    
    % Plot the two points
    plot([x1, x2], [y1, y2], 'ro');  % Red circles
    
    % Plot the circle arc
    theta = linspace(start_angle, end_angle, 100);
    arc_x = center_x + radius * cos(theta);
    arc_y = center_y + radius * sin(theta);
    plot(arc_x, arc_y, 'b-');  % Blue line
    
%     % Set plot aspect ratio to equal
%     axis equal;
%     
%     % Label the plot
%     xlabel('X');
%     ylabel('Y');
%     title('Circle Arc between Two Points');
%     
%     % Mark the center point
%     hold on;
    plot(center_x, center_y, 'gx', 'MarkerSize', 10);  % Green cross
%     hold off;
end

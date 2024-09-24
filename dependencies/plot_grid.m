function plot_grid(layer_type, sampling_points)

figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [700 100 500 500]);

if strcmp(layer_type, 'double')

    sampling_points_inner = sampling_points(:, 1:end/2);
    sampling_points_outer = sampling_points(:, end/2+1:end);

    plot3(sampling_points_outer(1, :), sampling_points_outer(2, :), sampling_points_outer(3, :), '.');
    hold on;
    plot3(sampling_points_inner(1, :), sampling_points_inner(2, :), sampling_points_inner(3, :), '.');
    hold off;

else
    plot3(sampling_points(1, :), sampling_points(2, :), sampling_points(3, :), '.');
end

grid on;
box on;
axis equal;

% find farthest sampling point
r = vecnorm(sampling_points, 2, 1);
axis([-1 1 -1 1 -1 1] * max(r)*1.2);

xlabel('x');
ylabel('y');
zlabel('z');

camproj('perspective');

end
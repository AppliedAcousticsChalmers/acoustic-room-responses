function [sdm_data] = quantize_sdm_angles(sdm_data, no_of_incidence_angles)

% get evenly distributed indicence angles (the actual number may differ slightly from no_of_incidence_angles)
points_target = points_on_sphere(no_of_incidence_angles);
points_target = points_target.';

fprintf('Quantizing the SDM data to %d directions.\n\n', size(points_target, 1));

% go to Cartesian coordinates
[points_sdm(:, 1), points_sdm(:, 2), points_sdm(:, 3)] = sph2cart(sdm_data.azi_rad, sdm_data.ele_rad, 1);

idx = zeros(size(points_sdm(:, 1)));

% find the closes target point for each point in the SDM data
for ii = 1 : size(points_sdm, 1)
    
    dist = sum((points_target - repmat(points_sdm(ii, :), size(points_target, 1), 1)).^2, 2);

    [~, idx(ii)] = min(dist);

end

% get the quantized SDM data in Cartesian coordinates
points_sdm_quantized = points_target(idx, :);

% figure;
% plot3(points_sdm_quantized(:, 1), points_sdm_quantized(:, 2), points_sdm_quantized(:, 3), 'r.');
% grid on;
% axis equal;

% convert back to SDM format
[sdm_data.azi_rad, sdm_data.ele_rad, ~] = cart2sph(points_sdm_quantized(:, 1), points_sdm_quantized(:, 2), points_sdm_quantized(:, 3));

end


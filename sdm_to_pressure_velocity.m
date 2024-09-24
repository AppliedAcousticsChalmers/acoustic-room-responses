% This script loads an SDM representation of an acoustic room response and
% converts it into sound pressure and velocity on a spatial grid by 
% imposing the digital samples of the pressure room response onto plane 
% waves and computing the sound pressure and velocity due to the set of 
% plane waves.
% 
% This script is inspired by the Chalmers Auralization Toolbox: 
% https://github.com/AppliedAcousticsChalmers/auralization-toolbox
%
% It considers three cases: 
% - The pressure is sought on a volumetric grid ('cubical_volume').
% - The pressure is sought on a double-layer surface that is either cubical or spherical.
% - Pressure and normal particle velocity are sought on a single-layer surface that is either cubical or spherical.
% 
% Note that the variables 'pressure' and 'velocity' have taps_pw/2 leading 
% and taps_pw/2 trailing zeros compared to 'sdm_data.p'.

clear;

addpath('dependencies/');

% SDM data file
%sdm_file = 'room_data/sdm_big_hall.mat';
sdm_file = 'room_data/sdm_listening_lab.mat';

% sampling grid
%grid_shape = 'cubical_volume'; L_cube = 5;
grid_shape = 'spherical_surface'; L_sphere = 25; % must be a square number
%grid_shape = 'cubical_surface'; L_cube = 5;

% L_cube is the number of sampling points per dimension (if cubical grid)
% L_sphere is the total number of sampling points (if spherical grid);    
%                               must be a square number in this case

layer_type = 'single'; % pressure and velocity are computed for surface grids
%layer_type = 'double'; % only pressure is computed

R     = .07;  % radius of sphere / half edge length of cube
delta = .001; % distance between layers (if double layer)

% the incidence angles in the SDM data will be quantized to this number of
% equally distributed directions (actual number may be slightly different)
no_of_incidence_angles = 200; 

fs  = 48000;
c   = 343;
rho = 1.2; % kg/m^3

% -------------------------------------------------------------------------
% ------------------ do not change anything below this line ---------------
% -------------------------------------------------------------------------

% ------------------ avoid confusion and syntax errors -------------------

if strcmp(grid_shape, 'cubical_volume')
    layer_type = '';
    L_sphere = NaN;
end

if strcmp(grid_shape, 'spherical_surface')
    L_cube = NaN;
end

if strcmp(grid_shape, 'cubical_surface')
    L_sphere = NaN;
end

if strcmp(layer_type, 'single')
    delta = NaN;
end

% ------------------------------ load data --------------------------------

fprintf('\n');

sdm_data = load(sdm_file);

% ------------------------- do some preparations --------------------------

% check if we need to compute the velocity
if (strcmp(layer_type, 'single') && ~strcmp(grid_shape, 'cubical_volume'))
    compute_velocity = 1;
else
    compute_velocity = 0;
end

% check if sampling rates match 
sdm_data = resample_sdm_data(sdm_data, fs);

taps    = length(sdm_data.p); % length of impulse responses
taps_pw = 1024; % duration of a plane wave

% limit the number of possible incidence angles
sdm_data = quantize_sdm_angles(sdm_data, no_of_incidence_angles);

fprintf('Be aware that the computations are going to take a while.\n\n');

% % ------------- quantize sdm_data to align with the HRTF angles -----------
% % download the HRTFs if they don't exist
% hrir_path = 'hrtfs/HRIR_L2702.sofa';
% download_hrtfs(hrir_path);
% 
% SOFAstart;
% hrirs_sofa = SOFAload(hrir_path);
% head_azimuth_deg = 0; % counterclock-wise
% % get HRIR incides for all incidence angles
% indices_hrirs = SOFAfind(hrirs_sofa, sdm_data.azi_rad/pi*180 - head_azimuth_deg, sdm_data.ele_rad/pi*180);
% 
% sdm_data.azi_rad = hrirs_sofa.SourcePosition(indices_hrirs, 1)/180*pi;
% sdm_data.ele_rad = hrirs_sofa.SourcePosition(indices_hrirs, 2)/180*pi;

% ---- get the spatial grid on which to compute the pressure/velocity ----
[output_1, output_2, normal_vector] = get_sampling_grid(grid_shape, layer_type, R, L_cube, L_sphere, delta);

% sort the output data
if ~strcmp(grid_shape, 'cubical_volume') && strcmp(layer_type, 'double')
    sampling_points_inner = output_1;
    sampling_points_outer = output_2;

    % combine inner and outer layer in one variable for simpler implementation
    sampling_points = [sampling_points_inner, sampling_points_outer];
else
    sampling_points = output_1;
    % output_2 is not used
end

% -------------------------------------------------------------------------

% convert incidence direction to propagation direction
azi_prop_rad =  sdm_data.azi_rad.' + pi;
ele_prop_rad = -sdm_data.ele_rad.';

% plane wave incidence directions in Cartesian coordinates
pw_prop_cart = [cos(azi_prop_rad) .* sin(pi/2 - ele_prop_rad);
                sin(azi_prop_rad) .* sin(pi/2 - ele_prop_rad);
                cos(pi/2 - ele_prop_rad)];

all_bins = (0:taps_pw/2).';

% -------------------- loop over all plane waves --------------------------

display_progress('Computing the sound field:');

pressure = zeros(taps+taps_pw, size(sampling_points, 2));

if compute_velocity
    velocity = zeros(taps+taps_pw, size(sampling_points, 2));
end

for index = 1 : taps
    
    display_progress(index/taps);

    % allocate memory
    pressure_pw = zeros(taps_pw/2+1, size(sampling_points, 2));

    if compute_velocity
        velocity_pw = zeros(taps_pw/2+1, size(sampling_points, 2));
    end

    % Dot product to compute distance between sampling points and the wave
    % front the moment the wave front passes the origin
    % https://mathinsight.org/distance_point_plane
    % https://se.mathworks.com/matlabcentral/answers/371665-distance-from-point-to-plane-plane-was-created-from-3d-point-data
    distances = sampling_points.' * pw_prop_cart(:, index);
            
    % convert to propagation time in samples
    delay_samples = distances/c * fs;

    % loop over all sampling points
    for l = 1 : size(sampling_points, 2)

        % delay_samples is the offset to the middle of the buffer
        pressure_pw(:, l) = sdm_data.p(index) .* exp(-1i .* all_bins/taps_pw .* 2*pi .* (delay_samples(l)+taps_pw/2));

        if compute_velocity
            velocity_pw(:, l) = 1/(rho*c) .* pressure_pw(:, l) .* dot(pw_prop_cart(:, index), normal_vector(:, l));
        end

    end

    % this comprises the current plane wave
    pressure_pw = ifft([pressure_pw; conj(flipud(pressure_pw(2:end-1, :)))], [], 1, 'symmetric');

    % add the plane wave to the rest of the sound field
    pressure(index:index+taps_pw-1, :) = pressure(index:index+taps_pw-1, :) + pressure_pw;

    if compute_velocity
        % this comprises the current plane wave
        velocity_pw = ifft([velocity_pw; conj(flipud(velocity_pw(2:end-1, :)))], [], 1, 'symmetric');

        % add the plane wave to the rest of the sound field
        velocity(index:index+taps_pw-1, :) = velocity(index:index+taps_pw-1, :) + velocity_pw;
    end

end

% ----------------------- window the irs, just in case --------------------
win      = hann(round_up_to_even(500*fs/48000)); 
fade_in  = win(1:end/2);
fade_out = win(end/2+1:end);

pressure(1:length(fade_in),          :, :) = pressure(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(pressure, 2), size(pressure, 3));
pressure(end-length(fade_out)+1:end, :, :) = pressure(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(pressure, 2), size(pressure, 3));

if compute_velocity
    velocity(1:length(fade_in),          :, :) = velocity(1:length(fade_in),          :, :) .* repmat(fade_in,  1, size(velocity, 2), size(velocity, 3));
    velocity(end-length(fade_out)+1:end, :, :) = velocity(end-length(fade_out)+1:end, :, :) .* repmat(fade_out, 1, size(velocity, 2), size(velocity, 3));
end

fprintf('\n\n');

% ----------------------------- plot data ---------------------------------

plot_grid(layer_type, sampling_points);

% plot pressure 
figure;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [1200 100 800 500]);

if compute_velocity
    subplot(1, 2, 1);
end

plot(pressure);
grid on;
xlabel('t (samples)');
title('Pressure at grid points');

if compute_velocity
    subplot(1, 2, 2);
    plot(velocity);
    grid on;
    xlabel('t (samples)');
    title('Velocity at grid points');
end


% ---------------------------- store data ---------------------------------
fprintf('Converting to single precision.\n\n'); % to save space
        
pressure = single(pressure);

if compute_velocity
    velocity = single(velocity);
end

if strcmp(grid_shape, 'cubical_volume')

    file_name = sprintf('sound_field_p_%s_%s_L%d.mat', grid_shape, sdm_data.room, size(sampling_points, 2));

    save(file_name, 'pressure', 'sampling_points', 'grid_shape', 'fs', 'layer_type', '-v7.3');

elseif strcmp(layer_type, 'double')

    pressure_inner = pressure(:, 1:end/2);
    pressure_outer = pressure(:, end/2+1:end); 

    file_name = sprintf('sound_field_pp_%s_%s_L%d.mat', grid_shape, sdm_data.room, size(sampling_points, 2)/2);

    save(file_name, 'pressure_inner', 'pressure_outer', 'sampling_points_inner', 'sampling_points_outer', 'grid_shape', 'fs', 'layer_type', '-v7.3');
    
elseif strcmp(layer_type, 'single')

    file_name = sprintf('sound_field_pv_%s_%s_L%d.mat', grid_shape, sdm_data.room, size(sampling_points, 2));

    save(file_name, 'pressure', 'velocity', 'sampling_points', 'normal_vector', 'grid_shape', 'fs', 'layer_type', '-v7.3');

else
    error('Something is wrong here.');
end

fprintf('Stored the data in the file ''%s''.\n\n', file_name);


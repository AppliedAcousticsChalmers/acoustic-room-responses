% This script loads an SDM representation of an acoustic room response and
% converts it into an ambisonic representation by imposing the digital
% samples of the pressure room response onto plane waves.

clear;

addpath('dependencies/');

% SDM data file
%sdm_file = 'room_data/sdm_big_hall.mat';
sdm_file = 'room_data/sdm_listening_lab.mat';

% the incidence angles in the SDM data will be quantized to this number of
% equally distributed directions (actual number may be slightly different)
no_of_incidence_angles = 200;

N = 32; % spherical harmonic order of the simulation 

fs = 48000;
c  = 343;

% -------------------------------------------------------------------------
% ------------------ do not change anything below this line ---------------
% -------------------------------------------------------------------------

% ------------------------------ load data --------------------------------

fprintf('\n');

sdm_data = load(sdm_file);

% ------------------------- do some preparations --------------------------

fprintf('Be aware that the computations are going to take a while.\n\n');

% check if sampling rates match 
sdm_data = resample_sdm_data(sdm_data, fs);

taps = length(sdm_data.p);

% limit the number of possible incidence angles
sdm_data = quantize_sdm_angles(sdm_data, no_of_incidence_angles);

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

% ------------------- Compute the SH coefficients -------------------------

display_progress('Computing the ambisonic signal:');

s_nm = zeros(length(sdm_data.p), (N+1)^2);

for index = 1 : taps
    
    display_progress(index/taps);
    
    for n = 0 : N
        for m = -n : n
            s_nm(index, n^2+n+m+1) = s_nm(index, n^2+n+m+1) + sdm_data.p(index) .* sphharm(n, m, pi/2-sdm_data.ele_rad(index), sdm_data.azi_rad(index), 'real'); 
        end
    end
    
end

fprintf('\n\n');

room = sdm_data.room;

% create a useful file name
sdm_path_components = strsplit(sdm_file, '/');
sdm_id              = sdm_path_components{end};
file_out_name       = ['ambisonics_' sdm_id(1:end-4) '.mat'];

fprintf('Storing the data in the file ''%s'' ... ', file_out_name);

save(['ambisonics_' sdm_id(1:end-4) '.mat'], 'fs', 's_nm', 'room', 'c');

fprintf('done.\n\n');




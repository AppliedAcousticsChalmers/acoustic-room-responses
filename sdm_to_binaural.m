% This script loads an SDM representation of an acoustic room response and
% converts it into a binaural signal by imposing the digital samples of the 
% pressure room response onto HRTFs.

clear;

addpath('dependencies/');

% SDM data file
%sdm_file = 'room_data/sdm_big_hall.mat';
sdm_file = 'room_data/sdm_listening_lab.mat';

head_azimuth_deg = 0; % counterclock-wise in degrees

% the incidence angles in the SDM data will be quantized to this number of
% equally distributed directions (actual number may be slightly different)
no_of_incidence_angles = 200; 

% -------------------------------------------------------------------------
% ------------------ do not change anything below this line ---------------
% -------------------------------------------------------------------------

% ------------------------------ load data --------------------------------

fprintf('\n');

sdm_data = load(sdm_file);

% download the HRTFs if they don't exist
hrir_path = 'hrtfs/HRIR_L2702.sofa';
download_hrtfs(hrir_path);

SOFAstart;
hrirs_sofa = SOFAload(hrir_path);

% ------------------------- do some preparations --------------------------
fs = double(hrirs_sofa.Data.SamplingRate);

% check if sampling rates match 
sdm_data = resample_sdm_data(sdm_data, fs);

taps = length(sdm_data.p);

% limit the number of possible incidence angles
sdm_data = quantize_sdm_angles(sdm_data, no_of_incidence_angles);

taps_hrirs = size(hrirs_sofa.Data.IR, 3);
brirs = zeros(taps + taps_hrirs - 1, 2);

% get HRIR incides for all incidence angles (this determines the actual
% qualitization of the incidence diretions!)
indices_hrirs = SOFAfind(hrirs_sofa, sdm_data.azi_rad/pi*180 - head_azimuth_deg, sdm_data.ele_rad/pi*180);

display_progress('Computing the binaural signal:');

% ------------- loop over all digital samples of the SDM data -------------
for counter = 1 : taps
    
    display_progress(counter/taps);
    
    brirs(counter : counter+taps_hrirs-1, :) = brirs(counter : counter+taps_hrirs-1, :) + sdm_data.p(counter) .* squeeze(hrirs_sofa.Data.IR(indices_hrirs(counter), :, :)).';
    
end

fprintf('\n\n');

% ----------------------------- audio example -----------------------------
[sig, fs_audio] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_audio);

out_gt = [fftfilt(brirs(:, 1), sig), fftfilt(brirs(:, 2), sig)];
out_gt = out_gt/max(abs(out_gt(:))) .* .99;

% create a useful file name
sdm_path_components = strsplit(sdm_file, '/');
sdm_id              = sdm_path_components{end};
file_out_name       = ['binaural_' sdm_id(1:end-4) '.wav'];

fprintf('Storing the binaural signal in the file ''%s'' ... ', file_out_name);

audiowrite(file_out_name, out_gt, fs);

fprintf('done.\n\n');

% This script loads an ambisonic soud field representation (for example, 
% created with sdm_data_to_ambisonics.m and computes the signals that the 
% microphones of a rigid spherical array would capture due to the sound
% field.

clear;

addpath('dependencies/');

% execute script sdm_data_to_ambisonics.m to obtain these data
ambisonics_file_name = 'ambisonics_sdm_listening_lab.mat';
%ambisonics_file_name = 'ambisonics_sdm_big_hall.mat';

% define microphone array (Eigenmike)
[azi_rad, col_rad, R, grid_weights_array, N_sma] = get_eigenmike_mic_pos();

% -------------------------------------------------------------------------
% ------------------ do not change anything below this line ---------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
fprintf('\n');

load(ambisonics_file_name);

% get SH order of the ambisonic room response
N = round(sqrt(size(s_nm, 2))) - 1;

% make even length
if rem(size(s_nm, 1), 2)
    s_nm = [s_nm; zeros(1, size(s_nm, 2))];
end

% frequency domain
S_nm = fft(s_nm);
S_nm = S_nm(1:end/2+1, :);

f = linspace(0, fs/2, size(S_nm, 1)).';
f(1) = 1e-30; % avoid NaN
k = 2*pi*f/c;
kR = k*R;

hankel_type = 2;

% microphone transfer functions
S = zeros(size(S_nm, 1), size(azi_rad, 2));

display_progress('Computing the microphone signals:');

for n = 0 : N

    % compute radial terms
    hankel_prime = 1/(2*n+1) * (n * sphbesselh(n-1, hankel_type, kR) - (n+1) * sphbesselh(n+1, hankel_type, kR));
    b_n          = -4*pi*1i^n .* 1i./kR.^2 .* 1./hankel_prime;

    for m = -n : n

        display_progress((n+1)^2/N^2);

        S = S + S_nm(:, n^2+n+m+1) .* b_n .* sphharm(n, m, col_rad, azi_rad);
   
    end
end

fprintf('\n\n');

% fix DC
S(1, :) = real(S(2, :));

s = ifft([S; conj(flipud(S(2:end-1, :)))], 'symmetric');

% ---------------------- create example recording -------------------------
[sig, fs_audio] = audioread('resources/drum_loop_48k.wav');

assert(fs == fs_audio);

s_drums = fftfilt(s, sig);

% normalize
s_drums = s_drums/max(abs(s_drums(:))) .* .99;

% -------------------------- store everything -----------------------------

file_out_name = ['mic_array_irs_' ambisonics_file_name(12:end-4) '.mat'];

fprintf('Storing the data in the files ''%s'' ... ', file_out_name);

save(file_out_name, 'fs', 's', 'room');

file_out_name = ['mic_array_recording_' ambisonics_file_name(12:end-4) '.mat'];

fprintf('and ''%s'' ... ', file_out_name);

% convert variable names to what render_sma_to_ambisonics.m requires
alpha_sma     = azi_rad; 
beta_sma      = col_rad;
N             = N_sma;
array_signals = s_drums;

save(file_out_name, 'array_signals', 'fs', 'R', 'alpha_sma', 'beta_sma', 'N', 'grid_weights_array', '-v7.3');

fprintf('done.\n\n');



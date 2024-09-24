function [azi_rad, col_rad, R, grid_weights, N] = get_eigenmike_mic_pos()
% azi_rad, col_rad: Azimuth and colatitude (zenith angle) in radians of the
%                   microphone positions
% R:                Radius of the microphone array
% grid_weights:     Quadrature weights of the microphones (normalized to
%                   \sum grid_weights = 1)
% N:                Maximum SH order supported by the grid

R = 0.042; % m

% col, azi
mic_angles = [ 69 0;
               90 32
              111 0;
               90 328;
               32 0;
               55 45;
               90 69;
               125 45;
               148 0;
               125 315;
               90 291;
               55 315;
               21 91;
               58 90;
               121 90;
               159 89;
               69 180;
               90 212;
               111 180;
               90 148;
               32 180;
               55 225;
               90 249;
               125 225;
               148 180;
               125 135;
               90 111;
               55 135;
               21 269;
               58 270;
               122 270;
               159 271 ];

col_rad = mic_angles(:, 1).'/180*pi;
azi_rad = mic_angles(:, 2).'/180*pi;

N = 4;
grid_weights = ones(size(azi_rad))/32;

% [x, y, z] = sph2cart(azi_rad, pi/2-col_rad, R);
%
% figure;
% plot3(x, y, z, '.');
% grid on;
% axis equal;

end




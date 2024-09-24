function [sdm_data] = resample_sdm_data(sdm_data, fs)

if (sdm_data.fs ~= fs)
    
    fprintf('Resampling SDM data from %d Hz to %d Hz ... ', sdm_data.fs, fs);
    
    sdm_data.p       = resample(sdm_data.p, fs, sdm_data.fs);
    sdm_data.azi_rad = resample(unwrap(sdm_data.azi_rad), fs, sdm_data.fs);
    sdm_data.ele_rad = resample(sdm_data.ele_rad, fs, sdm_data.fs);
    sdm_data.fs      = fs;
    
    fprintf('done.\n\n');

end

end


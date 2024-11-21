function [loc_beam_time] = get_beamtime(tdoa_measured_select,datetime_select,beam_threshold)
% function get_beamtime.m finds the time when the tdoa tracks cross the
% beam (if they cross it).

% Find valid (non-NaN) measurements
    valid_indices = find(~isnan(tdoa_measured_select));
    valid_tdoa = tdoa_measured_select(valid_indices);
    valid_times = datetime_select(valid_indices);

   
    if length(valid_tdoa) >= 2
        % First check if any values are within beam threshold 
        [min_val, min_idx] = min(abs(valid_tdoa));

        if min_val < beam_threshold
            % Found value within threshold
            loc_beam_time = datetime(valid_times(min_idx), 'ConvertFrom','datenum');
        else
            % Check for zero crossings
            zero_crossings = find(valid_tdoa(1:end-1) .* valid_tdoa(2:end) <= 0);

            if ~isempty(zero_crossings)
                % Linear interpolation to find more precise zero crossing time
                idx = zero_crossings(1);  % Use first zero crossing
                t1 = valid_times(idx);
                t2 = valid_times(idx+1);
                v1 = valid_tdoa(idx);
                v2 = valid_tdoa(idx+1);

                % Interpolate to find exact zero crossing time
                zero_fraction = abs(v1) / (abs(v1) + abs(v2));
                zero_time = t1 + (t2 - t1) * zero_fraction;

                loc_beam_time = datetime(zero_time, 'ConvertFrom','datenum');
            else
                % No zero crossing found
                loc_beam_time = datetime(7e+05, 'ConvertFrom','datenum');
                %I cannot assign Nan value to this- needs to be datetime array, so
                %the times where the tracks dont cross the beam will be
                %'14-Jul-1916 00:00:00'.
            end
        end
    else
        % Not enough valid measurements
        loc_beam_time = datetime(7e+05, 'ConvertFrom','datenum');
    end
end
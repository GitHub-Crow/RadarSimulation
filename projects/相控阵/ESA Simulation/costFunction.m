function cost = costFunction(complex_wgts, ele_wgts_voltages,...
                             lower_bound, upper_bound,...
                             auto_exempt_main_beam, enforce_symmetry)

 Nele = numel(complex_wgts);
 if enforce_symmetry
     complex_wgts(Nele/2 + 1 : Nele) = flipud(complex_wgts(1 : Nele/2));
 end
 
 Npattern_points = numel(ele_wgts_voltages);
 % calculate array factor by fft directly
 pattern = abs(fftshift(fft(complex_wgts, Npattern_points))) .* ele_wgts_voltages;
 
 pattern_norm = pattern/max(pattern);
 lower_excess = max(lower_bound - pattern_norm, 0);
 upper_excess = max(pattern_norm - upper_bound, 0);
 
 if auto_exempt_main_beam
     [~,index_of_max_value] = max(pattern_norm);
     spikes = diff(sign(diff(pattern_norm)));
     % find loction of main beam
     null_indices = find(spikes == 2) + 1; % lower points
     right_null_indices = null_indices(null_indices > index_of_max_value(1));
     left_null_indices = null_indices(null_indices < index_of_max_value(1));
     
     if isempty(right_null_indices)
         index_of_right_null = Npattern_points;
     else
         index_of_right_null = right_null_indices(1);
     end
     if isempty(left_null_indices)
         index_of_left_null = 1;
     else
         index_of_left_null = left_null_indices(end);
     end
     % exempt main beam from cost function
     lower_excess(index_of_left_null : index_of_right_null) = 0;
     upper_excess(index_of_left_null : index_of_right_null) = 0;
 end
 
 cost = 10*log10(sum(lower_excess) + sum(upper_excess) + eps);
end
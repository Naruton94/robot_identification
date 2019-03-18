d = [];
for i = 1:6
    d1 = (1 - norm(tau_pos(:, i) - tau_filt(:, i)) / norm(tau_pos(:, i)));
    d = [d, d1];
end
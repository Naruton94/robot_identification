N = 3;
% fake_err = binornd (2, 0.6 ,size(tau_pos));
% x = sign(rand(size(fake_err))-0.5);

% fake_err = fake_err .* x;
% tau_fake = tau_pos + fake_err;

tau_origin = [];

for i = 1: 6
    fake_err1 = unifrnd (-1,1, size(tau_pos));
    fake_err2 = randn(size(tau_pos));
    tau_fake = tau_pos + tau_pos .* fake_err1 * 0.2 + tau_pos .* fake_err2 * 0.2;
    dd = smooth(tau_fake(:, i), 89);
    dd = medfilt1(dd, 55);
    plot(dd)
    hold on
    plot(tau_filt(:, i))
    tau_origin = [tau_origin, dd];
end
% save tau_origin tau_origin

xlswrite("tau_pos.xlsx", tau_origin);


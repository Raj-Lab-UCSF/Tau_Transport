function f = unscale_outputs_f(f_scaled)
    f = (f_scaled .* 3.216671e-09) + -2.330046e-10; % Inverse Z-score from dataset NN is trained on for frac = 0.7
end
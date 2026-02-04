function w = unscale_outputs_w(w_scaled)
    w = (w_scaled .* 2.583219e-06) + 2.096271e-06; % Inverse Z-score from dataset NN is trained on for frac = 0.7
end
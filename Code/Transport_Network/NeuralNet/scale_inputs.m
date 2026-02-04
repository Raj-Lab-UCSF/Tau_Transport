function X_scaled = scale_inputs(X)
    X_scaled = zeros(size(X));
    X_scaled(:,1) = (X(:,1) - 7.254e-03) ./ 1.1345e-02;
    X_scaled(:,2) = (X(:,2) - 3.8941e-02) ./ 2.4867e-02;
    X_scaled(:,3) = (X(:,3) - 55.049541) ./ 25.953783;
    X_scaled(:,4) = (X(:,4) - 54.975716) ./ 26.067346;
    X_scaled(:,5) = (X(:,5) - 1.003e-03) ./ 7.09e-04;
    X_scaled(:,6) = (X(:,6) - 1.001e-03) ./ 7.07e-04;
end
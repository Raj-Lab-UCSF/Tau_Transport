function [net] = get_nn_model_pytorch(filepath, input_size)
    net = importNetworkFromPyTorch(filepath, 'PyTorchInputSizes', input_size);
    X_init = dlarray(rand(input_size), 'UUU');
    net = initialize(net, X_init);
end
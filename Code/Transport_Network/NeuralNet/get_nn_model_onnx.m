function [net] = get_nn_model_onnx (filepath, input_size)
    net = importNetworkFromONNX(filepath);
    X_init = dlarray(rand(input_size), 'UUU');
    net = initialize(net, X_init);
end
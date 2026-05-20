import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import mean_squared_error, r2_score

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, TensorDataset

from collections import OrderedDict

from scipy.stats import pearsonr

import pickle

import sys

from my_nn import NeuralNetwork
from nn_hyp_utils import model_fit

data_path = sys.argv[1]

target = sys.argv[2]
target = [target]

data_save_file = sys.argv[3]

random_state = int(sys.argv[4])
hidden_size_in = int(sys.argv[5])
hidden_layers_in = int(sys.argv[6])
epochs_in = int(sys.argv[7])
batch_size = int(sys.argv[8])

val_data_path = sys.argv[9]
model_save_file = sys.argv[10]

fold_num = int(sys.argv[11])

learn_rate_in = float(sys.argv[12])

# Select device to run job on

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('Using GPU')
else:
    device = torch.device('cpu')
    print('GPU not available, using CPU instead')

# Fit Scaling Function using training dataset

fit_data_filepath = data_path

scale_data = pd.read_csv(fit_data_filepath)

features = ['gamma1', 'lambda1', 'delta', 'epsilon', 'NRow', 'NCol']

X_scale = scale_data[features].values
y_scale = scale_data[target].values

PredictorScaler = StandardScaler()
TargetScaler = StandardScaler()

PredScaleFit = PredictorScaler.fit(X_scale)
TargetScaleFit = TargetScaler.fit(y_scale)

# Define model architecture and instantiate neural network

input_size = 6
output_size = 1

hidden_size = hidden_size_in
hidden_layers = hidden_layers_in

epoch_chunk_size = 50

epochs_n = epochs_in

size = (int(epochs_n/epoch_chunk_size), 6)

#output_data = np.zeros(size)

'''
model = NeuralNetwork(input_size, hidden_size, hidden_layers, output_size)
model = model.to(device)

criterion = nn.MSELoss()
optimizer = optim.SGD(model.parameters())
'''

# Prepare Training Dataset

data = pd.read_csv(data_path)

features = ['gamma1', 'lambda1', 'delta', 'epsilon', 'NRow', 'NCol']

X = data[features].values
y = data[target].values

X_train_scaled = PredScaleFit.transform(X)
y_train_scaled = TargetScaleFit.transform(y)

k_folds = 10
kf = KFold(n_splits=k_folds, shuffle=True, random_state=random_state)

size = (int(epochs_n/epoch_chunk_size), 2, 6)

output_data = np.zeros(size)

# Prepare Validation Dataset

data_val = pd.read_csv(val_data_path)

X_val = data_val[features].values
y_val = data_val[target].values

X_val_scaled = PredScaleFit.transform(X_val)
X_val_pt = torch.tensor(X_val_scaled, dtype=torch.float32).to(device)

# Train Model

for fold, (train_index, test_index) in enumerate(kf.split(X_train_scaled, y_train_scaled)):

    # Skip incorrect folds - Only run one fold per job

    if (fold + 1) != fold_num:
        continue

    # Create train-test data splits

    X_train, X_test = X_train_scaled[train_index], X_train_scaled[test_index]
    y_train, y_test = y_train_scaled[train_index], y_train_scaled[test_index]

    X_train_pt = torch.tensor(X_train, dtype=torch.float32).to(device)
    X_test_pt = torch.tensor(X_test, dtype=torch.float32).to(device)
    y_train_pt = torch.tensor(y_train, dtype=torch.float32).to(device)

    # Initialize the model, loss criteria, and optimization technique
            
    model = NeuralNetwork(input_size, hidden_size, hidden_layers, output_size)
    model = model.to(device)

    criterion = nn.MSELoss()
    optimizer = optim.SGD(model.parameters(), lr=learn_rate_in)

    data_index = 0

    min_metric = float('inf')
    best_model = model

    for epoch_i in range(epochs_n):

        dataset = TensorDataset(X_train_pt, y_train_pt)

        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        for batch_i, (inputs, targets) in enumerate(dataloader):

            # Forward pass
            outputs = model(inputs)
            loss = criterion(outputs, targets)

            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Print progress
        print(f'Epoch [{epoch_i+1}/{epochs_n}], Loss: {loss.item():.4f}')

        if (epoch_i + 1) % epoch_chunk_size == 0:

            # predict fold test data and get error metrics

            preds_test_scaled = model(X_test_pt)
            preds_test = TargetScaleFit.inverse_transform(preds_test_scaled.detach().numpy())

            y_test_unscaled = TargetScaleFit.inverse_transform(y_test)
            metrics_test = model_fit(y_test_unscaled.flatten(), preds_test.flatten(), verbose=False)
        
            # predict validation data and get error metrics

            preds_val_scaled = model(X_val_pt)
            preds_val = TargetScaleFit.inverse_transform(preds_val_scaled.detach().numpy())

            metrics_val = model_fit(y_val.flatten(), preds_val.flatten(), verbose=False)
        
            output_data[data_index, 0, :] = metrics_test
            output_data[data_index, 1, :] = metrics_val
            data_index += 1

            with open(data_save_file, 'wb') as f:
                pickle.dump(output_data, f)

            if metrics_test[4] < min_metric:
                best_model = copy.deepcopy(model)
                min_metric = metrics_test[4]

                dummy_input = torch.randn((1,6), dtype=torch.float32)
                traced_model = torch.jit.trace(best_model.forward, dummy_input)

                model_filename = model_save_file + ".pt"

                traced_model.save(model_filename)          
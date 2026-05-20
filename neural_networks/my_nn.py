import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

import torch
import torch.nn as nn
import torch.optim as optim

from collections import OrderedDict

class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, hidden_layers, output_size):
        super().__init__()

        layers = []
        layers.append(('layer_first', nn.Linear(input_size, hidden_size)))
        layers.append(('relu_1', nn.ReLU()))
        for i in range(hidden_layers-1):
            layers.append((f'layer_hidden_{i+1}', nn.Linear(hidden_size, hidden_size)))
            layers.append((f'relu_{i+2}', nn.ReLU()))
        layers.append(('layer_last', nn.Linear(hidden_size, output_size)))

        model = nn.Sequential(OrderedDict(layers))

        self.linear_relu_stack = model

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits
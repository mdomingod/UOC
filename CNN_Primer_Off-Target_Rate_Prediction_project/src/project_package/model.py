# src/project_package/model.py
import torch
import torch.nn as nn
import torch.nn.functional as F
from .config import Config

class PordleCNN(nn.Module):
    def __init__(self, cfg: Config):
        super().__init__()
        f = cfg.n_filters
        self.conv1 = nn.Conv2d(cfg.img_channels, f, kernel_size=5,      padding='same')
        self.conv2 = nn.Conv2d(f,                f, kernel_size=5,      padding='same')
        self.pool1 = nn.MaxPool2d(kernel_size=2, stride=2)
        self.conv3 = nn.Conv2d(f, f, kernel_size=(3, 5), padding='same')
        self.conv4 = nn.Conv2d(f, f, kernel_size=(3, 5), padding='same')
        self.pool2 = nn.MaxPool2d(kernel_size=2, stride=2, ceil_mode=True)
        self.dropout = nn.Dropout(cfg.dropout_rate)
        self.fc1 = nn.Linear(2 * 17 * f, cfg.fc_units)
        self.fc2 = nn.Linear(cfg.fc_units, cfg.fc_units)
        self.fc_out = nn.Linear(cfg.fc_units, 1)
        self.allow_negative = cfg.allow_negative_output

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = self.pool1(x)
        x = F.pad(x, (0, 0, 0, 1))
        x = self.dropout(x)
        x = F.relu(self.conv3(x))
        x = F.relu(self.conv4(x))
        x = self.pool2(x)
        x = self.dropout(x)
        x = x.flatten(start_dim=1)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc_out(x).squeeze(-1)
        return x if self.allow_negative else F.relu(x)

    def count_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
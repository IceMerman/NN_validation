import scipy.io
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# Load data from .mat file
mat = scipy.io.loadmat('data.mat')
X = mat['X'] # Features (power demand, voltage magnitude, angle, etc.)
y = mat['y'] # Classification labels (safe or not safe)

# Convert data to PyTorch tensors
X = torch.tensor(X, dtype=torch.float32)
y = torch.tensor(y, dtype=torch.float32)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Define the neural network architecture
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(X.shape[1], 64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.sigmoid(self.fc3(x))
        return x

# Create the neural network
net = Net()

# Define the loss function and optimizer
criterion = nn.BCELoss()
optimizer = optim.SGD(net.parameters(), lr=0.1)

# Train the neural network
for epoch in range(1000):
    optimizer.zero_grad()
    outputs = net(X_train)
    loss = criterion(outputs.view(-1), y_train)
    loss.backward()
    optimizer.step()

# Predict on the testing set
y_pred = (net(X_test).view(-1) > 0.5).float()

# Evaluate accuracy
acc = accuracy_score(y_test, y_pred)
print("Accuracy: {:.2f}%".format(acc * 100))

import scipy.io
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score

# Load data from .mat file
mat = scipy.io.loadmat('data.mat')
X = mat['X'] # Features (power demand, voltage magnitude, angle, etc.)
y = mat['y'] # Classification labels (safe or not safe)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Train an SVM classifier
clf = SVC()
clf.fit(X_train, y_train)

# Predict on the testing set
y_pred = clf.predict(X_test)

# Evaluate accuracy
acc = accuracy_score(y_test, y_pred)
print("Accuracy: {:.2f}%".format(acc * 100))

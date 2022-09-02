import math
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tensorflow.keras import Model
from tensorflow.keras import Sequential
from tensorflow.keras.optimizers import Adam
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.layers import Dense, Dropout
from sklearn.model_selection import train_test_split
from tensorflow.keras.losses import MeanSquaredLogarithmicError


TRAIN_DATA_PATH = './matlab/data models/4 vars/ide_set_all.csv'
TEST_DATA_PATH = './matlab/data models/4 vars/val_set_all.csv'

TARGET_NAME = 'integral'
print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")
# x_train = features, y_train = target
train_data = pd.read_csv(TRAIN_DATA_PATH)
train_data = train_data.drop('priority', axis=1)
train_data = train_data.drop('max_delta', axis=1)
train_data = train_data.drop('pi', axis=1)

test_data = pd.read_csv(TEST_DATA_PATH)
test_data = test_data.drop('priority', axis=1)
test_data = test_data.drop('pi', axis=1)
test_data = test_data.drop('max_delta', axis=1)

x_train, y_train = train_data.drop(TARGET_NAME, axis=1), train_data[TARGET_NAME]
x_test, y_test = test_data.drop(TARGET_NAME, axis=1), test_data[TARGET_NAME]

union = pd.concat([x_train, x_test], ignore_index=True)
print(union)
def scale_datasets(x_train, x_test):
    # Standard Scale test and train data
    # Z - Score normalization
    standard_scaler = StandardScaler()
    x_train_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_train),
        columns=x_train.columns
    )
    #print("mean: " + str(standard_scaler.mean_))
    #print("var: " + str(standard_scaler.var_))
    #print(x_train_scaled)
    #data_orig = standard_scaler.inverse_transform(x_train_scaled)
    #print(data_orig)


    x_test_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_test),
        columns=x_test.columns
    )
    x_test_orig = pd.DataFrame(
        standard_scaler.inverse_transform(x_test),
        columns=x_test.columns
    )
    #print("mean: " + str(standard_scaler.mean_))
    #print("var: " + str(standard_scaler.var_))


    return x_train_scaled, x_test_scaled, x_test_orig


x_train_scaled, x_test_scaled, x_test_orig = scale_datasets(x_train, x_test)

hidden_units1 = 160
hidden_units2 = 480
hidden_units3 = 256
learning_rate = 0.01


# Creating model using the Sequential in tensorflow
# Dense = fully connected layer, with num of outputs; Dropout helps avoid overfitting
def build_model_using_sequential():
    mdl = Sequential([
        Dense(hidden_units1, kernel_initializer='normal', activation='relu'),
        Dropout(0.2),
        Dense(hidden_units2, kernel_initializer='normal', activation='relu'),
        Dropout(0.2),
        Dense(hidden_units3, kernel_initializer='normal', activation='relu'),
        Dense(1, kernel_initializer='normal', activation='linear')
    ])
    return mdl


# build the model
model = build_model_using_sequential()

# loss function
msle = MeanSquaredLogarithmicError()
# model settings for learning
# Adam is a stochastic gradient descent method based on adaptive estimation of first-order and second-order moments
model.compile(
    loss=msle,
    optimizer=Adam(learning_rate=learning_rate),
    metrics=[msle]
)
# train the model
history = model.fit(
    x_train_scaled.values,
    y_train.values,
    epochs=500,
    batch_size=936,
    validation_split=0.2
)


def plot_history(history, key):
    plt.plot(history.history[key])
    plt.plot(history.history['val_' + key])
    plt.xlabel("Epochs")
    plt.ylabel(key)
    plt.legend([key, 'val_' + key])
    plt.grid(True)
    plt.show()


# Plot the history
plot_history(history, 'mean_squared_logarithmic_error')


print(x_test)


output_predicted = model.predict(x_test_scaled.values)

x_test['prediction'] = model.predict(x_test_scaled)

#output_predicted = pd.DataFrame(output_predicted, columns = ['pred integral'])
#print(x_test)
#data_orig = pd.concat([data_orig, output_predicted], axis=1)
#print(data_orig)

x_test.to_csv("pino.csv", index=False, float_format='%g', encoding="utf-8")
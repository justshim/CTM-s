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
from tensorflow.keras.losses import *


TRAIN_DATA_PATH = '../../data/opti_data_all_priority.csv'
TEST_DATA_PATH = '../../data/opti_data_all_validazione.csv'

TARGET_NAME_1 = 'integral'
TARGET_NAME_2 = 'pi'
print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")
# x_train = features, y_train = target
train_data = pd.read_csv(TRAIN_DATA_PATH)
train_data = train_data.drop(['max_delta'], axis=1)

test_data = pd.read_csv(TEST_DATA_PATH)
test_data = test_data.drop(['max_delta'], axis=1)

y_train = train_data[[TARGET_NAME_1, TARGET_NAME_2]]
y_test = test_data[[TARGET_NAME_1, TARGET_NAME_2]]
#print(y_test)

y_train['pi']=y_train['pi']*100
y_test['pi']=y_test['pi']*100

x_train = train_data.drop([TARGET_NAME_1, TARGET_NAME_2], axis=1)
x_test = test_data.drop([TARGET_NAME_1, TARGET_NAME_2], axis=1)

print(x_train)
#print(y_test)
def scale_datasets(x_train, x_test):
    # Standard Scale test and train data
    # Z - Score normalization
    standard_scaler = StandardScaler()
    x_train_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_train),
        columns=x_train.columns
    )
    x_test_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_test),
        columns=x_test.columns
    )
    return x_train_scaled, x_test_scaled


x_train_scaled, x_test_scaled = scale_datasets(x_train, x_test)

hidden_units1 = 15
hidden_units2 = 10
hidden_units3 = 150
learning_rate = 0.1


# Creating model using the Sequential in tensorflow
# Dense = fully connected layer, with num of outputs; Dropout helps avoid overfitting
def build_model_using_sequential():
    mdl = Sequential([
        Dense(15, input_shape=(15,), activation='relu'),
        Dense(hidden_units1, kernel_initializer='he_normal', activation='relu'),
        Dropout(0.2),
        Dense(hidden_units2, kernel_initializer='he_normal', activation='relu'),
        Dropout(0.2),
        #Dense(hidden_units3, kernel_initializer='he_normal', activation='relu'),
        #Dropout(0.2),
        Dense(2)
    ])
    return mdl


# build the model
model = build_model_using_sequential()

# loss function
mae = MeanAbsoluteError()
msle = MeanSquaredLogarithmicError()
mse = MeanSquaredError()
# model settings for learning
# Adam is a stochastic gradient descent method based on adaptive estimation of first-order and second-order moments
model.compile(
    loss=msle,
    #optimizer=Adam(learning_rate=learning_rate)
    optimizer=Adam(),
    metrics=['accuracy']
)
# train the model
history = model.fit(
    x_train_scaled.values,
    y_train.values,
    epochs=60,
    batch_size=256,
    #validation_data=(x_test_scaled, y_test),
    validation_split=0.2
)


#output_predicted = model.predict(x_test_scaled.values)

model.save('SavedModel/modellovdue', save_format="h5")

def plot_history(history):
    plt.figure(1)
    plt.subplot(211)
    plt.plot(history.history['acc'])
    plt.plot(history.history['val_acc'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.grid(True)

    plt.subplot(212)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.grid(True)
    plt.show()


# print("\nEvaluate")
# result = model.predict(x_test_scaled.values)
# print("Results: ")
# print(result)
# dict(zip(model.metrics_names, result))

# Plot the history
plot_history(history)

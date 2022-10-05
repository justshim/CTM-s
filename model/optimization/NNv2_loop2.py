import csv
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

INPUT_PATH = '../../data/input_all_all.csv'
OUTPUT_PATH = '../../data/output_all_all.csv'

print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")

x_train = pd.read_csv(INPUT_PATH, sep=";", header=None)
y_train = pd.read_csv(OUTPUT_PATH, sep=";", header=None)

#print(x_train)
#print(y_train)


def scale_datasets(x_t):
    # Standard Scale test and train data
    # Z - Score normalization
    standard_scaler = StandardScaler()
    x_t_s = pd.DataFrame(
        standard_scaler.fit_transform(x_t),
        columns=x_t.columns
    )
    return x_t_s


x_train_scaled = scale_datasets(x_train)

for b in range(400, 800, 16):
    if b == 0:
        continue

    hidden_units1 = 55
    hidden_units2 = 0
    hidden_units3 = 0
    # learning_rate = 0.001


    def input_flatten(df):
        output = pd.DataFrame()
        for index in range(0, len(df.index), 15):
            matrix = df.loc[index:index + 14, :].to_numpy()
            flat = matrix.flatten()
            aa = pd.DataFrame(flat).transpose()
            output = pd.concat([output, aa])
        return output

    input = input_flatten(x_train_scaled)


    # Creating model using the Sequential in tensorflow
    # Dense = fully connected layer, with num of outputs; Dropout helps avoid overfitting
    def build_model_using_sequential():
        mdl = Sequential([
            Dense(90, input_shape=(90,), activation='relu'),
            Dense(hidden_units1, kernel_initializer='he_normal', activation='relu'),
            Dropout(0.2),
            # Dense(hidden_units2, kernel_initializer='he_normal', activation='relu'),
            # Dropout(0.2),
            # Dense(hidden_units3, kernel_initializer='he_normal', activation='relu'),
            # Dropout(0.2),
            Dense(3)
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
        # optimizer=Adam(learning_rate=learning_rate)
        optimizer=Adam(),
        metrics=[msle]
    )
    # train the model
    history = model.fit(
        input,
        y_train.values,
        epochs=250,
        batch_size=b,
        # validation_data=(x_test_scaled, y_test),
        validation_split=0.2
    )

    # output_predicted = model.predict(x_test_scaled.values)

    model.save('SavedModel/modellovdue', save_format="h5")


    def plot_history(history):
        plt.figure(1)
        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'val'], loc='upper left')
        plt.grid(True)
        plt.show()


    min_loss = min(history.history['val_loss'])
    index_min = np.argmin(history.history['val_loss'])
    print(min_loss)
    print(index_min)

    # print("\nEvaluate")
    # result = model.predict(x_test_scaled.values)
    # print("Results: ")
    # print(result)
    # dict(zip(model.metrics_names, result))

    # Plot the history
    #plot_history(history)

    with open('../../data/nn2.csv', 'a', encoding='UTF8', newline='') as f:
        writer = csv.writer(f, delimiter=';')
        row = [b, min_loss, index_min]
        writer.writerow(row)

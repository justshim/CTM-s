import pandas as pd
import tensorflow as tf
from tensorflow.keras import Model
from sklearn.preprocessing import StandardScaler

TEST_DATA_PATH = '../../data/CTM_param_out_A2_15cells.CSV'

print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")

test_data = pd.read_csv(TEST_DATA_PATH, sep=";")
test_data = test_data.drop('t', axis=1)


def scale_datasets(x_test):
    # Standard Scale test and train data
    # Z - Score normalization
    standard_scaler = StandardScaler()
    x_test_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_test),
        columns=x_test.columns
    )
    return x_test_scaled

def input_flatten(df):
    output = pd.DataFrame()
    for index in range(0, len(df.index), 15):
        matrix = df.loc[index:index + 14, :].to_numpy()
        flat = matrix.flatten()
        aa = pd.DataFrame(flat).transpose()
        output = pd.concat([output, aa])
    return output


x_test_scaled = scale_datasets(test_data)
input = input_flatten(x_test_scaled)


#load NN model
model = tf.keras.models.load_model('SavedModel/modellovdue', compile=False)


output = model.predict(input)
print(output)


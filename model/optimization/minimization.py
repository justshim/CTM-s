import pandas as pd
import tensorflow as tf
from tensorflow.keras import Model
from sklearn.preprocessing import StandardScaler

TEST_DATA_PATH = '../../data/opti_data_all_validazione.csv'

TARGET_NAME_1 = 'integral'
TARGET_NAME_2 = 'pi'
print("GPU test: " + str(tf.test.is_gpu_available()) + "\n")

test_data = pd.read_csv(TEST_DATA_PATH)
#test_data = test_data.drop('priority', axis=1)
test_data = test_data.drop('max_delta', axis=1)

x_test = test_data.drop([TARGET_NAME_1, TARGET_NAME_2], axis=1)


def scale_datasets(x_test):
    # Standard Scale test and train data
    # Z - Score normalization
    standard_scaler = StandardScaler()
    x_test_scaled = pd.DataFrame(
        standard_scaler.fit_transform(x_test),
        columns=x_test.columns
    )
    return x_test_scaled


x_test_scaled = scale_datasets(x_test)


#load NN model
model = tf.keras.models.load_model('SavedModel/modellovdue', compile=False)


x_test[['integ_pred', 'pi_pred']] = model.predict(x_test_scaled)
print(x_test)
x_test.to_csv("../../data/pino.csv", index=False, decimal=',', float_format='%g', encoding="utf-8", sep=';')

sol_ammissibili = x_test.loc[(x_test['j'] <= x_test['i']+2)]

min_integ = sol_ammissibili['integ_pred'].min()
max_pi = sol_ammissibili['pi_pred'].max()

print("min integral: " + str(min_integ))
print("max pi: " + str(max_pi))

upperbound_integ = (min_integ*0.2)+min_integ
lowerbound_pi = max_pi-(max_pi*0.1)
print("upperbound integral: " + str(upperbound_integ))
print("lowerbound pi: " + str(lowerbound_pi))
sol_ottime_sciccherie = sol_ammissibili.loc[(sol_ammissibili['integ_pred'] <= upperbound_integ) & (sol_ammissibili['pi_pred'] >= lowerbound_pi)]
print(sol_ottime_sciccherie)


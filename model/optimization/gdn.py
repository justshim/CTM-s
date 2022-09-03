import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel

TRAIN_DATA_PATH = '../../data/opti_data_all.csv'
TEST_DATA_PATH = '../../data/opti_data_all_validazione.csv'

TARGET_NAME_1 = 'integral'
TARGET_NAME_2 = 'pi'

# x_train = features, y_train = target
train_data = pd.read_csv(TRAIN_DATA_PATH)
train_data = train_data.drop(['priority', 'max_delta'], axis=1)

test_data = pd.read_csv(TEST_DATA_PATH)
test_data = test_data.drop(['priority', 'max_delta'], axis=1)

y_train = train_data[[TARGET_NAME_1, TARGET_NAME_2]]
y_train = y_train.head(10000)
y_test = test_data[[TARGET_NAME_1, TARGET_NAME_2]]
y_test = y_test.head(10000)

y_train['pi'] = y_train['pi'] * 100
y_test['pi'] = y_test['pi'] * 100

x_train = train_data.drop([TARGET_NAME_1, TARGET_NAME_2], axis=1)
x_train = x_train.head(10000)
x_test = test_data.drop([TARGET_NAME_1, TARGET_NAME_2], axis=1)
x_test = x_test.head(10000)

X, y = x_train, y_train
print("Sto facendo le robe GDN.")
# kernel = (DotProduct() + WhiteKernel())
gpr = GaussianProcessRegressor(kernel=None, random_state=0).fit(X, y)
score = gpr.score(X, y)
print("Score: " + str(score))
pred = gpr.predict(X.head(10), return_std=True)
print("Pred: " + str(pred))

import pandas as pd
from sklearn.model_selection import train_test_split
from os.path import splitext
import sys


def merge_lid_parts(_lid_filePath, _datasetPath):
    data = pd.read_csv(_datasetPath, header=None)

    df = pd.read_csv(_lid_filePath, header=None)

    df.columns = ['lid', 'oid']
    df = df.sort_values(by=['lid'])
    q = df.quantile([0.00, 0.25, 0.50, 0.75, 1.00])
    col = 'lid'

    q1 = df[((df[col] >= q[col][0.00]) & (df[col] < q[col][0.25]))]
    q2 = df[((df[col] >= q[col][0.25]) & (df[col] < q[col][0.50]))]
    q3 = df[((df[col] >= q[col][0.50]) & (df[col] < q[col][0.75]))]
    q4 = df[((df[col] >= q[col][0.75]) & (df[col] <= q[col][1.00]))]

    train, test = train_test_split(q1, test_size=0.1)
    data_train = data.iloc[train['oid']]
    data_test = data.iloc[test['oid']]
    train_name = splitext(_datasetPath)[0] + '_Q1_train.csv'
    test_name = splitext(_datasetPath)[0] + '_Q1_test.csv'
    data_train.to_csv(train_name, index=False, header=False)
    data_test.to_csv(test_name, index=False, header=False)

    train, test = train_test_split(q2, test_size=0.1)
    data_train = data.iloc[train['oid']]
    data_test = data.iloc[test['oid']]
    train_name = splitext(_datasetPath)[0] + '_Q2_train.csv'
    test_name = splitext(_datasetPath)[0] + '_Q2_test.csv'
    data_train.to_csv(train_name, index=False, header=False)
    data_test.to_csv(test_name, index=False, header=False)

    train, test = train_test_split(q3, test_size=0.1)
    data_train = data.iloc[train['oid']]
    data_test = data.iloc[test['oid']]
    train_name = splitext(_datasetPath)[0] + '_Q3_train.csv'
    test_name = splitext(_datasetPath)[0] + '_Q3_test.csv'
    data_train.to_csv(train_name, index=False, header=False)
    data_test.to_csv(test_name, index=False, header=False)

    train, test = train_test_split(q4, test_size=0.1)
    data_train = data.iloc[train['oid']]
    data_test = data.iloc[test['oid']]
    train_name = splitext(_datasetPath)[0] + '_Q4_train.csv'
    test_name = splitext(_datasetPath)[0] + '_Q4_test.csv'
    data_train.to_csv(train_name, index=False, header=False)
    data_test.to_csv(test_name, index=False, header=False)


if __name__ == '__main__':
    lid_filePath = sys.argv[1]
    datasetPath = sys.argv[2]
    merge_lid_parts(lid_filePath, datasetPath)

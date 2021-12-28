# -*- coding: utf-8 -*-
import pandas
import pandas as pd
import numpy as np 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
def RFclass():
    #名称要对应
    X = pd.read_csv("Methylation_zen.csv",sep=",")
    record = pd.read_csv("a.csv",sep=",")
    Y = pd.read_csv("Vinblastine_Methylation.csv",sep=",")
    Y = Y["LN_IC50"]
    list = []
    for i in range(0,len(record["x"])):
        list.append(record["x"][i]-1)
#    for j in range(0,len(X[X.columns[1]])):
#        f = 0
#        for i in range(0,len(record["x"])):
#            if(j == (record["x"][i]-1)):
#                f= 1
#        if(f == 0):
#            list.append(i)
#    X = X.drop(list)
    print(len(record["x"]))
    #print(list)
    X = X.loc[list,X.columns]
    X = np.transpose(X) 
    X_train,X_test,y_train,y_test = train_test_split(X, Y,test_size=0.2)
    count = 0
    crct = 0
    neigh = RandomForestClassifier(criterion='gini' )
    neigh.fit(X_train, y_train)
    res = neigh.predict(X_test)
    for i in range(len(res)):
        count += 1
        if res[i] == y_test.iloc[i]:
            crct += 1
    number = crct/count
    list=[number]
    nu = pd.DataFrame(list)
    nu.to_csv("number.csv")
    print (crct,' Correct in total : ',count,'\t',crct/count * 100,'%')
    a = neigh.feature_importances_
    return a

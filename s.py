# -*- coding: utf-8 -*-
"""
@author : syq
@date   : 2022/10/1 0:16
"""
import pandas as pd
import numpy as np
import math
from numpy import array

def cal_weight(x):
    """
    @param x:作为计算数据输入时的表格
    @return: 返回值中x代表正向标准化过后的矩阵，w代表权重（weight），d代表冗余（redundancy）
    """
    # standardization
    x = x.apply(lambda x: ((x - np.min(x)) / (np.max(x) - np.min(x))))
    print(x)
    rows = x.index.size
    cols = x.columns.size
    k = 1.0 / math.log(rows)

    # information entropy
    x = array(x)
    D = np.full_like(x, None)

    for i in range(0, rows):
        for j in range(0, cols):
            if x[i][j] == 0:
                Dij = 0.0
            else:
                p = x[i][j] / x.sum(axis=0)[j]
                Dij = math.log(p) * p
            D[i][j] = Dij
    D = pd.DataFrame(D)
    E = D

    # Calculate redundancy
    d = 1 - (-k) * E.sum(axis=0)

    # Calculate the weight of each index
    w = [[None] * 1 for i in range(cols)]
    for j in range(0, cols):
        wj = d[j] / sum(d)
        w[j] = wj
    w = pd.DataFrame(w,columns=['weight'])
    d = pd.DataFrame(d,columns=['redundancy'])

    #return pd.concat([w,d],axis=1)
    return w


if __name__ == '__main__':
    df = pd.read_excel(r"C:\Users\93757\Desktop\data.xlsx").iloc[:,1:]
    result = cal_weight(df)
    raw_name=pd.read_excel(r"C:\Users\93757\Desktop\data.xlsx").iloc[:,0]
    data=pd.read_excel(r"C:\Users\93757\Desktop\data.xlsx").columns
    data=pd.DataFrame(data[1:],columns=['Gene'])

    weight=result["weight"].values.tolist()
    for i in range(len(weight)):
        df.iloc[:,i]=df.iloc[:,i]*weight[i]

    pathway_sum=df.sum(axis=1)
    pathway_sum= pd.DataFrame(pathway_sum,columns=['pathway_sum'])
    data_sum=pd.concat([raw_name,pathway_sum],axis=1)
    data_weight=pd.concat([data,result],axis=1)


    #df['通路活性']=df.sum(axis=1)
    #df["通路活性"].to_excel(r"C:\Users\93757\Desktop\1.xlsx",index=raw_name,encoding="utf_8_sig")
    data_weight.to_excel(r"C:\Users\93757\Desktop\gene_weight.xlsx")
    data_sum.to_excel(r"C:\Users\93757\Desktop\pathway_sum.xlsx",index=False,encoding="utf_8_sig")



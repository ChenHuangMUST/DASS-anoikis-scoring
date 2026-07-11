# DASS anoikis MLP modeling code
from sklearn.neural_network import MLPClassifier
import numpy as np

import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

# %% Cell 2
expr = pd.read_csv("D:/wuyingyi/FNN/ano_score/Anoikis/pancancer_fpkm_mRNA_01A.txt",sep="\t",index_col = 0)

# %% Cell 3
tran_expr = expr.T

# %% Cell 4
ac = pd.read_csv("D:/wuyingyi/FNN/ano_score/Anoikis/ac.txt",sep="\t",index_col = 0)

# %% Cell 5
ac = ac.replace({'ANOact': 1,'ANOint1':1,'ANOint2':2,'ANOint3':2,'ANOint4' :3,'ANOrep': 3})

# %% Cell 6
tag_ac = ac.iloc[:,0]

# %% Cell 7
X_train, X_test, y_train, y_test = train_test_split(tran_expr,tag_ac, test_size=0.3, random_state=42)

# %% Cell 8
gene_select_autoencoder3 = pd.read_csv("D:/wuyingyi/FNN/ano_score/Anoikis/feature_importance.txt",sep="\t",index_col = 0)

# %% Cell 9
selected_genes2 = gene_select_autoencoder3[gene_select_autoencoder3['relative_importance'] >= 0.95].index.tolist()

# %% Cell 10
gene_lst = list(set(selected_genes2) & set(X_test.columns))

# %% Cell 11
X_train_subset = X_train[gene_lst]
X_test_subset = X_test[gene_lst]

# %% Cell 12
param_grid = {
        'hidden_layer_sizes': [(50,),(100,),(200,),(300,),(50,50),(100,100),(200,200),(200,300),(50,50,50),(100,100,100)],
        'max_iter': [5000,10000,150000,500000],
        'learning_rate_init':[0.1,0.01,0.001,0.0001]
}

# %% Cell 13
model = MLPClassifier(hidden_layer_sizes=(100,), activation='tanh', max_iter=50000,learning_rate_init=0.0001)

# %% Cell 14
grid_search = GridSearchCV(model,param_grid=param_grid,cv=5,n_jobs=-1)
grid_search.fit(X_train_subset, y_train)

# %% Cell 15
print("Best parameters found: ",grid_search.best_params_)
print("Best score: ",grid_search.best_score_)

# %% Cell 16
model = MLPClassifier(hidden_layer_sizes=(200,), activation='tanh', max_iter=100000,learning_rate_init=0.0001)

# %% Cell 17
model.fit(X_train_subset, y_train)

# %% Cell 18
train_proba = model.predict_proba(X_train_subset)
test_proba = model.predict_proba(X_test_subset)

# %% Cell 19
train_proba2 = model.predict(X_train_subset)
test_proba2 = model.predict(X_test_subset)

# %% Cell 20
auc = roc_auc_score(y_train, train_proba, multi_class='ovo',average='macro')
print("training auc:", auc)

# %% Cell 21
auc = roc_auc_score(y_test, test_proba, multi_class='ovo',average='macro')
print("testing auc:", auc)

# %% Cell 22
### Output predicted probabilities ####
df = pd.DataFrame (train_proba,index=X_train.index, columns=['S1','S2','S3'])
df.to_csv('D:/wuyingyi/FNN/ano_score/Anoikis/figure/figure5/ae_train_proba.csv',sep="\t")

df = pd.DataFrame (test_proba,index=X_test.index, columns=['S1','S2','S3'])
df.to_csv('D:/wuyingyi/FNN/ano_score/Anoikis/figure/figure5/ae_test_proba.csv',sep="\t")

# %% Cell 23
### Output predicted classes ####
df2 = pd.DataFrame (train_proba2,index=X_train.index)
df2.to_csv('D:/wuyingyi/FNN/ano_score/Anoikis/figure/figure5/ae_train_predict.csv',sep="\t")

df2 = pd.DataFrame (test_proba2,index=X_test.index)
df2.to_csv('D:/wuyingyi/FNN/ano_score/Anoikis/figure/figure5/ae_test_predict.csv',sep="\t")

# %% Cell 24
import joblib

# %% Cell 25
joblib.dump(model, 'D:/wuyingyi/FNN/ano_score/code/model.pkl')


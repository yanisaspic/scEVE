import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

data = pd.read_csv("classification.csv", index_col=0)
cm = pd.crosstab(data.ground, data.label)
fmt = ".2f"
cm = cm.div(cm.sum(axis=1), axis=0)
sns.heatmap(cm, annot=True, fmt=fmt)
plt.show()
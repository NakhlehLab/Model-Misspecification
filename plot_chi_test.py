import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy
df = pd.read_csv("test_chi2_gtest2.txt", names=['chi2', 'p-chi2', 'gtest', 'p-gtest'])
print(df)
fig, axes = plt.subplots(2, 2, figsize=(10,8))
sns.histplot(data=df, x="chi2", stat="probability", ax=axes[0][0])
x = np.arange(0, 400, .001)
axes[0][0].plot(x, scipy.stats.chi2.pdf(x, df=104), color='r', lw=2)
sns.histplot(data=df, x="p-chi2", stat="probability", ax=axes[0][1])
sns.histplot(data=df, x="gtest", stat="probability", ax=axes[1][0])
sns.histplot(data=df, x="p-gtest", stat="probability", ax=axes[1][1])
axes[1][0].plot(x, scipy.stats.chi2.pdf(x, df=104), color='r', lw=2)
plt.savefig("./test_null2.pdf")
plt.show()
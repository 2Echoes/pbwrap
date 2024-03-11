
#### TESTING ###
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pbwrap.plot.superplot as sp

N = 1000
letters = np.array(['type1','type2','wt', 'control'])
exp = np.array(['exp1', 'exp2', 'exp3'])
labels = letters[np.random.randint(0, len(letters), N)]
exp_num = exp[np.random.randint(0, len(exp), N)]


super_df = pd.DataFrame({
    'id' : np.arange(N),
    'lvl3' : exp_num,
    'lvl2' : np.random.randint(0,5,N),
    'lvl1' : labels,
    'data' : np.random.rand(N)
    })

super_df_lvl1 = super_df.groupby(['lvl1'])['data'].apply(list)
super_df_lvl2 = super_df.groupby(['lvl1','lvl2'])['data'].apply(list)
super_df_lvl3 = super_df.groupby(['lvl1','lvl2','lvl3'])['data'].apply(list)

print(super_df)
print(super_df_lvl1)
print(super_df_lvl2)
print(super_df_lvl3)


figsize= (10,10)
fig = plt.figure(figsize=figsize)
ax = fig.gca()
ax = sp.distribution_super_plot(
    super_df_lvl3, 
    ax,
    title= 'SUPERPLOT !',
    xlabel= 'Label',
    ylabel= 'distribution measure'
    )
plt.show()

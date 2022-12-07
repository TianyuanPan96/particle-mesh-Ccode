# %%
import matplotlib.pyplot as plt
%matplotlib inline
print(-1 % 3)
import os
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))

print(cwd, dir_path)
plt.plot([1,2,3], [0,2,5])

# %%

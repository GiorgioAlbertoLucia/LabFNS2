import pandas as pd
import matplotlib.pyplot as plt

# Read CSV file into DataFrame
df = pd.read_csv('TCT/data/input/dati1.csv', comment='#')

# Extract columns of interest
x = df['Attenuation']
df['difference'] = -(df['Asign']-df['Anoise'])
y_err = df['Anoise_err']

df.to_csv('TCT/data/output/dati1_with_difference.csv')

# Plot data with error bars
plt.errorbar(x, df['difference'], yerr=y_err, fmt='o')
plt.xlabel('Attenuation')
plt.ylabel('Difference (pWb)')
plt.show()

   
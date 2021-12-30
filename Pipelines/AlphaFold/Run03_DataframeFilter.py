#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
infile_name = 'C:/Dev/Github/LeucipPipelines/Pipelines/AlphaFold/Csv/pdb_data_geometry.csv'
outfile_name = 'C:/Dev/Github/LeucipPipelines/Pipelines/AlphaFold/Csv/pdb_data_geometry_filter1.csv'
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import pandas as pd


# LOAD ##################################################################
print('load csv',infile_name)
df_geometry = pd.read_csv(infile_name)

# FILTER AS YOU CHOOSE
print('filter csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df_geometry = df_geometry.query('`' + 'C:N+1' + '` >=' + str(2))
#df_geometry = df_geometry.query('Probability >= 90')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# SAVE ##################################################################
print('save csv',infile_name)
df_geometry.to_csv(outfile_name, index=False)
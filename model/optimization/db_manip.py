import pandas as pd

path_data = "../../data/opti_data_i-j-delta-beta_validazione.csv"
path_cells = "../../../../CTMs-identification/fnc/extracted_data/CTM_param_out_validazione.xls"

data = pd.read_csv(path_data)
cells = pd.read_excel(path_cells, sheet_name='Cells parameters')
cells = cells.drop('t', axis=1)
cells = cells.rename(columns={'ID': 'i'})

final = pd.merge(data, cells, on='i', how='left')
cells = cells.rename(columns={'i': 'j'})
final = pd.merge(final, cells, on='j', how='left')

final.to_csv("../../data/opti_data_all_validazione.csv", index=False, float_format='%g', encoding="utf-8")



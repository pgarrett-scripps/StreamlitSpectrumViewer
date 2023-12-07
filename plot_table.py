import pandas as pd
import streamlit as st

df = pd.read_csv(r"C:\Users\Ty\Downloads\2023-12-06T00-57_export.csv")

seq = 'TALLDAAGVASLLTTAEVVVTEIPK'

st.dataframe(df)

forward_ions = 'abc'
reverse_ions = 'xyz'

ion_types = ['b', 'y']
charges = [1, 1]

data = {'AA': list(seq)}
for ion_type, charge in zip(ion_types, charges):
    ion_df = df[(df['ion_type'] == ion_type) & (df['charge'] == charge)]
    if ion_type in forward_ions:
        ion_df = ion_df.sort_values(by='mz')
    else:
        ion_df = ion_df.sort_values(by='mz', ascending=False)
    data['+'*charge + ion_type] = ion_df['mz'].tolist()


df = pd.DataFrame(data)
st.dataframe(df, hide_index=True)
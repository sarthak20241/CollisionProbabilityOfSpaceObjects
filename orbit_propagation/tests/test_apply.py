import pandas as pd

# Sample DataFrame
df = pd.DataFrame({
    'A': [1, 2, 3],
    'B': [4, 5, 6],
    'C': [7, 8, 9]
})

# Define a function that takes a row as input and returns a value
def process_row(row):
    return row['A'] + row['B'] + row['C']

# Apply the function to each row
result = df.apply(process_row, axis=1)

print(result)
import pandas as pd

arr1 = [
	[1,2,3],
	[4,5,6]
]

arr2 = [
	[5,6,7],
	[8,9,10]
]

df1 = pd.DataFrame(arr1, columns=['a', 'b', 'c'])
df2 = pd.DataFrame(arr2, columns=['p','q','r'])

print(df1)
print(df2)

df3 = pd.concat([df1, df2])
df4 = pd.concat([df1,df2],axis=1)

print(df3)
print(df4)
from test_ESPT import run_test
from OUP.Visualize import Visualize

df, mean_df, tim = run_test()
Visualize.sample_point_visualize(df)
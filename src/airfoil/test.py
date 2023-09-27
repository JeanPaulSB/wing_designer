from naca import Naca
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

profile = Naca("2408", 1, 1000)


profile.set_equations()
profile.compute()
df = pd.DataFrame(
    {
        "Alpha (deg)": profile.thinAirfoilTheory((-2, 10))[0],
        "Cl": profile.thinAirfoilTheory((-2, 10))[1],
    }
)
print(df)

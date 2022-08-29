import pandas as pd
import numpy as np

# Import and format data sets

# Aubert 2012 Data

AubertData = pd.read_csv(r'~/PathToData')  # Add path to Aubert Data


AubertData_gran_base = pd.DataFrame(AubertData, columns=['Granulo']).to_numpy()[:, 0]
AubertData_age_base = pd.DataFrame(AubertData, columns=['age']).to_numpy()[:, 0]

AubertData_age = np.delete(AubertData_age_base, np.argwhere(np.isnan(AubertData_gran_base)))
AubertData_gran = np.delete(AubertData_gran_base, np.argwhere(np.isnan(AubertData_gran_base)))

# Andreu 2022 Data

AndreuData = pd.read_csv(r'~/PathToData')  # Add path to Andreu Data

AndreuData_age_base = AndreuData.iloc[:, 1].to_numpy()[8:]
AndreuData_gran_base = AndreuData.iloc[:, 14].to_numpy()[8:]

AndreuData_age = np.delete(AndreuData_age_base, np.argwhere(AndreuData_gran_base == "-"))
AndreuData_gran = np.delete(AndreuData_gran_base, np.argwhere(AndreuData_gran_base == "-"))

AndreuData_gran = AndreuData_gran.astype("float64")
AndreuData_age = AndreuData_age.astype("float64")

# Alder 2018 Data

AlderData_gran = np.array(
    [])

AlderData_age = np.array(
    [])


AlderData_age = np.delete(AlderData_age, np.argwhere(AlderData_gran == None))
AlderData_gran = np.delete(AlderData_gran, np.argwhere(AlderData_gran == None))
AlderData_age = AlderData_age.astype('float64')
AlderData_gran = AlderData_gran.astype('float64')

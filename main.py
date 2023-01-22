from scipy.io import loadmat

# Load .mat data into Python
mat = loadmat('Demanda/BaseDatosNE_31_10_2022New_England.mat')

popG = mat['poperationG']
popL = mat['poperationL']

# Datos de la corrida
opf = mat['OPF']

# Voltage Angle
# Voltage Magnitude
volt_ang = opf['VA']
volt_mag = opf['VM']

# Conversión a matrices

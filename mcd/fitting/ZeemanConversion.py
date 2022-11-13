from scipy import constants as const
import numpy as np
import matplotlib

ev_list = [0.087265039,
0.029256444,
0.043754813,
0.052807003,
0.064074696,
]


B_list = [6.053056286,
0.909465149,
1.55709258,
2.859458022,
4.336181186,
]

e=const.e #charge of electron (C)
m_e=const.m_e #mass of electron (kg)
m_list = []

for ev, B in zip(ev_list,B_list):
    w_c=ev/1000/const.physical_constants['Planck constant over 2 pi in eV s'][0] #cyclotron resonance frequency from Planck constant in eV/Hz
    effective_mass=e*B/w_c/2/m_e/const.pi #effective mass (m*/m_e)
    # print(effective_mass)
    m_list.append(effective_mass)

m_avg = np.mean(m_list)
m_std = np.std(m_list)

print(m_avg, u"\u00B1", m_std)
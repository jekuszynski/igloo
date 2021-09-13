import sys
import numpy as np
from pint.unit import Unit
from scipy import constants as const
from pint import UnitRegistry

def mhz_to_nm(mhz, lambd):
    probe_freq = 1/lambd
    mhz_convert = mhz/const.c/1000
    nm = 1 / (probe_freq - mhz_convert)
    return nm

def ev_to_nm(eV):
    return 1240/eV

def nm_to_ev(nm):
    return 1240/nm

def wavenum_to_nm(wavenum):
    return 1/(wavenum*(1e-7))

ureg = UnitRegistry()
quantity = ureg.Quantity
try: 
    num, start, end = input('What would you like to convert? Format: "# [start_unit] [end_unit]"\n').split()
except ValueError:
    num, start, end = 20, 'm', 'nm'
_quant = quantity(num, str(start))
# print(start_quant)
print('Starting with ' + str(_quant))
old_quant = str(str(num) + ' * ' + start)


## what the fuck 

_quant(old_quant).to(end)
print(_quant)
print('Ending with ' + str(_quant))


sys.exit()


unit1 = str(unit1)
unit2 = str(unit2)

quantity = ureg.Quantity
quantity(num,unit1)

print('Started with ')
quantity
new_quantity = quantity.to(unit2)
print('Ended with ' + new_quantity)

if __name__ == '__main__':
    print(wavenum_to_nm(1.57e4))

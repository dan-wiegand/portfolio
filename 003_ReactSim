import numpy as np
from scipy.integrate import odeint


temp_reactor = 600 # Units: K
pressure_reactor = 2e06 # Units: Pa
volume_reactor = (0.001)*(3) # For volume: (cross-section area)*(height) in units of m^3
gas_constant = 8.3145 #
total_concentration = ((pressure_reactor)*(volume_reactor)) / ((gas_constant)*(temp_reactor))
Ct = (total_concentration/volume_reactor) # mol/m^3
print("Reactor total concentration: " + str(Ct))
print("The assertion will fail here because the np.isclose expecting something way too accurate to what is defined.")
print(str(Ct) + " vs " + str(400.9079860285491))

F_A = (12/720) # Initial feed of cellulose mol/s (converted using MW of cellulose as 20% of slurry)
F_CO = 0 # Initial feed of CO is null as this is one of the products
F_H2 = 0 # Initial feed of H2 is null as this is one of the products
F_H20 = (48/18) # Initial feed of water in mol/s (converted useing MW of cellulose as 40% of slurry)

molar_flow_rates = np.array([F_A, F_CO, F_H2, F_H20])

def fvector(molar_flow_rates, volume):
    """
    A function that returns the volume derivative of `yvector`,
    where `yvector` is a vector with 3 elements.
    """
    yvector = molar_flow_rates

    # Definitions for the rate of cellulose consumption
    k = 10  # s s^-1
    a = 0.001  # mole/m^3

    Ft = np.sum(yvector) # Total molar flow rate
    C_A = Ct * yvector[0]/Ft  # Concentration of cellulose
    C_CO = Ct * yvector[1]/Ft # Concentration of CO

    r_A = -k * C_A / (1 + a * C_CO) # Rate of cellulose comsumption
    r_CO = -24 * r_A # Rate of CO production
    r_H2 = -24 * r_A # Rate of H2 production
    r_H20 = 0 # Rate of H20 consumption/production, this is null

    # Return the array containing the rates, does not include H20
    yvector_derv = np.array([r_A, r_CO, r_H2])

    return print(yvector_derv)

fvector(molar_flow_rates, 1)


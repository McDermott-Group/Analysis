from antennalib import Blackbox as Bbox

Mode = Bbox()

f_A = 80e9
C_A = 7e-15
R_A = 656

f_Q = 4.178e9
C_Q = 58.2e-15
L_Q = 21e-9

AntennaInput = [f_A, C_A, R_A]
QBInput = [f_Q, C_Q, L_Q]

Mode.Input_update(AntennaInput=AntennaInput, QBInput=QBInput)
# print("g/2pi (GHz)=", Mode.g/(2*3.14161e9))
# print("T_Q(Purcell)=", Mode.QB_mode["T1"])
# print("Antenna=", Mode.Antenna_mode)
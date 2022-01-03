import matplotlib.pyplot as plt
import numpy as np

### generate a scatter circle with xy coordinates
n = 361 # number of points
# n = 25+1 # number of points
n_list = np.arange(0, n, 1)
r = 1 # radius
x = r * np.cos(np.linspace(0, 2 * np.pi, n))
y = r * np.sin(np.linspace(0, 2 * np.pi, n))
# print('len(n_list)=', len(n_list))
# print('len(x)=', len(x))


plt.rcParams["figure.figsize"] = (8, 8)
# plt.scatter(x[:10], y[:10])
# plt.show()

### for each xy or n index, get an amplitude
A = 0.15 # amplitude
offset = A * np.sin(np.linspace(0, 2 * np.pi, n))
# plt.plot(offset)

### each xy scatter point's amp with offset
AwithOffset1 = [r-offset[i] for i in range(n)]
AwithOffset2 = [r+offset[i] for i in range(n)]
# plt.plot(AwithOffset)

x1 = []
y1 = []
x2 = []
y2 = []

for i in range(n):
    x1_new = x[i]/r*AwithOffset1[i]
    y1_new = y[i]/r*AwithOffset1[i]
    x1.append(x1_new)
    y1.append(y1_new)
    x2_new = x[i]/r*AwithOffset2[i]
    y2_new = y[i]/r*AwithOffset2[i]
    x2.append(x2_new)
    y2.append(y2_new)


"""Voltage"""
step = 36
length_ref = 0.3
for i in n_list[0::step]:
    if y1[i] >= 0.1 or y1[i] <= -0.1:
        length = np.sqrt((x2[i]-x[i])**2+(y2[i]-y[i])**2)
        scale = (length/length_ref)**(0.4)
        # plt.arrow(x1[i], y1[i], (x2[i]-x1[i])*scale, (y2[i]-y1[i])*scale,
        print('scale=', scale)
        plt.arrow(x[i], y[i], (x2[i]-x[i])*scale, (y2[i]-y[i])*scale,
                   lw=4*scale, head_width=0.02*scale, head_length=0.02*scale)
plt.plot(x, y, '--', color="grey")
plt.plot(x1, y1, '--', color="blue")
plt.plot(x2, y2, '--', color="red")
plt.xlim([-1.5, 1.5])
plt.ylim([-1.5, 1.5])
plt.axis('off')
# plt.savefig('AntennaVoltage.eps', format='eps')
# plt.savefig('AntennaVoltage.pdf', format='pdf', dpi=1200)
# plt.imsave('AntennaVoltage.pdf')
plt.show()

# """Current"""
# step = 36
# start = 6
# end = 6
# length_ref = 0.08
# plt.plot(x, y, '--', color="grey")
# for i in n_list[0::step]:
#     if y1[i] >= 0.1 or y1[i] <= -0.1:
#         length = np.sqrt((x2[i]-x[i])**2+(y2[i]-y[i])**2)
#         scale = (length/length_ref)
#         # plt.arrow(x1[i], y1[i], (x2[i]-x1[i])*scale, (y2[i]-y1[i])*scale,
#         print('scale=', scale)
#         if y[i]>0:
#             plt.arrow(x[i+end], y[i+end], (x[i-start]-x[i+end]), (y[i-start]-y[i+end]),
#                        lw=1.5*scale, head_width=0.02*scale, head_length=0.02*scale)
#         else:
#             plt.arrow(x[i-start], y[i-start], (x[i+end]-x[i-start]), (y[i+end]-y[i-start]),
#                        lw=1.5*scale, head_width=0.02*scale, head_length=0.02*scale)
# plt.plot(x1, y1, '--', color="blue")
# plt.plot(x2, y2, '--', color="red")
# plt.xlim([-1.5, 1.5])
# plt.ylim([-1.5, 1.5])
# plt.axis('off')
# # plt.savefig('AntennaCurrent.eps', format='eps')
# plt.savefig('AntennaCurrent.svg', format='svg')
# # plt.savefig('AntennaCurrent.pdf', format='pdf', dpi=1200)
# plt.show()
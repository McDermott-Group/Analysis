"""
This is for test the analyze_QPTunneling_pomegranate.py robustness against noise.
"""
from .analyze_QPTunneling_pomegranate import *
from .heatedmap import *

t_meas = 0.1  # units ms
t_QP = 10.0  # units ms
p_QP = 1.0 - exp(-t_meas/t_QP)

p_0g_pool = np.linspace(100, 80, 11)/100
p_1e_pool = np.linspace(60, 100, 11)/100
error = np.zeros((len(p_0g_pool), len(p_1e_pool)))

for r in range(len(p_0g_pool)):
    for c in range(len(p_1e_pool)):
        error_1D = []
        readout_fidelity = [p_0g_pool[r], p_1e_pool[c]]
        print('read_fidelity=', readout_fidelity)
        for i in range(100):
            Hidden_Signal = generate_hidden_signal(p_QP=[p_QP, p_QP])
            Observed_Signal = hidden_to_observed_signal(Hidden_Signal, readout_fidelity=readout_fidelity)
            Recovered_Signal = observed_to_recovered_signal(Observed_Signal)
            error_1D.append(((transitions_count(Recovered_Signal) - transitions_count(Hidden_Signal)) /
                             transitions_count(Hidden_Signal)) * 100)
        error[r][c] = average(error_1D).round(decimals=2)

fig, ax = plt.subplots(figsize=(10, 10))
im, cbar = heatmap(error, p_0g_pool, p_1e_pool, ax=ax,
                   cmap="YlGn", cbarlabel="error_2D [f_0g=95%, f_1e=75%]")
texts = annotate_heatmap(im, valfmt="{x:.0f} %")
fig.tight_layout()
plt.show()

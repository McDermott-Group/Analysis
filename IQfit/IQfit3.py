import matplotlib.pyplot as plt
import matplotlib
from lmfit import minimize, Parameters, Parameter, Model
import lmfit
import numpy as np
from dataChest import *

class IQ_fit(object):
    def __init__(self, path, res_freq):  #path without the file name in ""
    
        self.path = path
        self.res_freq = res_freq
    
    def IQ_fit_single_power(self, file_index=1, power=None, Qc_guess=None, Qi_guess=None, att=0, save=False, Singleshotlike=False, CavitySpectroscopylike=False, ind=1, auto_res_freq=True):  #index of the file to be analyzed from the bottom (index starts from 1)
        d = dataChest(self.path)
        dataset_name = d.ls()[0][-file_index]
        
        
        if "SingleShot" in dataset_name or Singleshotlike is True:
            d.openDataset(dataset_name)
            
            dep_variables = d.getVariables()[1]
            for ii in range(0, len(dep_variables)):
                var_name = dep_variables[ii][0]
                if "S" in var_name and "Phase" not in var_name:
                    s_param_type = var_name
            # ind = 1
            
            data=d.getData()
            # print(data)
            S_21_dB = np.asarray(d.getData()[:, 2])[ind:-ind]
            f = np.asarray(d.getData()[:, 0])[ind:-ind]
            phase = np.asarray(d.getData()[:, 1])[ind:-ind]*2.*np.pi/360. 
            phase = phase - phase.mean()
            S_21 = 10**(S_21_dB/20)
            S_21 = S_21
            print("S_21.max()=", S_21.max())
            S_21 = S_21/S_21.max()
            if auto_res_freq is True:
                min_frequency = f[np.argmin(S_21_dB, axis=0)]
                print("I am here!")
            else:
                min_frequency = self.res_freq*1e9
            S_21 = S_21*np.exp(1j*phase)
            
            if Qc_guess is None:
                
                Qc_guess = 2.55e5
                
                if Qi_guess is None:
                    Qi_guess = 5e5
                                
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    #params.add('Q_i', value= Q_i_guess, min=0.80*Q_i_guess, max=1.2*Q_i_guess)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)
                    print("f_0=", params['f_0'].value)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", params['Q_c'].value)
                    print("f_0=", params['f_0'].value)
                    

                    Qi_guess = params['Q_i'].value
                    print("Qi_guess for next", Qi_guess)
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, params['Q_i'].value, params['Q_c'].value)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                else:
                    Q_i = Qi_guess
                    
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", params['Q_c'].value)
                    # print("Q_c=", Q_c)
                    
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, Q_i, params['Q_c'].value)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                    
            else:
                Q_c = Qc_guess
                
                if Qi_guess is None:
                    Qi_guess = 5e5
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", Q_c)
                    
                    Qi_guess = params['Q_i'].value
                    
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, params['Q_i'].value, Q_c)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                else:
                    Q_i = Qi_guess
                    def func(params, f, np_func):
                        A = params['A'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", Q_c)
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, Q_i, Q_c)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                    
                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()

           
        if "CavitySpectroscopy" in dataset_name or CavitySpectroscopylike is True:
            d.openDataset(dataset_name)
            
            dep_variables = d.getVariables()[1]
            for ii in range(0, len(dep_variables)):
                var_name = dep_variables[ii][0]
                if "S" in var_name and "Phase" not in var_name:
                    s_param_type = var_name
            # ind = 500
           
            # if power is None:
                # raise ValueError("Enter the NA power at which you want to do the IQ fit.")
            # else:
                # except Exception as e:
                    # print("This is error:")
           
            data=d.getData()
            correct_rows = np.where(data[:,0]==power)
            first_index = correct_rows[0].min()
            last_index = correct_rows[0].max()
            S_21_dB = np.asarray(data[first_index:last_index, 2])[ind:-ind]
            f=np.asarray(data[first_index:last_index, 1])[ind:-ind]
            phase=np.asarray(data[first_index:last_index, 3]*2.*np.pi/360.)[ind:-ind]
            
            
            
            phase = phase - phase.mean()
            S_21 = 10**(S_21_dB/20)
            S_21 = S_21
            print("S_21.max()=", S_21.max())
            S_21 = S_21/S_21.max()
            if auto_res_freq is True:
                min_frequency = f[np.argmin(S_21_dB, axis=0)]
            else:
                min_frequency = self.res_freq*1e9
            S_21 = S_21*np.exp(1j*phase)
            
            if Qc_guess is None:
                
                Qc_guess = 2.55e5
                
                if Qi_guess is None:
                    Qi_guess = 5e5
                                
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    #params.add('Q_i', value= Q_i_guess, min=0.80*Q_i_guess, max=1.2*Q_i_guess)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", params['Q_c'].value)
                    

                    Qi_guess = params['Q_i'].value
                    print("Qi_guess for next", Qi_guess)
                    
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, params['Q_i'].value, params['Q_c'].value)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                    

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                else:
                    Q_i = Qi_guess
                    
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", params['Q_c'].value)
                    
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, Q_i, params['Q_c'].value)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                    

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                    
            else:
                Q_c = Qc_guess
                
                if Qi_guess is None:
                    Qi_guess = 5e5
                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", Q_c)
                    
                    Qi_guess = params['Q_i'].value
                    
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, params['Q_i'].value, Q_c)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                    

                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                    
                    #############
                    ### Use this block to make good formatted figures for paper
                    #############
                    # s_param_type = "S21"
                    # fig, ax = plt.subplots(dpi=1200)
                    # plt.rcParams['font.size'] = '18'
                    # ax.tick_params(labelsize='medium')
                    # # plt.figure(221)
                    # ax.plot(f/1e9, np.absolute(S_21),".", markersize="9")
                    # ax.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)), linewidth=3)
                    # plt.title(r'$|%s|$' % s_param_type)
                    # plt.xlabel('Frequency (GHz)', fontsize=22)
                    # plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type), fontsize=22)
                    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
                    # plt.tight_layout()
                    # # plt.savefig("S21_amps.pdf")
                    
                    # fig1, ax1 = plt.subplots(dpi=1200)
                    # plt.rcParams['font.size'] = '18'
                    # ax.tick_params(labelsize='medium')
                    # # plt.figure(222)
                    # ax1.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi),".", markersize="9")
                    # ax1.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi), linewidth=3)
                    # plt.title(r'$\angle %s$' % s_param_type)
                    # plt.xlabel('Frequency (GHz)', fontsize=22)
                    # plt.ylabel(r'$\angle %s\ (degrees)$' % s_param_type, fontsize=22)
                    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%.5f'))
                    # plt.tight_layout()
                    # # plt.savefig("S21_phase.pdf")

                    # plt.figure(dpi=1200)
                    # plt.plot(np.real(1./S_21), np.imag(1./S_21),".", markersize="9")
                    # plt.plot(func(params, f, np.real), func(params, f, np.imag), linewidth=3)
                    # plt.title(r"$%s^{-1}$" % s_param_type)
                    # plt.xlabel(r'Re[$%s^{-1}$]' % s_param_type, fontsize=22)
                    # plt.ylabel(r'Im[$%s^{-1}$]' % s_param_type, fontsize=22)
                    # plt.xlim(0.5,3.5)
                    # plt.tight_layout()
                    # # plt.savefig("S21_IQ.pdf")
                    # plt.show()
                else:
                    Q_i = Qi_guess
                    def func(params, f, np_func):
                        A = params['A'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    
                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", Q_c)
                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as fn:
                            np.savetxt(fn,np.c_[(power-att, Q_i, Q_c)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                    
                    
                    s_param_type = "S21"
                    plt.subplot(221)
                    plt.plot(f/1e9, np.absolute(S_21))
                    plt.plot(f/1e9, 1/np.absolute(func(params, f, np.real)+1j*func(params, f, np.imag)))
                    plt.title('|%s|' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('%s/max[%s]' % (s_param_type, s_param_type))

                    plt.subplot(222)
                    plt.plot(f/1e9, 360.*np.angle(1./S_21)/(2.0*np.pi))
                    plt.plot(f/1e9, 360.*np.angle(func(params, f, np.real)+1j*func(params, f, np.imag))/(2.0*np.pi))
                    plt.title('<%s' % s_param_type)
                    plt.xlabel('Frequency (GHz)')
                    plt.ylabel('<%s (degrees)' % s_param_type)

                    plt.subplot(223)
                    plt.plot(np.real(1./S_21), np.imag(1./S_21),".")
                    plt.plot(func(params, f, np.real), func(params, f, np.imag))
                    plt.title("1/%s" % s_param_type)
                    plt.xlabel('Re[1/%s]' % s_param_type)
                    plt.ylabel('Im[1/%s]' % s_param_type)
                    plt.show()
                    
                    
            
    def IQ_fit_multi_power(self, file_index=0, powers=None, Qc_guess=None, Qi_guess=None, att=0, save=False, ind=1, auto_res_freq=True):
        d = dataChest(self.path)
        dataset_name = d.ls()[0][-file_index]
        
        d.openDataset(dataset_name)
            
        dep_variables = d.getVariables()[1]
        for ii in range(0, len(dep_variables)):
            var_name = dep_variables[ii][0]
            if "S" in var_name and "Phase" not in var_name:
                s_param_type = var_name
        # ind = 500
       
        # if power is None:
            # raise ValueError("Enter the NA power at which you want to do the IQ fit.")
        # else:
            # except Exception as e:
                # print("This is error:")
        Qi_fitted = []
        Qc_fitted = []
        Qi_fit_error = []
        Qc_fit_error = []
        if Qc_guess is None:
            Qc_guess = 2.55e5
            if Qi_guess is None:
                Qi_guess = 5e5
                for index, power in enumerate(powers):
                    print(power)
                    data=d.getData()
                    correct_rows = np.where(data[:,0]==power)
                    first_index = correct_rows[0].min()
                    last_index = correct_rows[0].max()
                    S_21_dB = np.asarray(data[first_index:last_index, 2])[ind:-ind]
                    f=np.asarray(data[first_index:last_index, 1])[ind:-ind]
                    phase=np.asarray(data[first_index:last_index, 3]*2.*np.pi/360.)[ind:-ind]
                    # phase = phase % 360

                    phase = phase - phase.mean()
                    S_21 = 10**(S_21_dB/20)
                    S_21 = S_21
                    print("S_21.max()=", S_21.max())
                    S_21 = S_21/S_21.max()
                    if auto_res_freq is True:
                        min_frequency = f[np.argmin(S_21_dB, axis=0)]
                    else:
                        min_frequency = self.res_freq*1e9
                    S_21 = S_21*np.exp(1j*phase)
                    print(min_frequency)


                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    # print("out2.params=", out2.params['Q_i'].valu)

                    # ci, trace = lmfit.conf_interval(mini, out2, sigmas=[1, 2],
                                                    # trace=True, verbose=False)


                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", params['Q_c'].value)
                    print("Q_i_error=", params['Q_i'].stderr)
                    print("Q_c_error=", params['Q_c'].stderr)

                    Qi_guess = params['Q_i'].value
                    print("Qi_guess for next", Qi_guess)

                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as f:
                            if index==0:
                                np.savetxt(f,np.c_[(power-att, params['Q_i'].value, params['Q_c'].value, params['Q_i'].stderr, params['Q_c'].stderr)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c        Q_i_error        Q_c_error")
                            else:
                                np.savetxt(f,np.c_[(power-att, params['Q_i'].value, params['Q_c'].value, params['Q_i'].stderr, params['Q_c'].stderr)], fmt="%10.5e", delimiter=' ')
                        
                    Qi_fitted.append(params['Q_i'].value)
                    Qc_fitted.append(params['Q_c'].value)
                    Qi_fit_error.append(params['Q_i'].stderr)
                    Qc_fit_error.append(params['Q_c'].stderr)
                    # Qc_fitted.append(Q_c)
                    print(Qi_fitted)
                    print(Qc_fitted)  
                    print(Qi_fit_error)
                    print(Qc_fit_error)
                
                plt.figure()
                # plt.plot(powers-att,Qi_fitted,".", color="blue", label="Q_i")
                # plt.plot(powers-att,Qc_fitted,".", color="red", label="Q_c")
                plt.errorbar(powers-att,Qi_fitted, yerr=Qi_fit_error, fmt=".", color="blue", label="Q_i", capsize=2)
                plt.errorbar(powers-att,Qc_fitted, yerr=Qc_fit_error, fmt=".", color="red", label="Q_c", capsize=2)
                plt.title("IQ fit %s GHz" % self.res_freq)
                plt.xlabel("NA source power")
                plt.ylabel("Q")
                plt.yscale("log")
                plt.legend()
                plt.savefig("%s_IQ_multi_power_plot.png" %self.res_freq)
                plt.show()
            else:
                Q_i = Qi_guess
                for index, power in enumerate(powers):
                    print(power)
                    data=d.getData()
                    correct_rows = np.where(data[:,0]==power)
                    first_index = correct_rows[0].min()
                    last_index = correct_rows[0].max()
                    S_21_dB = np.asarray(data[first_index:last_index, 2])[ind:-ind]
                    f=np.asarray(data[first_index:last_index, 1])[ind:-ind]
                    phase=np.asarray(data[first_index:last_index, 3]*2.*np.pi/360.)[ind:-ind]
                    # phase = phase % 360

                    phase = phase - phase.mean()
                    S_21 = 10**(S_21_dB/20)
                    S_21 = S_21
                    print("S_21.max()=", S_21.max())
                    S_21 = S_21/S_21.max()
                    if auto_res_freq is True:
                        min_frequency = f[np.argmin(S_21_dB, axis=0)]
                    else:
                        min_frequency = self.res_freq*1e9
                    S_21 = S_21*np.exp(1j*phase)
                    print(min_frequency)


                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_c = params['Q_c'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_c', value= Qc_guess, min=1e2, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    # print("out2.params=", out2.params['Q_i'].valu)

                    # ci, trace = lmfit.conf_interval(mini, out2, sigmas=[1, 2],
                                                    # trace=True, verbose=False)


                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", params['Q_c'].value)
                    print("Q_c_error=", params['Q_c'].stderr)
                    # print("Q_c=", Q_c)

                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as f:
                            if index==0:
                                np.savetxt(f,np.c_[(power-att, Q_i, params['Q_c'].value, params['Q_c'].stderr)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c        Q_c_error")
                            else:
                                np.savetxt(f,np.c_[(power-att, Q_i, params['Q_c'].value, params['Q_c'].stderr)], fmt="%10.5e", delimiter=' ')
                        
                    Qi_fitted.append(Q_i)
                    Qc_fitted.append(params['Q_c'].value)
                    Qc_fit_error.append(params['Q_c'].stderr)
                    print(Qi_fitted)
                    print(Qc_fitted)  
                    print(Qc_fit_error)
                
                plt.figure()
                # plt.plot(powers-att,Qi_fitted,".", color="blue", label="Q_i")
                # plt.plot(powers-att,Qc_fitted,".", color="red", label="Q_c")
                plt.errorbar(powers-att,Qi_fitted, fmt=".", color="blue", label="Q_i", capsize=2)
                plt.errorbar(powers-att,Qc_fitted, yerr=Qc_fit_error, fmt=".", color="red", label="Q_c", capsize=2)
                plt.title("IQ fit %s GHz" % self.res_freq)
                plt.xlabel("NA source power")
                plt.ylabel("Q")
                plt.yscale("log")
                plt.legend()
                plt.savefig("%s_IQ_multi_power_plot.png" %self.res_freq)
                plt.show()
    
        else:
            Q_c = Qc_guess
            if Qi_guess is None:
                Qi_guess = 5e5
                for index, power in enumerate(powers):
                    print(power)
                    data=d.getData()
                    correct_rows = np.where(data[:,0]==power)
                    first_index = correct_rows[0].min()
                    last_index = correct_rows[0].max()
                    S_21_dB = np.asarray(data[first_index:last_index, 2])[ind:-ind]
                    f=np.asarray(data[first_index:last_index, 1])[ind:-ind]
                    phase=np.asarray(data[first_index:last_index, 3]*2.*np.pi/360.)[ind:-ind]
                    # phase = phase % 360

                    phase = phase - phase.mean()
                    S_21 = 10**(S_21_dB/20)
                    S_21 = S_21
                    print("S_21.max()=", S_21.max())
                    S_21 = S_21/S_21.max()
                    if auto_res_freq is True:
                        min_frequency = f[np.argmin(S_21_dB, axis=0)]
                    else:
                        min_frequency = self.res_freq*1e9
                    S_21 = S_21*np.exp(1j*phase)
                    print(min_frequency)


                    def func(params, f, np_func):
                        A = params['A'].value
                        Q_i = params['Q_i'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('Q_i', value= Qi_guess, min=0.5e3, max=15e6)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    # print("out2.params=", out2.params['Q_i'].valu)

                    # ci, trace = lmfit.conf_interval(mini, out2, sigmas=[1, 2],
                                                    # trace=True, verbose=False)


                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", params['Q_i'].value)
                    print("Q_c=", Q_c)
                    print("Q_i_error", params['Q_i'].stderr)

                    Qi_guess = params['Q_i'].value
                    print("Qi_guess for next", Qi_guess)

                    
                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as f:
                            if index==0:
                                np.savetxt(f,np.c_[(power-att, params['Q_i'].value, Q_c, params['Q_i'].stderr)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c        Q_i_error")
                            else:
                                np.savetxt(f,np.c_[(power-att, params['Q_i'].value, Q_c, params['Q_i'].stderr)], fmt="%10.5e", delimiter=' ')
                        
                    Qi_fitted.append(params['Q_i'].value)
                    Qc_fitted.append(Q_c)
                    Qi_fit_error.append(params['Q_i'].stderr)
                    print(Qi_fitted)
                    print(Qc_fitted)
                    print(Qi_fit_error)
                
                plt.figure()
                # plt.plot(powers-att,Qi_fitted,".", color="blue", label="Q_i")
                # plt.plot(powers-att,Qc_fitted,".", color="red", label="Q_c")
                plt.errorbar(powers-att,Qi_fitted, yerr=Qi_fit_error, fmt=".", color="blue", label="Q_i", capsize=2)
                plt.errorbar(powers-att,Qc_fitted, fmt=".", color="red", label="Q_c", capsize=2)
                plt.title("IQ fit %s GHz" % self.res_freq)
                plt.xlabel("NA source power")
                plt.ylabel("Q")
                plt.yscale("log")
                plt.legend()
                plt.savefig("%s_IQ_multi_power_plot.png" %self.res_freq)
                plt.show()
    
            else:
                Q_i = Qi_guess
                for index, power in enumerate(powers):
                    print(power)
                    data=d.getData()
                    correct_rows = np.where(data[:,0]==power)
                    first_index = correct_rows[0].min()
                    last_index = correct_rows[0].max()
                    S_21_dB = np.asarray(data[first_index:last_index, 2])[ind:-ind]
                    f=np.asarray(data[first_index:last_index, 1])[ind:-ind]
                    phase=np.asarray(data[first_index:last_index, 3]*2.*np.pi/360.)[ind:-ind]
                    # phase = phase % 360

                    phase = phase - phase.mean()
                    S_21 = 10**(S_21_dB/20)
                    S_21 = S_21
                    print("S_21.max()=", S_21.max())
                    S_21 = S_21/S_21.max()
                    if auto_res_freq is True:
                        min_frequency = f[np.argmin(S_21_dB, axis=0)]
                    else:
                        min_frequency = self.res_freq*1e9
                    S_21 = S_21*np.exp(1j*phase)
                    print(min_frequency)


                    def func(params, f, np_func):
                        A = params['A'].value
                        phi = params['phi'].value
                        phi_0 = params['phi_0'].value
                        f_0 = params['f_0'].value
                        phi_v = params['phi_v'].value
                        dx = (f-f_0)/f_0  
                        return np_func(A*(1.0 + (Q_i/Q_c)*np.exp(1j*phi_0)*(1./(1.+2.*1j*Q_i*dx)))* np.exp(-1j*(dx*f_0*phi_v+phi))) #(f-f_0)

                    def objective(params):
                        resid = np.sqrt((np.real(1./S_21) - func(params, f, np.real))**2+(np.imag(1./S_21) - func(params, f, np.imag))**2)
                        return resid
                        
                    df = f[-1]-f[0]
                    phi_v_guess = ((phase[-1]-phase[0])/df)#*min_frequency
                    print("phi_v_guess=", phi_v_guess)
                    alpha_guess = (np.mean(np.absolute(S_21)[-20:-1])-np.mean(np.absolute(S_21)[0:20]))/(np.mean(f[-20:-1])-np.mean(f[0:20]))
                    print("alpha_guess=", alpha_guess)
                    phi_guess = phase[np.argmin(S_21_dB, axis=0)]
                    print("phi_guess=", phi_guess)
                    
                    # create a set of Parameters
                    params = Parameters()
                    params.add('A', value= 1.0, min=.1, max=1.5)
                    params.add('f_0', value= min_frequency, min=min_frequency-df/20., max=min_frequency+df/20.)
                    params.add('phi_v', value= phi_v_guess, min=1000.*phi_v_guess, max=phi_v_guess*.001)
                    params.add('phi_0', value= 0, min=-2*np.pi, max=2*np.pi)
                    params.add('phi', value= phi_guess, min=-2*np.pi, max=2*np.pi)



                    # run the global fit to all the data sets
                    mini = lmfit.Minimizer(objective, params, nan_policy='omit')
                    out1 = mini.minimize(method='powell') #powell cobyla

                    params = out1.params

                    # then solve with Levenberg-Marquardt using the
                    # Nelder-Mead solution as a starting point
                    out2 = mini.minimize(method='emcee', params=out1.params) #emcee least_squares


                    lmfit.report_fit(out2.params, min_correl=0.5)
                    # print("out2.params=", out2.params['Q_i'].valu)

                    # ci, trace = lmfit.conf_interval(mini, out2, sigmas=[1, 2],
                                                    # trace=True, verbose=False)


                    params = out2.params

                    param_names = ['A', 'Q_i', 'Q_c', 'f_0', 'phi_v', 'phi_0', 'phi_0']

                    print("Q_i=", Q_i)
                    print("Q_c=", Q_c)

                    if save is True:
                        fname = "IQ_data" + str(self.res_freq) + ".txt"        
                        with open(fname,"a") as f:
                            if index==0:
                                np.savetxt(f,np.c_[(power-att, Q_i, Q_c)], fmt="%10.5e", delimiter=' ', header="power        Q_i        Q_c")
                            else:
                                np.savetxt(f,np.c_[(power-att, Q_i, Q_c)], fmt="%10.5e", delimiter=' ')
                        
                    Qi_fitted.append(Q_i)
                    Qc_fitted.append(Q_c)
                    print(Qi_fitted)
                    print(Qc_fitted)  
                
                plt.figure()
                plt.plot(powers-att,Qi_fitted,".", color="blue", label="Q_i")
                plt.plot(powers-att,Qc_fitted,".", color="red", label="Q_c")
                plt.title("IQ fit %s GHz" % self.res_freq)
                plt.xlabel("NA source power")
                plt.ylabel("Q")
                plt.yscale("log")
                plt.legend()
                plt.savefig("%s_IQ_multi_power_plot.png" %self.res_freq)
                plt.show()
    
    # def photons_stored_in_resonator(self, file_index = 1, IQdata_filepath, att_outsidefridge = 0, att_infridge = 0, save=False):
        # E_stored = []
        # powers = np.loadtxt(IQdata_filepath, usecols=[0])
        # Qi = np.loadtxt(IQdata_filepath, usecols=[1])
        # Qc = np.loadtxt(IQdata_filepath, usecols=[2])
        
        # i = 0
        # for index, power in enumerate(powers):
            # print(self.res_freq)
            # print(power)
            # d = dataChest(self.path)
            # dataset_name = d.ls()[0][-file_index]
            # d.openDataset(dataset_name)

            # dep_variables = d.getVariables()[1]
            # for ii in range(0, len(dep_variables)):
                # var_name = dep_variables[ii][0]
                # if "S" in var_name and "Phase" not in var_name:
                    # s_param_type = var_name
            # ind = 300

            # data=d.getData()
            # correct_rows = np.where(data[:,0]==power+att_outsidefridge)
            # first_index = correct_rows[0].min()
            # last_index = correct_rows[0].max()
            # S_21_dB = np.asarray(data[first_index:last_index, 3])[ind:-ind]

            # S_21_dB_High = np.average(S_21_dB[0:200])
            # S_21_db_Low = S_21_dB[np.argmin(S_21_dB)]

            # S_21_dB_dip = S_21_db_Low - S_21_dB_High
            # print("S_21_dB_dip=",S_21_dB_dip)
            
            # E = (2*Qi[i]**2)*(10**(((power+att_outsidefridge) - att_outsidefridge - att_infridge + S_21_dB_dip)/10)*1e-3)/(2*np.pi*self.res_freq*1e9*Qc[i])
            # print("Q_i=",Qi[i])
            # print("Q_c=",Qc[i])
            # print("E=",E)
            # E_stored.append(E)

            # n_photons = E/(6.626e-34 * self.res_freq*1e9)
            # print("n_photons",n_photons)

            # fname = "photon_data" + str(self.res_freq) + ".txt"     
            # with open(fname, "a") as f:
                # if index==0:
                    # np.savetxt(f, np.c_[(power, Qi[i], Qc[i], E, n_photons)], fmt="%1.4e", delimiter=' ', header="power        Q_i        Q_c        Energy_stored        photons")
                # else:
                    # np.savetxt(f, np.c_[(power, Qi[i], Qc[i], E, n_photons)], fmt="%1.4e", delimiter=' ')
            
            # i = i + 1
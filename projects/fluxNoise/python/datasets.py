template = {
              'Q': 'Q2',
           'date': '05-05-20',
       'files_P1': range(8,168+1),
       'files_T1': [],
       'files_QP': [],
   'files_charge': [],
      'thresh_P1': 0.05,
      'thresh_T1': 0.44,
    'path_charge': None
}

q4_0430_T1 = {
              'Q': 'Q4',
           'date': '04-30-20',
       'files_P1': range(0,231+1),
       'files_T1': range(0,231+1),
    'recalibrate': True,
      'thresh_P1': 0.031 
}
q3_0501_T1 = {
              'Q': 'Q3',
           'date': ['05-01-20','05-02-20'],
       'files_P1': [range(666,686+1),range(0,196+1)],
       'files_T1': [range(1,21+1),range(0,196+1)],
      'thresh_P1': 0.017,
      # 'thresh_T1': 0.18,
    'recalibrate': True
}
q1_0503_T1 = {
              'Q': 'Q1',
           'date': '05-03-20',
       'files_P1': range(0,273+1),
       'files_T1': range(0,273+1),
      'thresh_P1': 0.039 
}
q2_0505_T1 = {
              'Q': 'Q2',
           'date': '05-05-20',
       'files_P1': range(8,335+1),
       'files_T1': range(4,331+1),
      'thresh_P1': 0.065 
}
    


q3_0428_QP = {
              'Q': 'Q3',
           'date': ['04-28-20','04-29-20'],
       'files_P1': [range(0,821+1),range(0,110+1)],
       'files_QP': [range(0,821+1),range(0,110+1)],
      'thresh_P1': 0.05,
    'path_charge': None
}
q4_0430_QP = {
              'Q': 'Q4',
           'date': '04-30-20',
       'files_P1': range(234,470+1),
       'files_QP': range(1,237+1),
      'thresh_P1': 0.05,
    'path_charge': 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q4\General\Parameter\cww1911obk_charge_offset.hdf5'
}
q3_0501_QP = {
              'Q': 'Q3',
           'date': '05-01-20',
       'files_P1': range(0,664+1),
       'files_QP': range(0,663+1),
      'thresh_P1': 0.05,
    'path_charge': 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q3\General\Parameter\cwx0615mbj_charge_offset.hdf5'
}
q1_0504_QP = {
              'Q': 'Q1',
           'date': '05-04-20',
       'files_P1': range(0,552+1),
       'files_QP': range(0,552+1),
      'thresh_P1': 0.05,
    'path_charge': 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1\General\Parameter\cxa1458apn_charge_offset.hdf5'
}


# # # # # # # # # # # # # # # # # # #
######### Q1-Q2-Q3-Q4 charge ##########
# # # # # # # # # # # # # # # # # # #
for _ in [1]:

    q1q2q3q4_charge_1 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
       'exclude_data': {'Q2':[(None,None)]},
          'path_base': 'fluxNoise',
        'path_charge': 'cvc0232imo'
    }

    q1q2q3q4_charge_2 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
       'exclude_data': {'Q4':[(None,None)]},
          'path_base': 'fluxNoise',
        'path_charge': 'cvd0430bxy'
    }

    q1q2q3q4_charge_3 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
       'exclude_data': {'Q4':[(None,None)]},
          'path_base': 'fluxNoise',
        'path_charge': 'cvd0802ltu'
    }

    q1q2q3q4_charge_4 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise',
        'path_charge': 'cvd1745asc'
    }

    q1q2q3q4_charge_5 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise',
        'path_charge': 'cve0104ppt'
    }

    q1q2q3q4_charge_6 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise',
        'path_charge': 'cvf0542fei'
    }

    q1q2q3q4_charge_7 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise2',
        'path_charge': 'cvj0203nun'
    }

    q1q2q3q4_charge_8 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise2',
        'path_charge': 'cvj1806ywn'
    }

    q1q2q3q4_charge_9 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise2',
        'path_charge': 'cvk0412mns'
    }

    q1q2q3q4_charge_10 = { # maybe remove bad part in Q4?
                  'Q': 'Q1Q2Q3Q4',
               'date': [''],
          'path_base': 'fluxNoise2',
        'path_charge': 'cvm0651mhe'
    }

    q1q2q3q4_charge_11 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': ['08-23-20','08-24-20'],
          'path_base': 'fluxNoise2',
        'path_charge': 'dbh2053vac'
    }

    q1q2q3q4_charge_12 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': ['08-24-20','08-25-20'],
          'path_base': 'fluxNoise2',
        'path_charge': 'dbi2358fmz'
    }

    q1q2q3q4_charge_13 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': ['08-25-20','08-26-20'],
          'path_base': 'fluxNoise2',
        'path_charge': 'dbj1930bdn'
    }

    q1q2q3q4_charge_14 = {
                  'Q': 'Q1Q2Q3Q4',
               'date': ['08-26-20','08-27-20'],
          'path_base': 'fluxNoise2',
       'exclude_data': {'Q4':[(None,None)]},
        'path_charge': 'dbl0020qwm'
    }

# # # # # # # # # # # # # # # # # # #
######### Q1 Charge, Q2 QP ##########
# # # # # # # # # # # # # # # # # # #
for _ in [1]:

    # # after T1 measurements, starting 6/13
    # cyo0514vli, cyp0520ynw, cyw0458jbl

    # 500 us duty cycle, calibration before only 
    for _ in [1]:

        q1q2_0701_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-01-20','07-02-20'],
          'thresh_charge': 0.98,
           'exclude_data': [(40,55),(148,163),(246,None)],
            'path_charge': 'czh0451zqa'
        }

        q1q2_0702_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-02-20','07-03-20'],
          'thresh_charge': 0.98,
           'exclude_data': [(69,107),(232,283),(363,383)],
            'path_charge': 'czi0122dsm'
        }

    # 500 us duty cycle, calibration before and after 
    for _ in [1]:

        q1q2_0703_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-03-20','07-04-20'],
          'thresh_charge': 0.98,
           'exclude_data': [(40,66),(118,143),(254,267),(315,327)],
            'path_charge': 'czi2109vra'
        }

        q1q2_0704_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-04-20','07-05-20','07-06-20'],
          'thresh_charge': 0.98,
           'exclude_data': [(71,91),(112,127),(137,173),(455,466),(504,519),(484,None)], # files got corrupted at end somehow
            'path_charge': 'czj1344rqj'
        }

        q1q2_0706_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-06-20'],
          'thresh_charge': 0.98,
           'exclude_data': [(25,38),(76,89),(187,None)],
            'path_charge': 'czl0540pcx'
        }

    # 30us duty cycle, calibration before and after
    for _ in [1]:
        q1q2_0707_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-07-20','07-08-20'],
          'thresh_charge': 0.94,
           'exclude_data': [(41,77),(126,150),(227,239),(264,289),(323,335),(347,383),(407,432),(520,532)],
            'path_charge': 'czn0345qcp'
        }

        q1q2_0708_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-08-20','07-09-20'],
          'thresh_charge': 0.94,
           'exclude_data': [(54,76),(110,123),(134,156),(167,179),(250,267),(343,356),(416,None)],
            'path_charge': 'czo0442fuk'
        }

        q1q2_0709_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-09-20','07-10-20'],
          'thresh_charge': 0.94,
           'exclude_data': [],
            'path_charge': 'czo2157jnz'
        }

        q1q2_0710_charge_QP = {
                      'Q': 'Q1Q2',
                   'date': ['07-10-20'],
          'thresh_charge': 0.94,
           'exclude_data': [],
            'path_charge': 'czp1936gto'
        }

        q1q2_0710_charge_QP_2 = { # first no mem errors
                      'Q': 'Q1Q2',
                   'date': ['07-10-20','07-11-20'],
          'thresh_charge': 0.94,
           'exclude_data': [],
            'path_charge': 'czq0101nen'
        }

# # # # # # # # # # # # # # # # # # #
######### Q1 Charge, Q2 T1 ##########
# # # # # # # # # # # # # # # # # # #
for _ in [1]:

    ######## 200us duty cycle, measurement delay 3us ########
    for _ in [1]:

        # cxs2305ngo 5/22, 5/23

        q1q2_0524_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-24-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxu0723thw'
        }

        q1q2_0524_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['05-24-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxu1901pfy'
        }

        q1q2_0525_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-25-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxv0747mel'
        }

        q1q2_0526_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-26-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxw1630xbt'
        }

        q1q2_0526_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['05-26-20','05-27-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxx0332bgf'
        }

        q1q2_0527_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-27-20'],
          'thresh_charge': 0.91,
           'exclude_data': [(915,None)],
            'path_charge': 'cxx1757kfx'
        }

        q1q2_0528_charge_T1 = { # shitty data set, blobs looked wrong
                      'Q': 'Q1Q2',
                   'date': ['05-28-20'],
          'thresh_charge': 0.4,
            'path_charge': 'cxy0502hxc'
        }

        q1q2_0528_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['05-28-20'],
          'thresh_charge': 0.91,
           'exclude_data': [(905,None)],
            'path_charge': 'cxy2033hub'
        }


    ######## 500us duty cycle, measurement delay 3us ########
    for _ in [1]:

        q1q2_0528_charge_T1_3 = {
                      'Q': 'Q1Q2',
                   'date': ['05-28-20','05-29-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxz0428jaq'
        }

        q1q2_0529_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-29-20','05-30-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cxz2301mxf'
        }

        q1q2_0530_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['05-30-20','05-31-20','06-01-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cyb0359hkw'
        }

        q1q2_0603_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['06-03-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cye0625zge'
        }


    ######## 500us duty cycle, measurement delay 9us ########
    for _ in [1]:

        q1q2_0603_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['06-03-20','06-04-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cye2324ghg'
        }

        q1q2_0604_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['06-04-20','06-05-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cyf2353zve'
        }

        q1q2_0605_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['06-05-20','06-06-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cyg2351lgx'
        }

        q1q2_0607_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['06-07-20','06-08-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cyi0501krn'
        }

        q1q2_0608_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['06-08-20','06-09-20'],
          'thresh_charge': 0.91,
            'path_charge': 'cyj0526hzn'
        }


    ######## 30us duty cycle, measurement delay 10us ########
    for _ in [1]:

        q1q2_0711_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-11-20','07-12-20','07-13-20','07-14-20'],
          'thresh_charge': 0.96,
            'path_charge': 'czq0701sxa'
        }

        q1q2_0714_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-14-20','07-15-20'],
          'thresh_charge': 0.96,
           'exclude_data': [(2256,None)],
            'path_charge': 'czt2252myq'
        }

        q1q2_0715_charge_T1 = { # might be able to get something out of the bad data?
                      'Q': 'Q1Q2',
                   'date': ['07-15-20','07-16-20'],
          'thresh_charge': 0.96,
           'exclude_data': [(150,None)],
            'path_charge': 'czu2317rrd'
        }

        q1q2_0719_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-19-20','07-20-20'],
          'thresh_charge': 0.96,
            'path_charge': 'czy1947dka'
        }

        q1q2_0720_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-20-20','07-21-20'],
          'thresh_charge': 0.96,
            'path_charge': 'daa0259pwl'
        }

        q1q2_0722_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-22-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dab0542etb'
        }

        q1q2_0722_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['07-22-20','07-23-20','07-24-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dab2110niu'
        }

        q1q2_0724_charge_T1 = {
                      'Q': 'Q1Q2',
                   'date': ['07-24-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dad0805sah'
        }

        q1q2_0724_charge_T1_2 = {
                      'Q': 'Q1Q2',
                   'date': ['07-24-20','07-25-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dae0110ldi'
        }

        q1q2_0725_charge_T1 = { # no events
                      'Q': 'Q1Q2',
                   'date': ['07-25-20','07-26-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dae2218buk'
        }


        #### adding a third qubit across the chip ####

        q1q2q4_0726_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['07-26-20','07-27-20'],
          'thresh_charge': 0.96,
            'path_charge': 'daf2131pel'
        }

        q1q2q4_0727_charge_T1 = { # no events
                      'Q': 'Q1Q2Q4',
                   'date': ['07-27-20','07-28-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dag1517jng'
        }

        q1q2q4_0728_charge_T1 = { # no events
                      'Q': 'Q1Q2Q4',
                   'date': ['07-28-20','07-29-20'],
          'thresh_charge': 0.98,
            'path_charge': 'dah2039nen'
        }

        q1q2q4_0729_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['07-29-20','07-30-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dai0722bug'
        }

        q1q2q4_0730_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['07-30-20','07-31-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dak0101aiu'
        }

        q1q2q4_0731_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['07-31-20','08-01-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dak2054kks'
        }

        q1q2q4_0801_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-01-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dal1542vaz'
        }

        q1q2q4_0802_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-02-20','08-03-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dam0513qzu'
        }

        q1q2q4_0810_charge_T1 = { # no events
                      'Q': 'Q1Q2Q4',
                   'date': ['08-10-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dau0741dpr'
        }

        q1q2q4_0811_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-11-20','08-12-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dav2150jrn'
        }

        q1q2q4_0812_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-12-20','08-13-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dax0031ufg'
        }

        q1q2q4_0815_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-15-20'],
          'thresh_charge': 0.96,
            'path_charge': 'daz2351ozl'
        }

        q1q2q4_0816_charge_T1 = {
                      'Q': 'Q1Q2Q4',
                   'date': ['08-16-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dba1050uzq'
        }


        #### going back to just two qubits: Q1 - Q4 ####

        q1q4_0803_charge_T1 = {
                      'Q': 'Q1Q4',
                   'date': ['08-03-20','08-04-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dan2149qxx'
        }

        q1q4_0805_charge_T1 = {
                      'Q': 'Q1Q4',
                   'date': ['08-05-20','08-06-20'],
          'thresh_charge': 0.96,
            'path_charge': 'dap0636vpd'
        }
        
        # # # # P1 measurement # # # # #

        q1q2_0806_charge_P1 = {
                      'Q': 'Q1Q2',
                   'date': ['08-06-20','08-07-20'],
          'thresh_charge': 0.96,
            'path_charge': 'daq2130haq'
        }
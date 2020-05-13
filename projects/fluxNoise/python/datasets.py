template = {
              'Q': 'Q2',
           'date': '05-05-20',
       'files_P1': range(8,168+1),
       'files_T1': [],
       'files_QP': [],
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

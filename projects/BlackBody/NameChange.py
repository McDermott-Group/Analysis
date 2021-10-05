# importing os module
import os
import copy


# Function to rename multiple files

base_path = "Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorT1P1_2021Aug17"
base_path_backslash = "Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorT1P1_2021Aug17\\"
os.chdir(base_path)

# for count, filename in enumerate(os.listdir('.')):
#     print('filename=', filename)
#     if "0mDAC" in filename and "_0mDAC" not in filename:
#         matlab_path=base_path_backslash+filename+"\MATLABData"
#         os.chdir(matlab_path)
#         for count_mat, filename_mat in enumerate(os.listdir('.')):
#             filename_mat_new = copy.deepcopy(filename_mat)
#             # print('filename_new=', filename_new)
#             # print('filename=', filename)
#             filename_mat_new = filename_mat_new.replace("mDAC", "0mDAC")
#             os.rename(filename_mat, filename_mat_new)
#
#         HDF5_path=base_path_backslash+filename+"\HDF5Data"
#         os.chdir(HDF5_path)
#         for count_HDF5, filename_HDF5 in enumerate(os.listdir('.')):
#             filename_HDF5_new = copy.deepcopy(filename_HDF5)
#             # print('filename_new=', filename_new)
#             # print('filename=', filename)
#             filename_HDF5_new = filename_HDF5_new.replace("mDAC", "0mDAC")
#             os.rename(filename_HDF5, filename_HDF5_new)


        # filename_new = copy.deepcopy(filename)
        # # print('filename_new=', filename_new)
        # print('filename=', filename)
        # filename_new=filename_new.replace("mDAC", "0mDAC")
        # print('filename_new=', filename_new)
        # print('filename=', filename)
        # os.rename(filename, filename_new)
        ### Go to matlab file change name
        ### Go to hdf5 file change name


for count, filename in enumerate(os.listdir('.')):
    # if "mDAC" in filename and "_0mDAC" not in filename and "00mDAC" not in filename:
    if "mDAC" in filename and "_0mDAC" not in filename:
        filename_new = copy.deepcopy(filename)
        # print('filename_new=', filename_new)
        # print('filename=', filename)
        filename_new=filename_new.replace("mDAC", "0mDAC")
        # print('filename_new=', filename_new)
        # print('filename=', filename)
        os.rename(filename, filename_new)
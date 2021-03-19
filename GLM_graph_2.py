'''This script reads the .txt file generated by JesÃºs Lopez's code extracting andrepresenting energy data'''import matplotlib.pyplot as pltimport numpy as npfile_path = "/Users/jaimemorandominguez/Desktop/TFG/Pruebas_GLM/Salida/GLM_output_data.txt"# Reading the .txt file and writing it in variablesGLM_data = np.loadtxt(file_path, delimiter=' ')# 1st column: flash time [s]# 2nd column: Event latitude [deg]# 3rd column: Event longitude [deg]# 4th column: Event ID# 5th column: Flash latitude [deg]# 6th column: Flash longitude [deg]# 7th column: Radiance########### REPRESENTATIONS ############ Integrated event grid and flash representationplt.figure()plt.scatter(GLM_data[:,2],GLM_data[:,1])plt.scatter(GLM_data[:,5],GLM_data[:,4],marker='x')plt.axis('equal')plt.grid('on')plt.title("Event grid")plt.xlabel('Longitude [deg]')plt.ylabel('Latitude [deg]')plt.show()# Radiance vs time graph representationplt.figure()plt.bar(GLM_data[:,0],GLM_data[:,6],width=0.005)plt.grid('on')plt.title("Radiance VS Time")plt.xlabel('Time (second of the day) [s]')plt.ylabel('Radiance [J]')plt.show()
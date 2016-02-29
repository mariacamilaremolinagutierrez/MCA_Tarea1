import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from matplotlib.backends.backend_pdf import PdfPages

opt = int(sys.argv[1]) 
#Las opciones son 9 10 11 12 13 14 15 16 17

if(opt==9):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_9_simpletica = plt.figure()
    fig_9_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-3.0,3.0)
    plt.ylim(-2.5,2.5)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)

    datos_rk4 = sys.argv[3]
    nombreSinDat_rk4=datos_rk4.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_rk4 = nombreSinDat_rk4.split('/')
    arregloStrings_rk4 = arreglo_strings_grande_rk4[1].split('_')
    a_rk4=(arregloStrings_rk4[2].strip('.dat'))

    data_rk4 = np.loadtxt(datos_rk4)

    fig_9_rk4 = plt.figure()
    fig_9_rk4.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_rk4[1]+' con a: '+a_rk4, fontsize=30)
    plt.xlim(-3.0,3.0)
    plt.ylim(-2.5,2.5)
    plt.scatter(data_rk4[:,1],data_rk4[:,2], s=1.5)

    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_9"+'.pdf')
    pp.savefig(fig_9_simpletica)
    pp.savefig(fig_9_rk4)
    pp.close()

if(opt==10):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_10_simpletica = plt.figure()
    fig_10_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-2.5,-0.8)
    plt.ylim(-0.8,0.8)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)

    datos_rk4 = sys.argv[3]
    nombreSinDat_rk4=datos_rk4.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_rk4 = nombreSinDat_rk4.split('/')
    arregloStrings_rk4 = arreglo_strings_grande_rk4[1].split('_')
    a_rk4=(arregloStrings_rk4[2].strip('.dat'))

    data_rk4 = np.loadtxt(datos_rk4)

    fig_10_rk4 = plt.figure()
    fig_10_rk4.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_rk4[1]+' con a: '+a_rk4, fontsize=30)
    plt.xlim(-2.9,-0.7)
    plt.ylim(-1.1,1.5)
    plt.scatter(data_rk4[:,1],data_rk4[:,2], s=1.5)

    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_10"+'.pdf')
    pp.savefig(fig_10_simpletica)
    pp.savefig(fig_10_rk4)
    pp.close()

if(opt==11):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_11_simpletica = plt.figure()
    fig_11_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-0.25,0.8)
    plt.ylim(-0.4,0.4)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)

    datos_rk4 = sys.argv[3]
    nombreSinDat_rk4=datos_rk4.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_rk4 = nombreSinDat_rk4.split('/')
    arregloStrings_rk4 = arreglo_strings_grande_rk4[1].split('_')
    a_rk4=(arregloStrings_rk4[2].strip('.dat'))

    data_rk4 = np.loadtxt(datos_rk4)

    fig_11_rk4 = plt.figure()
    fig_11_rk4.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_rk4[1]+' con a: '+a_rk4, fontsize=30)
    plt.xlim(-0.25,0.8)
    plt.ylim(-0.4,0.4)
    plt.scatter(data_rk4[:,1],data_rk4[:,2], s=1.5)

    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_11"+'.pdf')
    pp.savefig(fig_11_simpletica)
    pp.savefig(fig_11_rk4)
    pp.close()

if(opt==12):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_12_simpletica = plt.figure()
    fig_12_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-2.0,2.0)
    plt.ylim(-2.5,2.5)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)


    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_12"+'.pdf')
    pp.savefig(fig_12_simpletica)
    pp.close()

if(opt==13):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_13_simpletica = plt.figure()
    fig_13_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(0.2,0.7)
    plt.ylim(-0.7,0.7)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)


    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_13"+'.pdf')
    pp.savefig(fig_13_simpletica)
    pp.close()

if(opt==14):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_14_simpletica = plt.figure()
    fig_14_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-0.08,0.08)
    plt.ylim(-0.1,0.1)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)


    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_14"+'.pdf')
    pp.savefig(fig_14_simpletica)
    pp.close()

if(opt==15):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_15_simpletica = plt.figure()
    fig_15_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-3.0,3.0)
    plt.ylim(-2.5,2.5)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)


    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_15"+'.pdf')
    pp.savefig(fig_15_simpletica)
    pp.close()

if(opt==16):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_16_simpletica = plt.figure()
    fig_16_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('q_3', fontsize=20)
    plt.ylabel('p_3', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    plt.xlim(-3.0,3.0)
    plt.ylim(-2.5,2.5)
    plt.scatter(data_simpl[:,1],data_simpl[:,2], s=1.5)


    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_16"+'.pdf')
    pp.savefig(fig_16_simpletica)
    pp.close()

if(opt==17):

    datos_simpl = sys.argv[2]
    nombreSinDat_simpl=datos_simpl.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_simpl = nombreSinDat_simpl.split('/')
    arregloStrings_simpl = arreglo_strings_grande_simpl[1].split('_')
    a_simpl=(arregloStrings_simpl[2].strip('.dat'))

    data_simpl = np.loadtxt(datos_simpl)
    #t q p energia

    fig_17_simpletica = plt.figure()
    fig_17_simpletica.set_size_inches(14.5,10.5)
    plt.xlabel('t', fontsize=20)
    plt.ylabel('E', fontsize=20)
    plt.title('Metodo '+arregloStrings_simpl[1]+' con a: '+a_simpl, fontsize=30)
    #plt.xlim(-2.5,-0.8)
    #plt.ylim(-0.8,0.8)
    plt.scatter(data_simpl[:,0],data_simpl[:,3], s=1.5)

    datos_rk4 = sys.argv[3]
    nombreSinDat_rk4=datos_rk4.strip('.dat')
    #Outputs/int_simpletico_0.353553.dat
    arreglo_strings_grande_rk4 = nombreSinDat_rk4.split('/')
    arregloStrings_rk4 = arreglo_strings_grande_rk4[1].split('_')
    a_rk4=(arregloStrings_rk4[2].strip('.dat'))

    data_rk4 = np.loadtxt(datos_rk4)

    fig_17_rk4 = plt.figure()
    fig_17_rk4.set_size_inches(14.5,10.5)
    plt.xlabel('t', fontsize=20)
    plt.ylabel('E', fontsize=20)
    plt.title('Metodo '+arregloStrings_rk4[1]+' con a: '+a_rk4, fontsize=30)
    #plt.xlim(-2.9,-0.7)
    #plt.ylim(-1.1,1.5)
    plt.scatter(data_rk4[:,0],data_rk4[:,3], s=1.5)

    path_Imagenes="Images/"
    pp = PdfPages(path_Imagenes+"Figura_17"+'.pdf')
    pp.savefig(fig_17_simpletica)
    pp.savefig(fig_17_rk4)
    pp.close()

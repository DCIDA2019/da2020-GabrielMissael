'''
Autor: Gabriel Missael Barco, 2020
Descripcion: Este modulo tiene funciones
que realizan y analizan cadenas de Markov
con el metodo de Montecarlo, para
ajustar parametros a un modelo dado.

REQUERIMIENTOS: Numpy, pandas, matplotlib.pyplot, matplotlib, getdist, cosmolopy.distance
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from getdist import plots, MCSamples
import getdist
import pandas as pd
import cosmolopy.distance as distance
from IPython.display import clear_output

def montecarlo_mc(N, datos, teoria, desviacion_teoria, modelo, desviacion_parametros, p_old, ln_likelihood, ln_prior, n_pasos=500, seed = 500):
    '''
    REQUIREMENTS: Numpy
    DESCRIPTION: Esta funcion aplica el metodo de montecarlo para obtener los valores
                 mas probables de determinados parametros de un modelo aplicado a
                 ciertos datos conocidos.
    IN: N = numero de cadenas deseadas
        datos = variable independiente (lista con W listas de [parametros])
        teoria = variable dependiente de los parametros (lista con W numeros),
        desviacion_teoria = desviacion de las y (lista, una desviavion por punto),
        modelo = funcion que aplica el modelo (nombre),
        desviacion_parametros = desviacion estandar para cada parametro (lista),
        p_old = aproximaciónes iniciales para cada parametro (lista),
        ln_likelihood = funcion del logaritmo del likelihood (nombre),
        ln_prior = funcion del logaritmo del prior (nombre),
        n_pasos = numero de pasos por cadena (int),
        seed = semilla para numeros aleatorios (int),
    NOTA: Las funciones de likelihood y prior deben aceptar:
        ln_likelihood(datos, p_old[i_cadena], teoria, desviacion_teoria, modelo)
        ln_prior(p_old[i_cadena], teoria, desviacion_teoria, modelo)
    OUT: N cadenas de markov con n_pasos puntos, de la forma [[parametros], posterior, bool]
    '''
    #Creamos lista donde estarán las cadenas
    markov_matrix = []
    
    #Corremos N veces para obtener n cadenas
    for i_cadena in range(N):
        
        #Para tener datos reproducibles, damos un seed para los numeros aleatorios
        np.random.seed(seed)

        #Inicializamos la cadena de markov
        markov_chain = []

        #agregamos el primer conjunto de parametros y su likelihood
        ln_post_old = (ln_likelihood(datos, p_old[i_cadena], teoria, desviacion_teoria, modelo)+
                       ln_prior(p_old[i_cadena], teoria, desviacion_teoria, modelo))
        markov_chain.append([p_old[i_cadena], ln_post_old, True])

        #Realizamos n iteraciones
        for paso in range(n_pasos):
            
            print(f'Cadena {i_cadena}, {paso}/n_pasos')
            clear_output(wait=True)
            #Lista para nuevo parametro
            p_new = []

            #Llenamos p_new con una distribución aleatoria normal
            for j, parametro in enumerate(p_old[i_cadena]):
                p_new.append(desviacion_parametros[j]*np.random.randn()+parametro)

            #Calculamos nuevo likelihood
            ln_post_new = (ln_likelihood(datos, p_new, teoria, desviacion_teoria, modelo)+
                           ln_prior(p_new, teoria, desviacion_teoria, modelo))

            #Guardamos nuevo punto, su likelihood y si lo aceptamos o no segun las condiciones
            if (ln_post_new>ln_post_old):
                markov_chain.append([p_new, ln_post_new, True])
                p_old[i_cadena] = p_new
                ln_post_old = ln_post_new
            elif  ln_post_new-ln_post_old > np.random.randn(): #Este exponencial no va, pero no pude implementar el prior
                markov_chain.append([p_new, ln_post_new, True])
                p_old[i_cadena] = p_new
                ln_ost_old = ln_post_new
            else:
                markov_chain.append([p_new, ln_post_new, False])
  
        #Agregamos la cadena de esta iteración a la lista de cadenas
        markov_matrix.append(markov_chain)
            
    else:
        return markov_matrix


def aceptacion(cadenas):
    '''
    DESCRIPTION: Obtiene fracción de aceptación de las cadenas de Markov
    IN: cadenas = lista de cadenas de markov
    OUT: lista con fraccion de aceptacion de cada cadena (en orden)
    '''
    #Obtengo el numero de puntos aceptados en cada cadena
    
    N = len(cadenas)
    
    aceptacion = np.zeros(N)
    for i, cadena in enumerate(cadenas):
        for p, post, a in cadena:
            aceptacion[i] += a
            
    #Calculamos la fraccion de aceptacion
    ratio = np.zeros(N)
    pruebas = len(cadenas[0])

    for i in range(N):
        ratio[i] = aceptacion[i]/pruebas
        print(f'Aceptación de cadena {i} = {ratio[i]}')
    
    return ratio


def gelman_rubin(cadenas, cut=0):
    
    '''
    DESCRIPTION: Se realiza el diagnostico de Gelman-Rubin para todos los parametros
    IN: cadenas = cadenas de markov
        cuts = cortes del burning
    OUT: Medida de convergencia de cada parametro
    '''
    
    #Sacamos los datos para manipularlos facilmente
    cad = []
    for i, cadena in enumerate(cadenas):
        cad.append([])
        for p, post, a in cadena:
            cad[i].append(p)
            
    cad = np.array(cad)
    
    #Realizamos el diagnostico Gelman-Rubin
    
    p = len(cad[0][0]) #numero de parametros
    m = len(cad) #numero de cadenas
    n = len(cad[0][cut:]) #elementos de la cadena sin el burning
    
    R = [] #lista de medida de convergencia
    
    for i in range(p):
        s_c = []
        media = []
        
        for j in range(m):
            s_c.append(0)
            media.append(np.average((cad[j][cut:].T)[i])) #media del parametro en esa cadena
            
            for x in cad[j][cut:][i]:
                s_c[j] += (x-media[j])**2 
            s_c[j] /= (n-1)
            
        s = sum(s_c)/m
        
        Bn = 0
        miu = np.average(np.array(media))
        
        for x in media:
            Bn += (x-miu)**2
        Bn /= (m-1)
        
        sigma = (n-1)/n*s+Bn
        
        R.append((sigma/s)**(1/2))
        print(f'Para el parametro {i}, R = {R[i]:.5}')
        
    return R


def resultados(datos):
    '''
    Esta funcion regresa las medias y los intervalos de confianza para los parametros
    IN: datos = Lista de puntos resultantes de las cadenas de markov
    OUT: Lista de elementos de la forma [media i, inf_i, sup_i], para cada parametro
    '''
    N = len(datos[0])
    parametro = [[] for x in range(N)]
    
    for punto in datos:
        for i, p in enumerate(punto):
            parametro[i].append(p)
            
    R = []
    for i in range(N):
        R.append([])
        R[i].append(np.percentile(parametro[i], 50))
        R[i].append(np.percentile(parametro[i], 16))
        R[i].append(np.percentile(parametro[i], 84))
    
    for i, parametro in enumerate(R):
        print(f'Parametro {i} = {parametro[0]:.4} ± ({parametro[0]-parametro[1]:.4}, {parametro[2]-parametro[0]:.4})')
    
    return R


def ver_cadenas(cadenas, i_parametro, p_reales=[], labels = []):
    '''
    DESCRIPTION: Muestra graficamente las cadenas de markov generadas, visualizando solo dos parametros
    IN : cadenas = lista con N cadenas, cada una con M elementos de la forma [[parametros], posterior, accepted]
         i_parametros = lista de dos numeros igual al indice de los parametros que se quieren mostrar
         p_reales = lista con valores reales de todos los parametro (en orden)
         labels = lista con nombres de los parametros (strings)
    OUT: Plot con cadenas de markov sobre dos paramtros
    '''
    # Extraemos los datos
    x = []
    y = []

    for i, cadena in enumerate(cadenas):
        x.append([])
        y.append([])
        for dato, post, accept in cadena:
            x[i].append(dato[i_parametro[0]])
            y[i].append(dato[i_parametro[1]])

    
    #Creamos figura
    fig, ax0 = plt.subplots(figsize=(6, 6))

    #Graficamos todas las cadenas enviadas
    for i in range(len(cadenas)):
        ax0.scatter(x[i], y[i], marker = 'o', label = 'Cadena '+str(i+1), alpha = 0.2)

    #Si se envian los parametros enviados, se muestran en la gráfica
    if (p_reales != []):
        plt.plot(p_reales[i_parametro[0]], p_reales[i_parametro[1]], marker = 's', label = 'Valor real', c = 'black')
        
    #Detalles de la gráfica
    if labels == []:
        ax0.set_ylabel(f'Parametro {i_parametro[1]}')
        ax0.set_xlabel(f'Parametro {i_parametro[0]}')
    else:    
        ax0.set_ylabel(labels[i_parametro[1]])
        ax0.set_xlabel(labels[i_parametro[0]])
        
    ax0.grid()
    plt.legend()
    plt.title(f'MCMC aplicado al modelo')

    plt.show()


def ver_parametros(cadena, i_parametros = [], cut = 0, labels = []):
    '''
    DESCRIPTION: Sirve para observar la variacion de los parametros a lo largo de los pasos de UNA SOLA CADENA
                 y poder decidir el corte del burning para cada cadena
    IN: cadenas = lista de una cadena de markov
        i_parametros = indice de parametros a mostrar (dafault = todos)
        cut = paso a partir del que se grafica (esto es util para elegir donde cortar el burn in)
        labels = lista con nombres de los parametros (strings)
    OUT: N plots de N parametros vs pasos de la cadenas
    '''
    #Saco listas de cada parametro de la cadena (en orden)
    N = len(cadena[0][0])
    plot_parametro = [[] for i in range(N)]
    pasos = [x for x in range(len(cadena))]

    for p, post, a in cadena:
        for i, parametro in enumerate(p):
            plot_parametro[i].append(parametro)

    #Si no se especifica, se grafican todos       
    if(i_parametros == []):
        i_parametros = [x for x in range(N)]
    
    #Creo figura con una subplot para cada parametro
    fig, axs = plt.subplots(len(i_parametros), 1, sharex=True, figsize = (10, N*2))
    fig.subplots_adjust(hspace=0)
    
    #Graficamos
    axs[0].set_title(f'Parametros a cada paso (Intervalo = [{cut},{len(cadena)}])')
    for i, index in enumerate(i_parametros):
        
        media = np.percentile(np.array(plot_parametro[index][cut:]), 50)
        axs[i].plot(pasos[cut:], plot_parametro[index][cut:])
        axs[i].plot(pasos[cut:], [media for x in range(len(pasos[cut:]))], linewidth=2.0)
        
        axs[i].grid()
        if labels != []:
            axs[i].set_ylabel(labels[index])
        else:
            axs[i].set_ylabel(f'Parametro {index}')
        
    axs[N-1].set_xlabel('Pasos')
    plt.show()


def hist2d(datos, i_param, bins=100, labels = []):
    '''
    DESCRIPTION: Esta funcion genera un histograma 2d de dos parametros dados
    IN: datos = Lista con puntos de parametros
        i_param = lista indices de los dos parametros a graficar
        bins = particiones del histograma
        labels = nombre de todos los parametros
    OUT = plot histograma 2d con colorbar inferno  
    '''
    
    #Obtenemos listas de ambos parametros
    x = []
    y = []
    for punto in datos:
        x.append(punto[i_param[0]])
        y.append(punto[i_param[1]])
    
    #Graficamos
    plt.hist2d(x, y, bins = bins, cmap = 'inferno')
    
    if labels == []:
        plt.ylabel(f'Parametro {i_param[1]}')
        plt.xlabel(f'Parametro {i_param[0]}')
    else: 
        plt.ylabel(labels[i_param[1]])
        plt.xlabel(labels[i_param[0]])
        
    plt.colorbar(label = 'Frecuencia')
    plt.show()


def triangulo(datos, labels = []):
    '''
    DESCRIPTION: Esta funcion realiza una grafica triangular para visualizar los resultados
    IN: datos = lista con puntos de datos (cada uno es una lista con valores para cada parametro)
        labels = nombre de los parametros (default p_i)
    OUT: Grafica triangular con histogramas 1d y 2d de todas la combinaciones
    '''
    g = plots.get_subplot_plotter(subplot_size=2)
    
    #Checamos si se asignó nombre a los parametros
    if labels!=[]:
        samples = MCSamples(samples=np.array(datos), labels = labels, names = labels)
    
    else:
        samples = MCSamples(samples=np.array(datos))
    
    #graficamos
    g.triangle_plot(samples, filled=True, title_limit=1)


def unir_cadenas(cadenas, cuts):
    '''
    DESCRIPTION: Esta funcion une los datos de todas las cadenas, ignorando los cortados por el burning
    IN: cadenas = Lista de cadenas de markov
        cuts = lista con paso de corte por burning de cada cadena
    OUT: lista con datos de las cadenas
    '''
    #Extraemos los datos de las cadenas ignorando el burning
    datos = []
    
    for i, cadena in enumerate(cadenas):
        for p, post, a in cadena[cuts[i]:]:
            datos.append(p)
            
    return datos




#-----------------------------------------------#
#---------- Interpolación de Lagrange ----------#
#-----------------------------------------------#

# Librerías necesarias 
import numpy as np # Expresiones y cosas matemáticas
import sympy as sym # Variables simbólicas para evaluar polinomios y funciones
import matplotlib.pyplot as plt # Plots

# X: Vector que contiene las abscisas de entrada
# Y: Ordenadas para las abscisas anteriores
def interpol(X, Y):

    # Preliminares
    x = sym.Symbol('x') # Declara 'x' variable simbólica
    L = [] # Vector que recoge los coeficientes de Lagrange
    P = [] # Vector que recoge los términos de mi Polinomio de Interpolación

    # Construcción
    for i in range(len(X)): # Recorrido para cada coeficiente del polinomio
        L.append(1) # Agrego 1's para editarles luego 

        # Coeficiente de Lagrange
        for j in range(len(X)): # Recorrido por cada abscisa x_j
            if i != j : # Condicional en la formula i distinto de j
                L[i] = L[i]*(x - X[j]) / (X[i] - X[j]) # Contruye el factorial

        # Coeficiente del polinomio de interpolación 
        P.append(0) # Agrego 0’s para editarles luego 
        P[i] = L[i] * Y[i] # Formulilla

    # Output
    polinomio = sum(P).expand() 
    eval = sym.lambdify(x, polinomio)

    # Gráfica
    dom = np.linspace(start = X[0], stop = X[-1], num = 100) # Dominio de la grafica 
    plt.plot(dom, eval(dom)) + plt.plot(X, Y, "o") # plot del polinomio 
    plt.show() # Lanza la interfaz para visualizar la grafica

    print("Coeficientes de Lagrange : \n " ,
        L,
        " \nTerminos polinomio de Interpolacion : \n " ,
        P,
        " \nPolinomio simplificado : \n " ,
        polinomio
    )

interpol(
    [0 , 1 , 2 , 3 , 4 , 5 , 6],
    [1 , 3 , 2 , 1 , 3 , 2 , 1]
)


#------------------------------------------------#
#---------- Interpolación de Chebyshev ----------#
#------------------------------------------------#

# Librerías adicionales necesarias 
from cmath import cos , e , exp , pi , sin 
x = sym.Symbol('x') # Declara 'x' variable simbólica

# fun: Funcion que se desea aproximar
# n: Grado del polinomio de aproximacion
# a: Extremo inferior del intervalo
# b: Extremo superior del intervalo

def cheby(fun, n, a, b): 
    # Previa
    evalf = sym.lambdify(x, fun, 'numpy') # evalf () evalua la funcion fun
    X = [] # Vector que recoge los nodos

    for i in range(n +1): # Recorrido
        X.insert(0, (b-a) / 2*cos(( (2*i+1)*pi ) / (2*n+2) ) + (a+b)/2) # Formulilla

    # Evaluación de la función en los nodos
    Y = evalf (np.array(X))

    # Obtención de los coeficientes c_j de Chebyshev
    C = [(1/( n +1)) * sum(Y) ] # c_0 primera iteración

    for j in range(n):
        c = 0 # Variable auxiliar para recoger las sumatorias
        for k in range(n+1):
            c += (2/(n+1)) * evalf(X[k]) * cos((j+1) * np.arccos(2*((X[k]-a)/(b-a))-1) ) # Formulilla
        C.append(c) # Guarda lo calculado

    # Cálculo de los polinomios T_j(x)
    T = [1 , x ] # Vector que recoge los polinomios T_j(x) con la construccion natural
    for i in range (n -1) : # Recorrido hasta el polinomio que necesito
        T.append(2*x*T[-1] - T[-2]) # Formulilla

    # Obtención del polinomio Chebyshev
    P = [] # Vector recogedor
    for i in range ( n +1) : # Recorrido
        P.append(0) # Añado 0's para editar esa posición luego 
        if i != 0: # T (0) = 1 (No lo puedo evaluar en x)
            T[i] = T[i].subs(x, 2*((x-a) / (b-a))-1) # Cálculos 
        else:
            P[i] = C[i]

        P[i] = C[i]* T[i] # Formulilla

    polinomio = sum(P).expand() 

    # Output
    print(f" Nodos: \n{X} \nImagen en esos nodos: \n{Y} \nPolinomio Interpolador: \n{polinomio} ")

cheby(1 / (1+x**2), 20, -5, 5)
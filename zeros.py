#----------------------------------------------------#
#---------- Método de bisección de Bolzano ----------#
#----------------------------------------------------#

# Librerías necesarias
import numpy as np

# [a, b] Es el intervalo inicial con el que trabajo
# fx : Expresión matemática de mi función (requiere anteponer lambda x: ...)
# fx(x0) : Evalua el polinomio en x0
# tol : Error que tolera
# error : Dimensión del intervalo que contiene la raíz (error)
# c : La mejor aproximación que tengo a mi raíz

def bisection(a , b, fx, tol):
    
    error = abs(b-a) # Estimación del error
    S = [] # Vector donde recojo todas mis soluciones ’c ’
    contador = 0 # Contador para el número de iteraciones
    
    while tol < error: # Mientras el error no cumpla con la tolerancia
        fa = fx(a)
        fb = fx(b)
        
        if 0 < np.sign(fa) * np.sign(fb):
            print("No hay condiciones")
            break
        
        c = (a+b) /2
        S.append(c)
        fc = fx(c)
        
        if np.sign(fa) * np.sign(fc) < 0: # Análisis de signos
            b = c
        else:
            a = c
        
        error = abs(b-a)
        contador += 1
    
    print('Vector de puntos medios : ', S,
        '\nIteraciones : ', contador,
        '\nc : ' , c,
        '\nf(c) : ' , fx(c),
        '\nError : ' , error
    )

bisection(1, 2, lambda x: x**2 - 3, 0.0001) # Ejemplo


#-------------------------------------------------#
#---------- Método de la posición falsa ----------#
#-------------------------------------------------#

# Librerías necesarias
import numpy as np

# [a,b] Es el intervalo inicial con el que trabajo
# fx : Expresión matemática de mi función (requiere anteponer lambda x: ...)
# fx(x0) : Evalua el polinomio en x0
# tol : Error que tolera
# error : Dimensión del intervalo que se recorta
# c : Corte de la recta y mejor aproximación al cero

def falseposition(a, b, fx, tol):

    error = abs (b - a )
    S = [] # Vector que recoge mis soluciones
    contador = 0 # Contador de iteraciones

    while tol < error : # Mientras la tolerancia sea menor que el error
        fa = fx(a)
        fb = fx(b)

        if 0 < np.sign(fa) * np.sign(fb):
            print("No hay condiciones")
            break

        c = b - ( fb*(a-b) / (fa-fb) ) # Fórmula
        S.append(c)
        fc = fx(c)

        if np.sign(fa) * np.sign(fc) < 0: # Análisis de signos
            error = abs(b-c) # El error se mide en el intervalo que estoy recortando
            b = c
        else:
            error = abs(c-a)
            a = c

        contador += 1

    print('Vector de soluciones : ', S,
        '\nIteraciones : ', contador,
        '\nc : ', c,
        '\nf(c) : ', fx(c),
        '\nError : ', error
    )

falseposition(-1, 0, lambda x: 230*x**4 + 18*x**3 + 9*x**2 - 221*x - 9, 10**(-6)) # Ejemplo


#-------------------------------------------------#
 #---------- Método de Newton - Rhapson ----------#
#-------------------------------------------------#

# Librerías necesarias
import numpy as np

# a0 : Punto inicial
# ai : Siguiente punto
# fx : Expresión matemática de mi función (requiere anteponer lambda x: ...)
# dfx : Derivada de fx 
# fx(x0) : Evalua el polinomio en x0
# tol : Error que tolera
# error : Distancia a la raíz (error)
# c : Solución
# f(c) : Función evaluada en la solución encontrada

def newrap(a0, fx, dfx, tol):
    
    S = [a0] # Vector que recoge las soluciones
    error = 2* tol # Un error inicial para darle condiciones al ciclo while
    contador = 0 # Para contar el número de iteraciones
    
    while tol < error:
        a0 = S[-1]
        ai = a0 - (fx(a0) / dfx(a0))
        S.append(ai)
        error = abs(ai-a0)
        contador += 1
    
    c = S[-1]
    
    print('Vector de soluciones : ', S,
        '\nIteraciones : ', contador,
        '\nc : ', c,
        '\nf(c) : ' , fx(c),
        '\nError : ', error
    )

newrap(0, lambda x: 230*x**4 + 18*x**3 + 9*x**2 - 221*x - 9, lambda x : 920*x**3 + 54*x**2 + 18*x - 221, 10**(-6)) # Ejemplo


#------------------------------------------#
#---------- Metodo de la Secante ----------#
#------------------------------------------#

# Librerias necesarias
import numpy as np

# a1, a0 : Valores iniciales
# ai : Siguiente punto
# fx : Expresion matematica de mi funcion ( requiere anteponer lambda x: ...)
# fx ( x0 ) : evalua el polinomio en x0
# tol : error que tolera
# error : distancia a la raiz ( error )
# c : solucion
# f ( c ) : funcion evaluada en la solucion encontrada

def secant(a0, a1, fx, tol):

    S = [ a0 , a1 ] # Vector que recoge las soluciones
    error = 2* tol # Un error inicial para darle condiciones al ciclo while
    contador = 0 # Para contar el numero de iteraciones

    while tol < error :
        a0 = S[-2]
        a1 = S[-1]
        ai = a1 - ( (fx(a1) * (a1-a0)) / (fx(a1)-fx(a0)) )
        S.append(ai)
        error = abs(ai-a1)
        contador += 1
    
    c = S[-1]

    print('Vector de soluciones : ', S,
        '\nIteraciones : ' , contador,
        '\nc : ', c,
        '\nf(c) : ', fx(c),
        '\nError : ', error
    )

secant(-1, 0, lambda x: 230*x**4 + 18*x**3 + 9*x**2 - 221*x - 9, 10**(-6)) # Ejemplo
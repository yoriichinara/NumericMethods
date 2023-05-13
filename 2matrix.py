#---------------------------------------------#
#---------- Sustitución hacia atras ----------#
#---------------------------------------------#

import numpy as np # librería útil para cosas matemáticas

def backsus(A, B): # def de la función que recibe matriz A, vector independiente B

    A = np.array(A); B = np.array(B) # Convierto esas listas a matrices numericas

    # matriz aumentada
    AB = np.concatenate((A, B), axis=1) # .concatenate() junta vectores, axis(dimensión de la la matriz)

    # parámetros para recorrer la matriz por filas y columnas
    tam = np.shape(AB) # .shape(AB) devuelve el tamaño de la matriz AB filas x columnas
    endrow = tam[0]-1 # filas
    endcol = tam[1]-1 # columnas

    X = np.zeros(tam[0]) # vector de ceros que recoge las soluciones
    
    for i in range(endrow, -1, -1): # i para recorrer filas
        suma = 0 # variable que recoge la suma que debo restar (en la fórmula)
        
        for j in range (i+1, endcol, 1): # j para recorrer columnas
            suma = suma + AB[i, j]*X[j] # calcula la suma que debo restar en mi iteracion pues evalua los x_j que ya he encontrado en la coordenada a_ij
        
        b = AB[i, endcol] # ultimo elemento del vector independiente b en la fila i
        X[i] = (b-suma) / AB[i, i] # formulilla para la variable x_i
    
    X = np.transpose([X]) # Para ponerlo en vertical que se ve mas bonito
    
    return X


#-----------------------------------------------------------------#
#---------- Gauss con pivoteo + sustitución hacia atras ----------#
#-----------------------------------------------------------------#

def gaussprince(A, B) : # def de la funcion que recibe matriz A, vector independiente B
    
    A = np.array(A); B = np.array(B) # Convierto esas listas a matrices numericas
    
    # Matriz aumentada
    AB = np.concatenate((A, B), axis=1) # .concatenate() junta vectores, axis(dimensión de la la matriz)

    # Parámetros para recorrer la matriz por filas y columnas
    tam = np.shape(AB) # .shape(AB) devuelve tamaño matriz AB fila x columna
    n = tam [0] # tamaño fila
    m = tam [1] # tamaño columna
 
    # Eliminación hacia adelante
    for i in range(0, n-1, 1): # Recorrido por filas pivote
        
        pivote = AB[i, i] # Asignación del pivote a_ii
        adelante = i+1 # Filas que se veran afectadas
        
        for k in range(adelante, n, 1): # Rango para recorrer los que estan debajo del pivote
        
            factor = AB[k, i] / pivote # Cálculo del factor para fila k columna i
            AB[k, :] = AB[k, :] - AB[i, :]*factor # Fórmulilla
        
        newA = AB[:, :-1] # Todo menos la última columna, matriz pivoteada
        newB = []
        
        for i in range(len(newA)): # len() mide la longitud del vector
            newB.append([AB[i, -1]]) # Necesito que B este de la forma [[b1], [b2], ...] si no la otra función me bota error por la dimensión --- .append(x) añaade a x al final de la lista
        
    X = backsus(newA, newB) # vector solucion sustitucion hacia atras
    print ("X : \n", X)

# Ejemplo 
gaussprince( 
    [[1, 0, 0, 0, 0, 0, 0], 
    [1, 1, 1**2, 1**3, 1**4, 1**5, 1**6], 
    [1, 2, 2**2, 2**3, 2**4, 2**5, 2**6],
    [1, 3, 3**2, 3**3, 3**4, 3**5, 3**6],  
    [1, 4, 4**2, 4**3, 4**4, 4**5, 4**6], 
    [1, 5, 5**2, 5**3, 5**4, 5**5, 5**6], 
    [1, 6, 6**2, 6**3, 6**4, 6**5, 6**6], 
    ], [[1], [3], [2], [1], [3], [2], [1]]
)


#--------------------------------------#
#---------- Factorización LU ----------#
#--------------------------------------#

from decimal import Decimal
import numpy as np

def lumatrix(A, B): # B es el vector de valores in dependie ntes

    # Previa
    B = np.transpose([B]) # Acomoda el B como vector vertical
    AB = np.concatenate ((A, B), axis=1) # Matriz aumentada

    # Pivoteo parcial por filas
    tamano = np.shape(AB) # Tamaño matriz aumentada nxm
    n = tamano[0] # num filas
    m = tamano[1] # num columnas

    # Eliminación gaussiana hacia adelante
    L = np.identity(n, dtype = float) # Inicia L como matriz identidad nxn

    for i in range(0, n-1, 1): # Recorrido por filas
        
        pivote = AB[i, i] # Cálcula pivote
        adelante = i+1 # Solo filas debajo del pivote
        
        for k in range(adelante, n, 1): # Cada elemento debajo del pivote
            
            factor = AB[k, i] / pivote # Factor
            AB[k, :] = AB[k, :] - AB[i, :]*factor # Fórmulilla
            L[k, i] = factor # Recoge los factores
        
    U = np.copy(AB[:, :m-1]) # Matriz U, toda la anterior menos el vector de val independientes 
    LB = np.concatenate((L, B), axis =1) # Matrizaumentada
    
    # sustitucion hacia adelante
    Y = np.zeros(n, dtype = float) # vector ceros que recoge las soluciones
    Y[0] = LB[0, n] # primer valor independiente que arranca el despeje

    for i in range (1, n, 1): # recorrido por filas
        suma = 0 # inicia lizacion suma
    
        for j in range (0 ,i ,1) : # rango columnas ( variables a reemplazar para luego despejar    
            suma = suma + LB[i, j]*Y[j] # calculo suma para el despeje
        
        b = LB[i, n] # valor independiente del vector B
        Y[i] = (b-suma) / LB[i, i] # formulilla
    
    Y = np.transpose([Y]) # dispone Y como un vector vertical

    # Resolver UX = Y
    UY = np.concatenate((U, Y), axis=1) # matrix aumentada
    
    # sustitucion hacia atras
    ultfila = n-1 # ultima fila
    ultcolumna = m-1 # ultima columna
    X = np.zeros(n, dtype = float) # vector ceros que recoge soluciones
    
    for i in range (ultfila, 0-1, -1): # recorrido filas
        suma = 0 # inicia variable
        
        for j in range ( i +1 , ultcolumna ,1) : # recorrido columnas
            suma = suma + UY[i, j]*X[j] # calcula la suma evaluando x_j encontrados
        
        b = UY[i, ultcolumna] # valor independiente vector B
        X[i] = (b-suma) / UY[i, i] # formulilla
    
    X = np.transpose([X]) # dispone X como vector
    
    print("solucion : \n ", X)
    
    return X

# Ejemplo
A = np.array([
        [1, 0, 0, 0, 0, 0, 0], # np.array() para verla y trabajarla como matrix
        [1, 1, 1, 1, 1, 1, 1],
        [1, 2, 4, 8, 16, 32, 64],
        [1, 3, 9, 27, 81, 243, 729],
        [1, 4, 16, 64, 256, 1024, 4096],
        [1, 5, 25, 125, 625, 3125, 15625],
        [1, 6, 36, 216, 1296, 7776, 46656]], dtype = float) 

b1 = lumatrix(A, np.array([1, 0, 0, 0, 0, 0, 0], dtype = float))
b2 = lumatrix(A, np.array([0, 1, 0, 0, 0, 0, 0], dtype = float))
b3 = lumatrix(A, np.array([0, 0, 1, 0, 0, 0, 0], dtype = float))
b4 = lumatrix(A, np.array([0, 0, 0, 1, 0, 0, 0], dtype = float))
b5 = lumatrix(A, np.array([0, 0, 0, 0, 1, 0, 0], dtype = float))
b6 = lumatrix(A, np.array([0, 0, 0, 0, 0, 1, 0], dtype = float))
b7 = lumatrix(A, np.array([0, 0, 0, 0, 0, 0, 1], dtype = float))

INV = np.concatenate((b1, b2, b3, b4, b5, b6, b7), axis=1)
Id = np.round(np.dot(A,INV))

print("INV = Array de los b_i :\n", INV, "\n A*INV redondeado : \n", Id)
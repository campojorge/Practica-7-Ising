# Practica-7-Ising


Tests para verificar que N-Fold Way reproduce la dinámica de Glauber:

1. Misma semilla con el mismo algoritmo debe dar los mismos resultados (se generan los mismos números aleatorios en el mismo orden)

Si para la misma semilla (fija) obtenemos distintos valores en el fichero es porque N-fold Way no es reproducible y por tanto no reproduce Glauber.

2. Mismas condiciones iniciales para los dos algoritmos deben dar resultados estadísticamente compatibles.

Se puede tomar un histograma del estado final (la  magnetización por ejemplo), ajustar a una función determinada y comporobar que los resultados son estadísticamente compatibles.(ya implementado el histograma de la distribucion de tiempos)

3. Comprobar que Q10 equivale a la suma de las probabilidades de Glauber.
4. Comparar el número de intentos hasta el primer flip en los dos algoritmos

Generamos una configuración (siempre la misma) un número grande de veces (p.ej 10^5), comprobamos que el número de pasos de montecarlo que ocurren hasta el primer flip en glauber se corresponde estadísticamente con el primer valor de m en N-fold Way.

5. Comparar las curvas <M(t),E(t)>; con t en MCSS (test muy sensible a la dinámica)
6. Definir un funcion correlacion C(t)=<1/N*sumatorio(sigma(0)sigma(t))>

Se trata de un funcion que da la correlacion entre la configuracion inicial y la que hay en un instante t. Cosas a comprobar: 1-Tiene que ser una funcion que vaya del 0 al 1. 2-Para H=J=0 la sabemos la curva teorica por lo que podemos comparar. (Ya implementado). 3- Las curvas deben de ser iguales para ambos metodos dentro de las incertidumbres.
Problema:Como hemos hecho las simulaciones hasta M=0, hay simulaciones que lo alcanzan antes y otras despues. Por ello, para hacer el promedio temporal tengo quedarme con el t_min que no es el mismo para glauber que par n-fold y por lo tanto no se observa bien las cosas. Hay que cambiar esto.

7. Se pueden definir otras funciones de correlacion(Como la correlacion espacial)
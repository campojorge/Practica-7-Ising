# Practica-7-Ising


Tests para verificar que N-Fold Way reproduce la dinámica de Glauber:

1. Misma semilla con el mismo algoritmo debe dar los mismos resultados (se generan los mismos números aleatorios en el mismo orden)

Si para la misma semilla (fija) obtenemos distintos valores en el fichero es porque N-fold Way no es reproducible y por tanto no reproduce Glauber.

2. Mismas condiciones iniciales para los dos algoritmos deben dar resultados estadísticamente compatibles.

Se puede tomar un histograma del estado final (la  magnetización por ejemplo), ajustar a una función determinada y comporobar que los resultados son estadísticamente compatibles

3. Comprobar que Q10 equivale a la suma de las probabilidades de Glauber.
4. Comparar el número de intentos hasta el primer flip en los dos algoritmos

Generamos una configuración (siempre la misma) un número grande de veces (p.ej 10^5), comprobamos que el número de pasos de montecarlo que ocurren hasta el primer flip en glauber se corresponde estadísticamente con el primer valor de m en N-fold Way.

5. Comparar las curvas <M(t),E(t)>; con t en MCSS (test muy sensible a la dinámica)
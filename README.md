---
title: "Modelos de Vectores Autoregresivos VAR"
output: html_document
---

<br><br><br>


<p style="text-align:center;">
*Al igual que en la entrada de Series de Tiempo, este documento es únicamente una guía que me sirve a mí para practicar. No pretendo que funcione como un curso pedagógico para nadie. Si decides que algo te sirve, adelante. Sin embargo, por favor recuerda ir a clases, porque seguramente me equivoqué en algo.*
</p>

<p style="text-align:center;">
-C
</p>

<br><br><br><br><br><br>


```{r}
#Vamos a necesitar cargar las siguientes librerías:
library(vars)
library(tseries)
library(forecast)
library(urca)
library(highcharter)
```

<br><br><br>

A diferencia de los modelos $AR(p)$ donde las variables endógenas que determinaban el comportamiento futuro de la serie eran retrasos de la serie misma. En los modelos de vectores autoregresivos ($VAR$) el comportamiento futuro de la serie  será determinado por un vector de variables endógenas. Es decir, el vector de variables endógenas no depende únicamente de los valores pasados de sí mismas, sino de los valores rezagados de otras variables, permitiendo así estimar un sistema de ecuaciones, a diferencia de los modelos ARIMA tradicionales que permiten estimar sólo una ecuación.
 

Como ejemplo, supón que tenemos tres variables de series de tiempo diferentes, denotadas por $x_{t,1},x_{t,2},x_{t,3}$. El vector autorregresivo de orden 1 $VAR(1)$ de la forma reducida de este modelo está dado por
\begin{align}
      x_{t,1}&=\alpha_{1}+\beta_{11}x_{t-1,1}+\beta_{12}x_{t-1,2}+\beta_{13}x_{t-1,3}+\varepsilon_{t,1}\\ \\
      x_{t,2}&=\alpha_{2}+\beta_{21}x_{t-1,1}+\beta_{22}x_{t-1,2}+\beta_{23}x_{t-1,3}+\varepsilon_{t,2} \tag{1.1} \\ \\
      x_{t,3}&=\alpha_{3}+\beta_{31}x_{t-1,1}+\beta_{32}x_{t-1,2}+\beta_{33}x_{t-1,3}+\varepsilon_{t,3}
\end{align}

Cada variable es una función lineal de los rezagos de un periodo de sí mismas y de las otras dos variables.

Este modelo se puede generalizar a $n$ variables y expresar matricialmente de la siguiente manera:
$$\underbrace{\begin{pmatrix}x_{t,1} \\ x_{t,2} \\  \vdots \\ x_{t,n}\end{pmatrix}}_{X_{t}} = \underbrace{\begin{pmatrix} \alpha_{1} \\ \alpha_{2} \\ \vdots \\ \alpha_{n}\end{pmatrix}}_A + \underbrace{\begin{pmatrix} \beta_{11} & \beta_{12} & \cdots & \beta_{1n}\\ \beta_{21} & \beta_{22} & \cdots & \beta_{2n} \\ \vdots & \vdots & \ddots & \vdots\\ \beta_{n1} & \beta_{n2} & \cdots & \beta_{nn}\end{pmatrix}}_{B} \cdot \underbrace{\begin{pmatrix} x_{t-1,1} \\ x_{t-1,2} \\  \vdots \\ x_{t-1,n} \end{pmatrix}}_{X_{t-1}} + \underbrace{\begin{pmatrix} \varepsilon_{t,1} \\ \varepsilon_{t,2} \\ \vdots \\ \varepsilon_{t,n} \end{pmatrix}}_{E_{t}} \tag{1.2}$$
Luego, podemos resumir el modelo general de la ecuación como $X_{t}=A+BX_{t-1}+E_{t}$ es decir, un sistema de ecuaciones en diferencias de orden 1 con un termino estocastico ($E_{t}$).


En general, podemos determinar un $AR(p)$ rezagado $p$ con $n$ variables periodos como
$$X_{t}=A+\sum_{i=1}^{p}B_{i}X_{t-i}+E_{t}$$
Donde cada $B_{i}$ es la matriz de $n\times n$ coeficientes asociadas al periodo $t-i$.


## Construcción del Modelo

Debemos seguir el procedimiento iterativo de Box-Jenkins:

1. $\textbf{Identificación}:$ En lugar de identificar un solo modelo ARIMA, debes identificar las relaciones entre las series de tiempo y determinar cuántos rezagos (lags) son necesarios para modelar adecuadamente las dependencias temporales. Esto implica analizar las funciones de autocorrelación y autocorrelación parcial de cada serie de tiempo, así como las funciones de autocorrelación cruzada entre las series. 

El primer paso en el procedimiento iterativo de Box-Jenkins para el análisis de series de tiempo implica la identificación de un modelo apropiado. Esto incluye la determinación de los rezagos (lags) necesarios y la comprensión de la estructura de la serie de tiempo. 


#### Cargar y preparar los datos:
```{r}
data <- as.data.frame(Canada)                        #Esta es una base de datos de la librería vars
colnames(data) <- c("Employment", "Productivity",
                    "Real Wage", "Unemployment Rate")
data <- ts(data, start = c(1980,1),frequency = 4)

#Gráfico de las Series
par(mfrow = c(2, 2))
ts.plot(data[,1], 
        main = colnames(data)[1],
        ylab = "Valor")
ts.plot(data[,2], 
        main = colnames(data)[2],
        ylab = "Valor")
ts.plot(data[,3], 
        main = colnames(data)[3],
        ylab = "Valor")
ts.plot(data[,4], 
        main = colnames(data)[4],
        ylab = "Valor")
```

#### Pruebas de estacionariedad:

Al igual que en los modelos univariados para utilizar la metodología $VAR$ para el modelamiento y posterior pronóstico de series de tiempo se debe verificar la estacionaridad de la serie. Esta prueba se puede hacer a través de un test de Dickey-Fuller (adf.test()) o con la prueba de Phillips-Perron (pp.test()) serie a serie.

```{r Empleo}
adf.test(data[,1])
pp.test(data[,1])
```
Ambas pruebas arrojan un $p-value$ muy elevado, por lo tanto podemos decir que la serie de empleo es estacionaria, por lo que habrá que diferenciar la serie en el siguiente paso.


```{r Productividad}
adf.test(data[,2])
pp.test(data[,2])
```


```{r Salario Real}
adf.test(data[,3])
pp.test(data[,3])
```

```{r Tasa de Desempleo}
adf.test(data[,4])
pp.test(data[,4])
```

Ya que ninguna serie paso el test de estacionalidad tenemos que diferenciar las series mediante la funcion diff()
```{r Primera diferencia}
data.diff <- diff(data[,c(1:4)])

#Empleo
adf.test(data.diff[,1])
pp.test(data.diff[,1])

#Productividad
adf.test(data.diff[,2])
pp.test(data.diff[,2])

#Salario Real
adf.test(data.diff[,3])
pp.test(data.diff[,3])

#Tasa de Desempleo
adf.test(data.diff[,4])
pp.test(data.diff[,4])

```

#### Gráfica de las series diferenciadas:

```{r}
par(mfrow = c(2, 2))
ts.plot(data.diff[,1], 
        main = colnames(data.diff)[1],
        ylab = "Valor")
ts.plot(data.diff[,2], 
        main = colnames(data.diff)[2],
        ylab = "Valor")
ts.plot(data.diff[,3], 
        main = colnames(data.diff)[3],
        ylab = "Valor")
ts.plot(data.diff[,4], 
        main = colnames(data.diff)[4],
        ylab = "Valor")
```


#### Estimación:

Una vez que tenemos las series estacionarias procedemos a hacer las estimaciones correspondientes.


```{r}
## Debemos conocer el numero de rezagos apropiados para nuestro modelo:
VARselect(data.diff, lag.max = 8,
          type = "const")         #Esta función nos va a permitir conocer el numero de rezagos apropiados para simular nuesto modelo.


## Generar el modelo mediante la función VAR():
modelo_var<-VAR(data, p=2)        #El valor p=1 es el número de rezagos lo que modelará un VAR(1)
summary(modelo_var)
```
La elección del $VAR(1)$ fue debido a que los valores $p$ asociados a las regresiones eran estadisticamente significativos cuando utilizabamos rezagos más elevados.
 

Gráficamos el modelo:
```{r}
plot(modelo_var)
```
```{r}
roots(modelo_var)
```
Al ser todas menores a 1 este modelo cumple el supuesto de estacionieradad.

#### Comprobación de la autocorrelación en los errores:


La comprobación de la autocorrelación en los errores es una parte importante del análisis de modelos de series temporales y regresión. La autocorrelación se refiere a la correlación de un valor en una serie temporal o una secuencia de datos con sus valores anteriores o posteriores en la misma serie. En el contexto de los modelos de regresión, la autocorrelación en los errores puede indicar que el modelo no está capturando adecuadamente la estructura temporal de los datos o que los errores no son independientes y, por lo tanto, no se cumplen los supuestos básicos del modelo


```{r}
serial.test(modelo_var,
            lags.pt = 16, 
            type = "PT.adjusted")
```

La función serial.test() se utiliza para realizar prueba de Ljung-Box para evaluar la autocorrelación en los residuos del modelo. Estas pruebas evalúan si hay correlación en los residuos en diferentes retrasos (lags).

El resultado de la prueba será un valor p (p-value) que indica la probabilidad de observar la autocorrelación en los residuos si los residuos fueran realmente independientes y no autocorrelacionados. En general, si el valor p es pequeño (por ejemplo, menor que 0.05), se considera evidencia de autocorrelación significativa en los residuos, lo que sugiere que el modelo puede no ser adecuado.

Por lo tanto, en este ejemplo, no hay autocorrelacion de los errores

#### Prueba de nomalidad:

```{r}
normality.test(modelo_var,
               multivariate.only = F)
```
* *Si p.value es menor que el nivel de significancia (0.05), se rechaza la hipótesis nula y se concluye que los residuos no siguen una distribución normal.*
* *Si p.value es mayor o igual al nivel de significancia, no hay suficiente evidencia para afirmar que los residuos no siguen una distribución normal.*

```{r}

residuos <- residuals(modelo_var)

par(mfrow = c(2, 1))

hist(residuos, main = "Histograma de Residuos", xlab = "Valor Residual", prob = TRUE, col = "lightblue")


mu <- mean(residuos)
sigma <- sd(residuos)
x <- seq(min(residuos), max(residuos), length = 100)
y <- dnorm(x, mean = mu, sd = sigma)
lines(x, y, col = "red", lwd = 2)


legend("topright", legend = "Densidad Normal", col = "red", lty = 1, lwd = 2)

# Crear un gráfico QQ para evaluar la normalidad
qqnorm(residuos)
qqline(residuos)
```

#### Pronostico:

```{r}
pronosticos <- predict(modelo_var,
                       n.ahead = 2)                    #El pronostico del modelo_var se realizará para                                                        #dos trimestres adelante para cada una de las series.

## Pronosticos para las series:

pronostico_e <- pronosticos$fcst$Employment            #Pronostico de la serie de empleo

pronostico_p <- pronosticos$fcst$Productivity          #Pronostico de la serie de productividad

pronostico_rw <- pronosticos$fcst$Real.Wage.           #Pronostico de la serie de salario real

pronostico_ur <- pronosticos$fcst$Unemployment.Rate.   #Pronostico de la serie de tasa de desempleo

## Intervalos de condianza de cada pronostico:

# Empleo
intervalos_inf_e <- pronosticos$fcst$Employment[,2]
intervalos_sup_e <- pronosticos$fcst$Employment[,3]

# Productividad
intervalos_inf_p <- pronosticos$fcst$Productivity[,2]
intervalos_sup_p <- pronosticos$fcst$Productivity[,3]

# Salario Real
intervalos_inf_rw <- pronosticos$fcst$Real.Wage[,2]
intervalos_sup_rw <- pronosticos$fcst$Real.Wage[,3]

# Tasa de desempleo
intervalos_inf_ur <- pronosticos$fcst$Unemployment.Rate[,2]
intervalos_sup_ur <- pronosticos$fcst$Unemployment.Rate[,3]

```






<br><br><br><br><br><br>



<p style="text-align:center;">
*Cuiudate bb*
</p>





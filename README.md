
# Replication of of Table 17.3 of Chatper 17 of S&W 2020

## R code

Run the following code chunks using *R*.

### Part 1

Estimation results for AR and ADL models.

```

source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/r_scripts/r_script_part_01.R", print.eval=TRUE)

```

### Part 2

Pseudo out-of-sample forecasting comparison and replication of Table 17.3.

```

source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/r_scripts/r_script_part_02.R", print.eval=TRUE)

```

### Part 3

Extension of Table 17.3 with unconditional mean forecast.

```

source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/r_scripts/r_script_part_03.R", print.eval=TRUE)

```

### Part 4

Extension of Table 17.3 with POSS for other AR and ADL model forecasts.

```

source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/r_scripts/r_script_part_04.R", print.eval=TRUE)

```

### Part 5

Extension of Part 2 to Part 4 with the estimation of the unobserved factors and DFM parameters.

Note, there are different strategies to estimate the unobserved factors and DFM parameters. In the textbook S&W propose to estimate the unobserved factors by the principal components of the observed variables. This strategy is a nonparametric method and is based on cross-sectional averaging. (see Stock and Watson (2016)). Parameteric methods are state-space methods.

```

source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/ch-17-06/dfm_script_01.R", print.eval=TRUE)

```

<!--  
* Stock, J. H., & Watson, M. W. (2016). Dynamic factor models, factor-augmented vector autoregressions, and structural vector autoregressions in macroeconomics. In *Handbook of macroeconomics* (Vol. 2, pp. 415-525). Elsevier. *Ch.: 2.3*.
* Kilian, L., & Lütkepohl, H. (2017). *Structural vector autoregressive analysis*. Cambridge University Press. *Ch.: 16*.
-->

# Next Steps

* Use actual PCA to construct/estimate unobserved factors.
* Replication of the construction of the unbobserved, i.e., estimated factors.
# jaws

Estimate sparse loadings (i.e., coefficients) of Principal Component Analysis, Logistic Factor Analysis, and other techniques in the context of Latent Variable Models. Generally, this can facilitate calculation of shrunken R^2 and related quantities that represent estimated latent variables more accurately. Using systematic variation driven by latent variables, this package also estimate covariance matrices of high-dimensional data when a number of rows (variables) is exceedingly larger than a number of observations (columns).

# Installation

```R
install.packages("devtools")
library("devtools")
install_github("ncchung/jaws")
```

# Gin

Gin is a collection of pure R tools for generating and manipulating Gaussian 
Process models (GPs). It is also a nice [spirit](https://en.wikipedia.org/wiki/Gin)
![gin](/figs/martini.png)

Gin is a blend of the words *G*aussian process *IN*ference.

## Outline

The current version includes the top-level functions:



## Installation

Gin is an R package, but is still in development. To set up from GitHub first 
install (if you haven't already) Hadley Wickham's devtools package.
```
   install.packages("devtools")
```
Now you can install gin straight from GitHub:
```
   devtools::install_github("svdataman/gin")
```
Now you're good to go.

## References

For more on GPs, the best reference is:

* C. E. Rasmussen & K. I. Williams, _Gaussian Processes for Machine Learning_, 
[online here](http://www.gaussianprocess.org/gpml/chapters/)
# mcglm 0.9.0

The `mcglm` package fits multivariate covariance generalized linear
models (Bonat and Jorgensen, 2016).

## Introduction

`mcglm` is an R package designed to fit Multivariate Covariance
Generalized Linear Models. It allows you to specify a distinct linear
predictor for each response variable, offering exceptional flexibility
for analyses involving multiple outcomes.

With `mcglm`, you can model a wide range of response types — continuous,
discrete (such as counts and binary), limited, and even zero inflated
responses, whether continuous or mixed.

Its main strength lies in the ability to capture complex relationships
between variables through multiple covariance structures, enabling more
realistic and robust multivariate modeling.

This package was developed as part of the Wagner’s Ph.D. thesis,
combining academic rigor with practical value for the statistical
modeling community.

## Download and install

### Linux/Mac

Use the `devtools` package (available from
[CRAN](http://cran-r.c3sl.ufpr.br/web/packages/devtools/index.html)) to
install automatically from this GitHub repository:

    library(devtools)
    install_github("bonatwagner/mcglm")

## Authors

-   [Wagner Hugo Bonat](https://www.linkedin.com/in/wagner-bonat)
    (author and main developer)

## Contributing

This R package is develop using
[`roxygen2`](https://github.com/r-lib/roxygen2) for documentation and
[`devtools`](https://github.com/r-lib/devtools) to check and build.
Also, we adopt the [Gitflow
worflow](https://nvie.com/posts/a-successful-git-branching-model/) in
this repository.

## Instructions for contributing

Please, see the [instructions for
contributing](https://github.com/bonatwagner/mcglm/blob/main/CONTRIBUTING.md)
to collaborate.

<!-- links -->

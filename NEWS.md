# LorenzRegression 2.0.0

**Breaking Changes**:

* The package structure has been standardized and builds upon existing
    libraries. Please review the updated function documentation for
    current usage.
* `Lorenz.Reg`: The function structure has changed. It now acts as a
    wrapper for the fitting functions `Lorenz.GA`, `Lorenz.FABS`and
    `Lorenz.SCADFABS`. It returns an object of class `''LR''`
    (non-penalized regression) or `''PLR''` (penalized regression). Each
    class has a set of designated methods.
* `Lorenz.boot`: This function performs bootstrap calculations for
    objects of class `''LR''` or `''PLR''`. It uses the `boot` function
    from the `boot` package. For penalized regression, it also computes
    an out-of-bag score for tuning parameter selection. The function
    returns the updated object with bootstrap results and adds the class
    `''LR_boot''` (non-penalized regression) or `''PLR_boot''`
    (penalized regression).
* `PLR.CV`: This function performs cross-validation for objects of
    class `''PLR''`. It computes a cross-validation score for tuning
    parameter selection. The folds are constructed using `vfold_cv` from
    the `rsample` package. The function returns the updated object with
    cross-validation results and adds the class `''PLR_cv''`.
* Method availability for classes `''LR''` and `''PLR''` is now
    documented in the `Lorenz.Reg` help page. Each method also has its
    own help page.

**New Features**:

* The `grid.arg` and `grid.value` arguments in `Lorenz.Reg` allow
    users to specify one tuning parameter for penalized regression and
    construct a grid for it. Fitting is repeated for each grid value,
    and optimal values are determined using available methods (among
    BIC, bootstrap and cross-validation).
* `diagnostic.PLR` provides diagnostic information for penalized
    Lorenz regression.
* Methods `fitted`, `explainedIneq`and `autoplot` were added for objects of classes `''LR''` and `''PLR''`.

**Technical changes**

* Modifications occurred in `Lorenz.GA` (and the underlying Rcpp
    functions) in order to ensure full reproducibility of results.
* Functions implying randomness have a `seed` argument. This argument
    sets a local seed, used for the generation of random objects within
    the function. the seed is reverted to its previous state after the
    operation. This ensures that the seed settings do not interfere with
    the global random state or other parts of the code.

# LorenzRegression 1.0.0

* Added a `NEWS.md` file to track changes to the package.

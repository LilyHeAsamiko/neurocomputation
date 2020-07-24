# first line: 286
def path_scores(solver, X, y, mask, alphas, l1_ratios, train, test,
                solver_params, is_classif=False, n_alphas=10, eps=1E-3,
                key=None, debias=False, Xmean=None,
                screening_percentile=20., verbose=1):
    """Function to compute scores of different alphas in regression and
    classification used by CV objects

    Parameters
    ----------
    X : 2D array of shape (n_samples, n_features)
        Design matrix, one row per sample point.

    y : 1D array of length n_samples
        Response vector; one value per sample.

    mask : 3D arrays of boolean
        Mask defining brain regions that we work on.

    alphas : list of floats
        List of regularization parameters being considered.

    train : array or list of integers
        List of indices for the train samples.

    test : array or list of integers
        List of indices for the test samples.

    l1_ratio : float in the interval [0, 1]; optional (default .5)
        Constant that mixes L1 and TV (resp. Graph-Net) penalization.
        l1_ratio == 0: just smooth. l1_ratio == 1: just lasso.

    eps : float, optional (default 1e-3)
        Length of the path. For example, ``eps=1e-3`` means that
        ``alpha_min / alpha_max = 1e-3``.

    n_alphas : int, optional (default 10).
        Generate this number of alphas per regularization path.
        This parameter is mutually exclusive with the `alphas` parameter.

    solver : function handle
       See for example tv.TVl1Classifier documentation.

    solver_params: dict
       Dictionary of param-value pairs to be passed to solver.
    """
    if l1_ratios is None:
        raise ValueError("l1_ratios must be specified!")

    # misc
    _, n_features = X.shape
    verbose = int(verbose if verbose is not None else 0)

    # Univariate feature screening. Note that if we have only as few as 100
    # features in the mask's support, then we should use all of them to
    # learn the model i.e disable this screening)
    do_screening = (n_features > 100) and screening_percentile < 100.
    if do_screening:
        X, mask, support = _univariate_feature_screening(
            X, y, mask, is_classif, screening_percentile)

    # crop the mask to have a tighter bounding box
    mask = _crop_mask(mask)

    # get train and test data
    X_train, y_train = X[train].copy(), y[train].copy()
    X_test, y_test = X[test].copy(), y[test].copy()

    # it is essential to center the data in regression
    X_train, y_train, _, y_train_mean, _ = center_data(
        X_train, y_train, fit_intercept=True, normalize=False,
        copy=False)

    # misc
    if isinstance(l1_ratios, numbers.Number):
        l1_ratios = [l1_ratios]
    l1_ratios = sorted(l1_ratios)[::-1]  # from large to small l1_ratios
    best_score = -np.inf
    best_secondary_score = -np.inf
    best_l1_ratio = l1_ratios[0]
    best_alpha = None
    best_init = None
    all_test_scores = []
    if len(test) > 0.:
        # do l1_ratio path
        for l1_ratio in l1_ratios:
            this_test_scores = []

            # make alpha grid
            if alphas is None:
                alphas_ = _space_net_alpha_grid(
                    X_train, y_train, l1_ratio=l1_ratio, eps=eps,
                    n_alphas=n_alphas, logistic=is_classif)
            else:
                alphas_ = alphas
            alphas_ = sorted(alphas_)[::-1]  # from large to small l1_ratios

            # do alpha path
            if best_alpha is None:
                best_alpha = alphas_[0]
            init = None
            for alpha in alphas_:
                # setup callback mechanism for early stopping
                early_stopper = _EarlyStoppingCallback(
                    X_test, y_test, is_classif=is_classif, debias=debias,
                    verbose=verbose)
                w, _, init = solver(
                    X_train, y_train, alpha, l1_ratio, mask=mask, init=init,
                    callback=early_stopper, verbose=max(verbose - 1, 0.),
                    **solver_params)

                # We use 2 scores for model selection: the second one is to
                # disambiguate between regions of equivalent Spearman
                # correlations
                score, secondary_score = early_stopper.test_score(w)
                this_test_scores.append(score)
                if (np.isfinite(score) and
                        (score > best_score
                         or (score == best_score and
                             secondary_score > best_secondary_score))):
                    best_secondary_score = secondary_score
                    best_score = score
                    best_l1_ratio = l1_ratio
                    best_alpha = alpha
                    best_init = init.copy()
            all_test_scores.append(this_test_scores)
    else:
        if alphas is None:
            alphas_ = _space_net_alpha_grid(
                X_train, y_train, l1_ratio=best_l1_ratio, eps=eps,
                n_alphas=n_alphas, logistic=is_classif)
        else:
            alphas_ = alphas
        best_alpha = alphas_[0]

    # re-fit best model to high precision (i.e without early stopping, etc.)
    best_w, _, init = solver(X_train, y_train, best_alpha, best_l1_ratio,
                             mask=mask, init=best_init,
                             verbose=max(verbose - 1, 0), **solver_params)
    if debias:
        best_w = _EarlyStoppingCallback(
            X_test, y_test, is_classif=is_classif, debias=debias,
            verbose=verbose)._debias(best_w)

    if len(test) == 0.:
        all_test_scores.append(np.nan)

    # unmask univariate screening
    if do_screening:
        w_ = np.zeros(len(support))
        if is_classif:
            w_ = np.append(w_, best_w[-1])
            w_[:-1][support] = best_w[:-1]
        else:
            w_[support] = best_w
        best_w = w_

    if len(best_w) == n_features:
        if Xmean is None:
            Xmean = np.zeros(n_features)
        best_w = np.append(best_w, 0.)

    all_test_scores = np.array(all_test_scores)
    return (all_test_scores, best_w, best_alpha, best_l1_ratio, alphas_,
            y_train_mean, key)

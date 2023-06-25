test_that("Test OptimalAdaptiveIFRP and AdaptiveIFRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)
  n_observations = nrow(returns)

  covariance_factors_returns = stats::cov(factors, returns)
  variance_returns = stats::cov(returns)
  mean_returns = colMeans(returns)
  mean_factors = colMeans(factors)

  penalty_parameters = c(
    0.,
    exp(seq(from=log(1.0e-8), to=log(1e-2), length.out=100)),
    seq(1e-2, 1e2, 100)
  )

  ifrp = IFRP(returns, factors, include_standard_errors = TRUE)

  expect_no_error(
    OptimalAdaptiveIFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE
  ))

  expect_error(OptimalAdaptiveIFRP(
    t(returns),
    factors,
    penalty_parameters
  ))
  expect_error(OptimalAdaptiveIFRP(
    returns,
    t(factors),
    penalty_parameters
  ))
  expect_error(OptimalAdaptiveIFRP(
    t(returns),
    t(factors),
    penalty_parameters
  ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters = c("r", "s")
  ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      weighting_type = 4
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      tuning_type = 4
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      include_standard_errors = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      one_stddev_rule = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      gcv_scaling_n_assets = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      gcv_identifiation_check = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      n_folds = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      n_train_observations = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      n_test_observations = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      roll_shift = "r"
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      t(returns),
      t(factors),
      penalty_parameters,
      check_arguments = "r"
    ))

  weights =  intrinsicFRP:::AdaptiveWeightsCpp(returns, factors, 'c')

  aifrp0 = AdaptiveIFRP(
    returns,
    factors,
    penalty_parameters = 0,
    weights
  )

  expect_equal(aifrp0, ifrp$risk_premia)

  aifrp1 = AdaptiveIFRP(
    returns,
    factors,
    penalty_parameters = .1,
    weights
  )

  aifrpm =  AdaptiveIFRP(
    returns,
    factors,
    penalty_parameters,
    weights
  )

  expect_length(aifrp0, n_factors)
  expect_length(aifrp1, n_factors)
  expect_equal(ncol(aifrpm), length(penalty_parameters))

  ### GCV
  for (weighting_type in c('c', 'b', 'a', 'n')) {
    for (one_stddev_rule in c(TRUE, FALSE)) {
      for (gcv_vr_weighting in c(TRUE, FALSE)) {
        for (gcv_scaling_n_assets in c(TRUE, FALSE)) {
          for (gcv_identification_check in c(TRUE, FALSE)) {

            adaptive_ifrp = OptimalAdaptiveIFRP(
              returns,
              factors,
              penalty_parameters,
              weighting_type = weighting_type,
              tuning_type = 'g',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              gcv_vr_weighting = gcv_vr_weighting,
              gcv_scaling_n_assets = gcv_scaling_n_assets,
              gcv_identification_check = gcv_identification_check
            )

            expect_length(adaptive_ifrp$risk_premia, n_factors)
            expect_length(adaptive_ifrp$standard_errors, n_factors)
            expect_length(adaptive_ifrp$penalty_parameter, 1)


            adaptive_ifrp0 = OptimalAdaptiveIFRP(
              returns,
              factors,
              penalty_parameters = 0.,
              weighting_type = weighting_type,
              tuning_type = 'g',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              gcv_vr_weighting = gcv_vr_weighting,
              gcv_scaling_n_assets = gcv_scaling_n_assets,
              gcv_identification_check = gcv_identification_check
            )

            expect_equal(
              matrix(adaptive_ifrp0$risk_premia, n_factors, 1),
              ifrp$risk_premia
            )

            if (weighting_type == 'c') {

              adaptive_ifrp1 = OptimalAdaptiveIFRP(
                returns,
                factors,
                penalty_parameters = .1,
                weighting_type = weighting_type,
                tuning_type = 'g',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                gcv_vr_weighting = gcv_vr_weighting,
                gcv_scaling_n_assets = gcv_scaling_n_assets,
                gcv_identification_check = gcv_identification_check
              )

              expect_equal(
                matrix(adaptive_ifrp1$risk_premia, n_factors, 1),
                aifrp1
              )

            }
          }
        }
      }
    }
  }

  ### CV
  for (weighting_type in c('c', 'b', 'a', 'n')) {
    for (one_stddev_rule in c(TRUE, FALSE)) {
      for (n_folds in c(5, 10)) {

        adaptive_ifrp = OptimalAdaptiveIFRP(
          returns,
          factors,
          penalty_parameters,
          weighting_type = weighting_type,
          tuning_type = 'c',
          include_standard_errors = TRUE,
          one_stddev_rule = one_stddev_rule,
          n_folds = n_folds
        )

        expect_length(adaptive_ifrp$risk_premia, n_factors)
        expect_length(adaptive_ifrp$standard_errors, n_factors)
        expect_length(adaptive_ifrp$penalty_parameter, 1)


        adaptive_ifrp0 = OptimalAdaptiveIFRP(
          returns,
          factors,
          penalty_parameters = 0.,
          weighting_type = weighting_type,
          tuning_type = 'c',
          include_standard_errors = TRUE,
          one_stddev_rule = one_stddev_rule,
          n_folds = n_folds
        )

        expect_equal(
          matrix(adaptive_ifrp0$risk_premia, n_factors, 1),
          ifrp$risk_premia
        )

        if (weighting_type == 'c') {

          adaptive_ifrp1 = OptimalAdaptiveIFRP(
            returns,
            factors,
            penalty_parameters = .1,
            weighting_type = weighting_type,
            tuning_type = 'c',
            include_standard_errors = TRUE,
            one_stddev_rule = one_stddev_rule,
            n_folds = n_folds
          )

          expect_equal(
            matrix(adaptive_ifrp1$risk_premia, n_factors, 1),
            aifrp1
          )

        }
      }
    }
  }

  ### RV
  for (weighting_type in c('c', 'b', 'a', 'n')) {
    for (one_stddev_rule in c(TRUE, FALSE)) {
      for (n_train_observations in c(120, 240)) {
        for (n_test_observations in c(120, 240)) {
          for (roll_shift in c(10, 12)) {

            adaptive_ifrp = OptimalAdaptiveIFRP(
              returns,
              factors,
              penalty_parameters,
              weighting_type = weighting_type,
              tuning_type = 'r',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              n_train_observations = n_train_observations,
              n_test_observations = n_test_observations,
              roll_shift = roll_shift
            )

            expect_length(adaptive_ifrp$risk_premia, n_factors)
            expect_length(adaptive_ifrp$standard_errors, n_factors)
            expect_length(adaptive_ifrp$penalty_parameter, 1)


            adaptive_ifrp0 = OptimalAdaptiveIFRP(
              returns,
              factors,
              penalty_parameters = 0.,
              weighting_type = weighting_type,
              tuning_type = 'r',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              n_train_observations = n_train_observations,
              n_test_observations = n_test_observations,
              roll_shift = roll_shift
            )

            expect_equal(
              matrix(adaptive_ifrp0$risk_premia, n_factors, 1),
              ifrp$risk_premia
            )

            if (weighting_type == 'c') {

              adaptive_ifrp1 = OptimalAdaptiveIFRP(
                returns,
                factors,
                penalty_parameters = .1,
                weighting_type = weighting_type,
                tuning_type = 'r',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                n_train_observations = n_train_observations,
                n_test_observations = n_test_observations,
                roll_shift = roll_shift
              )

              expect_equal(
                matrix(adaptive_ifrp1$risk_premia, n_factors, 1),
                aifrp1
              )


            }
          }
        }
      }
    }
  }

  # AdaptiveIFRP
  for (weighting_type in c('c', 'b', 'a', 'n')) {

    weights = intrinsicFRP:::AdaptiveWeightsCpp(
      returns, factors, weighting_type
    )

    expect_no_error(
      AdaptiveIFRP(
        returns,
        factors,
        penalty_parameters,
        weights
      )
    )

    expect_equal(
      AdaptiveIFRP(
        returns,
        factors,
        penalty_parameters = 0.,
        weights
      ),
      aifrp0
    )

    expect_equal(
      AdaptiveIFRP(
        returns,
        factors,
        penalty_parameters = 1e7,
        weights
      ),
      matrix(0., n_factors, 1)
    )

  }

})

test_that("Test OracleTFRP", {

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

  tfrp = TFRP(returns, factors, include_standard_errors = TRUE)

  expect_no_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE
  ))

  expect_error(OracleTFRP(
    t(returns),
    factors,
    penalty_parameters
  ))
  expect_error(OracleTFRP(
    returns,
    t(factors),
    penalty_parameters
  ))
  expect_error(OracleTFRP(
    t(returns),
    t(factors),
    penalty_parameters
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters = c("r", "s")
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      weighting_type = 4
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      tuning_type = 4
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      one_stddev_rule = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      gcv_scaling_n_assets = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      n_folds = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      n_train_observations = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      n_test_observations = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      roll_shift = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      check_arguments = "r"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      check_arguments = "r",
      target_level_kp2006_rank_test = -1.
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      check_arguments = "r",
      target_level_kp2006_rank_test = 1.5
  ))

  expect_equal(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      gcv_identification_check = FALSE,
      target_level_kp2006_rank_test = 0.05
    )$risk_premia,
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      gcv_identification_check = TRUE,
      target_level_kp2006_rank_test = 0.5
    )$risk_premia
  )

  weights =  intrinsicFRP:::AdaptiveWeightsCpp(returns, factors, 'c')

  oracle_tfrp0 = OracleTFRP(
    returns,
    factors,
    penalty_parameters = 0.
  )

  expect_equal(oracle_tfrp0$risk_premia, tfrp$risk_premia)

  oracle_tfrp1 = OracleTFRP(
    returns,
    factors,
    penalty_parameters = .1
  )

  oracle_tfrpm =  OracleTFRP(
    returns,
    factors,
    penalty_parameters
  )

  expect_length(oracle_tfrp0$risk_premia, n_factors)
  expect_length(oracle_tfrp1$risk_premia, n_factors)

  ### GCV
  for (weighting_type in c('c', 'b', 'a', 'n')) {
    for (one_stddev_rule in c(TRUE, FALSE)) {
      for (gcv_scaling_n_assets in c(TRUE, FALSE)) {
        for (gcv_identification_check in c(TRUE, FALSE)) {
          for (relaxed in c(TRUE, FALSE)) {

            oracle_tfrp = OracleTFRP(
              returns,
              factors,
              penalty_parameters,
              weighting_type = weighting_type,
              tuning_type = 'g',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              gcv_scaling_n_assets = gcv_scaling_n_assets,
              gcv_identification_check = gcv_identification_check,
              relaxed = relaxed
            )

            expect_length(oracle_tfrp$risk_premia, n_factors)
            expect_length(oracle_tfrp$standard_errors, n_factors)
            expect_length(oracle_tfrp$penalty_parameter, 1)


            oracle_tfrp00 = OracleTFRP(
              returns,
              factors,
              penalty_parameters = 0.,
              weighting_type = weighting_type,
              tuning_type = 'g',
              include_standard_errors = TRUE,
              one_stddev_rule = one_stddev_rule,
              gcv_scaling_n_assets = gcv_scaling_n_assets,
              gcv_identification_check = gcv_identification_check,
              relaxed = relaxed
            )

            expect_equal(
              matrix(oracle_tfrp00$risk_premia, n_factors, 1),
              tfrp$risk_premia
            )

            if (weighting_type == 'c') {

              oracle_tfrp1 = OracleTFRP(
                returns,
                factors,
                penalty_parameters = .1,
                weighting_type = weighting_type,
                tuning_type = 'g',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                gcv_scaling_n_assets = gcv_scaling_n_assets,
                gcv_identification_check = gcv_identification_check
              )

              expect_equal(
                oracle_tfrp1$risk_premia,
                oracle_tfrp1$risk_premia
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
        for (relaxed in c(TRUE, FALSE)) {

          oracle_tfrp = OracleTFRP(
            returns,
            factors,
            penalty_parameters,
            weighting_type = weighting_type,
            tuning_type = 'c',
            include_standard_errors = TRUE,
            one_stddev_rule = one_stddev_rule,
            n_folds = n_folds,
            relaxed = relaxed
          )

          expect_length(oracle_tfrp$risk_premia, n_factors)
          expect_length(oracle_tfrp$standard_errors, n_factors)
          expect_length(oracle_tfrp$penalty_parameter, 1)


          oracle_tfrp00 = OracleTFRP(
            returns,
            factors,
            penalty_parameters = 0.,
            weighting_type = weighting_type,
            tuning_type = 'c',
            include_standard_errors = TRUE,
            one_stddev_rule = one_stddev_rule,
            n_folds = n_folds,
            relaxed = relaxed
          )

          expect_equal(
            matrix(oracle_tfrp00$risk_premia, n_factors, 1),
            tfrp$risk_premia
          )

          if (weighting_type == 'c') {

            oracle_tfrp11 = OracleTFRP(
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
              oracle_tfrp11$risk_premia,
              oracle_tfrp1$risk_premia
            )

          }
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
            for (relaxed in c(TRUE, FALSE)) {

              oracle_tfrp = OracleTFRP(
                returns,
                factors,
                penalty_parameters,
                weighting_type = weighting_type,
                tuning_type = 'r',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                n_train_observations = n_train_observations,
                n_test_observations = n_test_observations,
                roll_shift = roll_shift,
                relaxed = relaxed
              )

              expect_length(oracle_tfrp$risk_premia, n_factors)
              expect_length(oracle_tfrp$standard_errors, n_factors)
              expect_length(oracle_tfrp$penalty_parameter, 1)


              oracle_tfrp00 = OracleTFRP(
                returns,
                factors,
                penalty_parameters = 0.,
                weighting_type = weighting_type,
                tuning_type = 'r',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                n_train_observations = n_train_observations,
                n_test_observations = n_test_observations,
                roll_shift = roll_shift,
                relaxed = relaxed
              )

              expect_equal(
                matrix(oracle_tfrp00$risk_premia, n_factors, 1),
                tfrp$risk_premia
              )

              if (weighting_type == 'c') {

                oracle_tfrp11 = OracleTFRP(
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
                  oracle_tfrp11$risk_premia,
                  oracle_tfrp1$risk_premia
                )

              }
            }
          }
        }
      }
    }
  }

})

test_that("Test OracleTFRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)
  n_observations = nrow(returns)

  # Calculating necessary statistics for the expected risk premia.
  covariance_factors_returns = stats::cov(factors, returns)
  variance_returns = stats::cov(returns)
  mean_returns = colMeans(returns)
  mean_factors = colMeans(factors)

  # Generating a range of penalty parameters for testing.
  penalty_parameters = c(
    0.,
    exp(seq(from=log(1.0e-8), to=log(1e-2), length.out=100)),
    seq(1e-2, 1e2, 100)
  )

  # Computing tradable factor risk premia for comparison.
  tfrp = TFRP(returns, factors, include_standard_errors = TRUE)

  # Testing basic functionality of OracleTFRP without errors including standard errors.
  expect_no_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE
  ))

  # Test if prewhite works
  expect_no_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE,
      hac_prewhite = TRUE
    ))

  # Testing error handling for incorrect dimensions (transposed matrices).
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

  # Testing error handling for invalid types or values of arguments.
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
  expect_error(
    OracleTFRP(
      c(),
      factors,
      penalty_parameters
  ))
  expect_error(
    OracleTFRP(
      returns,
      c(),
      penalty_parameters
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = "c"
  ))
  expect_error(
    OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE,
      hac_prewhite = "c"
  ))

  # Checking consistency of risk premia across different parameter settings.
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

  # Weights calculation based on the correlation between factors and returns.
  weights = 1. / colMeans(stats::cor(returns, factors)^2)

  # Comparing OracleTFRP with zero penalty parameters to TFRP.
  oracle_tfrp0 = OracleTFRP(
    returns,
    factors,
    penalty_parameters = 0.
  )

  expect_equal(oracle_tfrp0$risk_premia, tfrp$risk_premia)

  # Testing OracleTFRP with a non-zero penalty parameter.
  oracle_tfrp1 = OracleTFRP(
    returns,
    factors,
    penalty_parameters = .1
  )

  # Testing OracleTFRP across a range of penalty parameters.
  oracle_tfrpm =  OracleTFRP(
    returns,
    factors,
    penalty_parameters
  )

  # Validating the length of risk premia vectors.
  expect_length(oracle_tfrp0$risk_premia, n_factors)
  expect_length(oracle_tfrp1$risk_premia, n_factors)

  #### Additional tests for different GCV, CV, and RV settings.

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

  # Expect error if n_folds > nrow(returns)
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'c',
    n_folds = nrow(returns) + 1
  ))

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

  # expect error if n_train_observations == 0
  expect_error(OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      tuning_type = 'r',
      n_train_observations = 0
  ))
  # expect error if n_train_observations > n_observations
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'r',
    n_train_observations = nrow(returns) + 1
  ))
  # expect error if n_test_observations = 0
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'r',
    n_train_observations = 10,
    n_test_observations = 0
  ))
  # expect error if n_test_observations > nrow(returns) / 2
  expect_error(OracleTFRP(
      returns,
      factors,
      penalty_parameters,
      tuning_type = 'r',
      n_train_observations = 10,
      n_test_observations = nrow(returns) + 1
  ))
  # expect error if n_train_observations + n_test_observations > n_observations
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'r',
    n_train_observations = floow(nrow(returns) / 2) + 2,
    n_test_observations = floow(nrow(returns) / 2) + 2
  ))
  # expect_error if roll_shift = 0
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'r',
    n_train_observations = 10,
    n_test_observations = 10,
    roll_shift = 0
  ))
  # expect_error if roll_shift > n_test_observations
  expect_error(OracleTFRP(
    returns,
    factors,
    penalty_parameters,
    tuning_type = 'r',
    n_train_observations = 10,
    n_test_observations = 10,
    roll_shift = 15
  ))

})

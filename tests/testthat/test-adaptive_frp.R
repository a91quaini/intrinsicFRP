test_that("Test OptimalAdaptiveFRP", {

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

  frp = FRP(returns, factors, include_standard_errors = TRUE)

  expect_no_error(
    OptimalAdaptiveFRP(
      returns,
      factors,
      penalty_parameters,
      include_standard_errors = TRUE
  ))

  expect_error(OptimalAdaptiveFRP(
    t(returns),
    factors,
    penalty_parameters
  ))
  expect_error(OptimalAdaptiveFRP(
    returns,
    t(factors),
    penalty_parameters
  ))
  expect_error(OptimalAdaptiveFRP(
    t(returns),
    t(factors),
    penalty_parameters
  ))
  expect_error(
    OptimalAdaptiveFRP(
      returns,
      factors,
      penalty_parameters = c("r", "s")
  ))
  expect_error(
    OptimalAdaptiveIFRP(
      returns,
      factors,
      penalty_parameters,
      weighting_type = 4
    ))
  expect_error(
    OptimalAdaptiveIFRP(
      returns,
      factors,
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

  afrp0 = OptimalAdaptiveFRP(
    returns,
    factors,
    penalty_parameters = 0.
  )

  expect_equal(afrp0$risk_premia, frp$risk_premia)

  afrp1 = OptimalAdaptiveFRP(
    returns,
    factors,
    penalty_parameters = .1
  )

  expect_length(afrp0$risk_premia, n_factors)
  expect_length(afrp1$risk_premia, n_factors)

  ### GCV
  for (weighting_type in c('c', 'b', 'a', 'n')) {
    for (one_stddev_rule in c(TRUE, FALSE)) {
      for (gcv_vr_weighting in c(TRUE, FALSE)) {
        for (gcv_scaling_n_assets in c(TRUE, FALSE)) {
          for (gcv_identification_check in c(TRUE, FALSE)) {
            for (target_level_kp2006_rank_test in c(0., 0.005)) {

              adaptive_frp = OptimalAdaptiveFRP(
                returns,
                factors,
                penalty_parameters,
                weighting_type = weighting_type,
                tuning_type = 'g',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                gcv_vr_weighting = gcv_vr_weighting,
                gcv_scaling_n_assets = gcv_scaling_n_assets,
                gcv_identification_check = gcv_identification_check,
                target_level_kp2006_rank_test = target_level_kp2006_rank_test
              )

              expect_length(adaptive_frp$risk_premia, n_factors)
              expect_length(adaptive_frp$standard_errors, n_factors)
              expect_length(adaptive_frp$penalty_parameter, 1)


              adaptive_frp0 = OptimalAdaptiveFRP(
                returns,
                factors,
                penalty_parameters = 0.,
                weighting_type = weighting_type,
                tuning_type = 'g',
                include_standard_errors = TRUE,
                one_stddev_rule = one_stddev_rule,
                gcv_vr_weighting = gcv_vr_weighting,
                gcv_scaling_n_assets = gcv_scaling_n_assets,
                gcv_identification_check = gcv_identification_check,
                target_level_kp2006_rank_test = target_level_kp2006_rank_test
              )

              expect_equal(
                matrix(adaptive_frp0$risk_premia, n_factors, 1),
                frp$risk_premia
              )

              if (weighting_type == 'c') {

                adaptive_frp1 = OptimalAdaptiveFRP(
                  returns,
                  factors,
                  penalty_parameters = .1,
                  weighting_type = weighting_type,
                  tuning_type = 'g',
                  include_standard_errors = TRUE,
                  one_stddev_rule = one_stddev_rule,
                  gcv_vr_weighting = gcv_vr_weighting,
                  gcv_scaling_n_assets = gcv_scaling_n_assets,
                  gcv_identification_check = gcv_identification_check,
                  target_level_kp2006_rank_test = target_level_kp2006_rank_test
                )

                expect_equal(
                  matrix(adaptive_frp1$risk_premia, n_factors, 1),
                  afrp1$risk_premia
                )

              }

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

        adaptive_frp = OptimalAdaptiveFRP(
          returns,
          factors,
          penalty_parameters,
          weighting_type = weighting_type,
          tuning_type = 'c',
          include_standard_errors = TRUE,
          one_stddev_rule = one_stddev_rule,
          n_folds = n_folds
        )

        expect_length(adaptive_frp$risk_premia, n_factors)
        expect_length(adaptive_frp$standard_errors, n_factors)
        expect_length(adaptive_frp$penalty_parameter, 1)


        adaptive_frp0 = OptimalAdaptiveFRP(
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
          matrix(adaptive_frp0$risk_premia, n_factors, 1),
          afrp0$risk_premia
        )

        if (weighting_type == 'c') {

          adaptive_frp1 = OptimalAdaptiveFRP(
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
            matrix(adaptive_frp1$risk_premia, n_factors, 1),
            afrp1$risk_premia
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

            adaptive_frp = OptimalAdaptiveFRP(
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

            expect_length(adaptive_frp$risk_premia, n_factors)
            expect_length(adaptive_frp$standard_errors, n_factors)
            expect_length(adaptive_frp$penalty_parameter, 1)


            adaptive_frp0 = OptimalAdaptiveFRP(
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
              matrix(adaptive_frp0$risk_premia, n_factors, 1),
              frp$risk_premia
            )

            if (weighting_type == 'c') {

              adaptive_frp1 = OptimalAdaptiveFRP(
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
                matrix(adaptive_frp1$risk_premia, n_factors, 1),
                afrp1$risk_premia
              )


            }
          }
        }
      }
    }
  }

})


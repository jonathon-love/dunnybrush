
testthat::test_that('it works!', {

  data(testdata, package='dunnybrush')

  testthat::expect_error({
    monolm('notadataframe', Y ~ A + B)
  }, class='ValueError', regexp='Data must be a data frame')

  testthat::expect_error({
    monolm(testdata, ~ A + B)
  }, class='ValueError', regexp='Dependent variable must be specified')

  testthat::expect_error({
    monolm(testdata, Y ~ A + B)
  }, class='ValueError', regexp='Only full factorial models are supported')

  # results <- monolm(testdata, Y ~ A * B)

  data <- ToothGrowth
  # data <- subset(data, data$dose != 0.5)
  data$dose <- as.factor(data$dose)

  results <- monolm(data, len ~ supp * dose)

})

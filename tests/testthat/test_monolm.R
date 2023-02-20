
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

  expect_equal(results[1, 'Mean.Sq'], 3.458)
  expect_equal(results[2, 'Mean.Sq'], 1.135)
  expect_equal(results[3, 'Mean.Sq'], 6.012)
  expect_equal(results[4, 'Mean.Sq'], 2.062)
  expect_equal(results[5, 'Mean.Sq'], 0.865)

  expect_equal(results[3, 'p.value'], 0.01, tolerance=.03)
  expect_equal(results[4, 'p.value'], 0.07, tolerance=.05)
  data <- ToothGrowth
  # data <- subset(data, data$dose != 0.5)
  data$dose <- as.factor(data$dose)

  results <- monolm(data, len ~ supp * dose)

})

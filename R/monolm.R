
#' @importFrom jmvcore decomposeTerms
#' @importFrom rlang abort
#' @export
monolm <- function(data, formula, nsample=1000, nstep=Inf) {

  if ( ! inherits(data, 'data.frame'))
    abort('Data must be a data frame', class='ValueError')

  parts <- as.character(formula)
  if (length(parts) != 3)
    abort('Dependent variable must be specified in formula', class='ValueError')

  lhs <- as.character(formula)[2]
  dep <- data[[lhs]]

  terms <- attr(terms(formula), 'term.labels')
  decomposed <- decomposeTerms(terms)
  components <- unique(unlist(decomposed))

  ff.rhs <- paste(components, collapse='*')
  ff.formula <- as.formula(paste0('~', ff.rhs))
  ff.terms <- attr(terms(ff.formula), 'term.labels')

  if ( ! identical(ff.terms, terms))
    abort('Only full factorial models are supported at this time', class='ValueError')

  if (length(components) < 2)
    abort('At least two factors must be specified', class='ValueError')

  if (length(components) > 2)
    abort('Only two-way models are supported at this time', class='ValueError')

  for (name in components) {
    factor <- data[[name]]
    if ( ! is.factor(factor))
      abort('Explanatory variables must be factors', class='ValueError')
  }

  if ( ! is.numeric(dep))
    abort('Dependent variable must be numeric', class='ValueError')

  me.rhs <- paste(components, collapse='+')
  me.formula <- as.formula(paste0('~', me.rhs))

  group <- vapply(1:nrow(data), function(i) paste0(data[i,components], collapse=' '), '')

  linear.model <- stats::aov(formula=formula, data=data)

  design <- unique(model.matrix(linear.model))
  rownames(design) <- NULL

  g <- character()

  # construct g from the design matrix, so it can be used to
  # reorder things
  for (rowNo in 1:nrow(design)) {
    colNo <- 2
    allValues <- integer()
    for (comp in components) {
      value = 1
      for (levelNo in 1:(nlevels(data[[comp]])-1)) {
        value = value + design[rowNo, colNo] * levelNo
        colNo = colNo + 1
      }
      allValues <- c(allValues, value)
    }
    g <- c(g, paste(allValues, collapse=' '))
  }

  mse <- mean(tapply(dep, group, var))

  y <- tapply(dep, group, mean)
  y <- y[g]  # reorder to match design matrix

  counts <- xtabs(~., data.frame(group))
  counts <- counts[g]  # reorder to match design matrix

  w <- diag(counts/mse)

  intercept.model <- design[,1]

  A <- data[[ components[1] ]]
  A.levels <- nlevels(A)
  A.model <- design[,1:A.levels]

  B <- data[[ components[2] ]]
  B.levels <- nlevels(B)
  B.model <- design[,c(1,(A.levels+1):(A.levels+B.levels-1))]

  add.levels <- A.levels+B.levels
  add.model <- design[,c(1:(A.levels+B.levels-1))]

  AB.levels <- ncol(design)-ncol(add.model)
  AB.model <- design

  linear.model.anova <- stats::anova(linear.model)

  # calculate residual SS
  resid.SS <- linear.model.anova$'Sum Sq'[4]; resid.df <- linear.model.anova$'Df'[4];
  # fit monotonic linear models
  m <- mLR(data=y, weights=w, design=intercept.model, nstep=nstep); intercept.fit <- m$fit;
  m <- mLR(data=y, weights=w, design=A.model, nstep=nstep); A.fit <- m$fit;
  m <- mLR(data=y, weights=w, design=B.model, nstep=nstep); B.fit <- m$fit;
  m <- mLR(data=y, weights=w, design=add.model, nstep=nstep); add.fit <- m$fit;
  A.SS <- mse*(intercept.fit-A.fit); A.df <- A.levels-1; A.msq <- A.SS/A.df; A.F <- A.msq/mse
  A.pvalue <- pf(A.F, A.df, resid.df, lower.tail=FALSE)
  B.SS <- mse*(intercept.fit-B.fit); B.df <- B.levels-1; B.msq <- B.SS/B.df; B.F <- B.msq/mse
  B.pvalue <- pf(B.F, B.df, resid.df, lower.tail=FALSE);
  add.SS <- mse*(intercept.fit-add.fit);  add.df <- A.df+B.df; add.msq <- add.SS/add.df; add.F <- add.msq/mse
  AB.SS <- mse*add.fit; AB.df <- ncol(design)-add.levels+1; AB.msq <- AB.SS/AB.df; AB.F <- AB.msq/mse

  # calculate p-values for additive and interaction effects by double bootstrap
  add.emp.F <- c(); AB.emp.F <- c()
  for (isample in 1:nsample) {
    Y.boot <- bootstrap(dep, group, g) # first bootstrap of observed data
    # find p-value for additive effect
    y <- tapply(Y.boot, group, mean)
    y <- y[g]

    m.int <- mLR(data=y, weights=w, design=intercept.model, nstep=nstep) # fit intercept model
    # shift data to consistent with intercept model
    Y.int <- Y.boot
    for (i in 1:length(y)){
      j <- group==g[i]
      Y.int[j] <- Y.boot[j] - y[i] + m.int$pred[i]
    }

    Y.int.boot <- bootstrap(Y.int, group, g) # second bootstrap of predicted data under intercept model
    mse.int <- mean(tapply(Y.int.boot, group, var));
    y.int <- tapply(Y.int.boot, group, mean)
    y.int <- y.int[g]

    w.int <- w
    m.int.boot <- mLR(data=y.int, weights=w.int, design=intercept.model, nstep=nstep) # fit intercept model
    m.add.boot <- mLR(data=y.int, weights=w.int, design=add.model, nstep=nstep) # fit additive model
    add.emp.F <- append(add.emp.F,(m.int.boot$fit-m.add.boot$fit)/add.df) # store additive empirical F under intercept model
    # find p-value for interaction effect
    m.add <- mLR(data=y, weights=w, design=add.model, nstep=nstep) # fit additive model
    # shift data to consistent with additive model
    Y.add <- Y.boot
    for (i in 1:length(y)){
      j <- group==g[i]
      Y.add[j] <- Y.boot[j] - y[i] + m.add$pred[i]
    }
    Y.add.boot <- bootstrap(Y.add, group, g) # second bootstrap of additive model
    mse.add <- mean(tapply(Y.add.boot, group, var))

    y.add <- tapply(Y.add.boot, group, mean)
    y.add <- y.add[g]

    w.add <- w
    m.add.boot <- mLR(data=y.add, weights=w.add, design=add.model, nstep=nstep) # fit additive model
    AB.emp.F <- append(AB.emp.F,m.add.boot$fit/AB.df) # store interaction empirical F under additive model
  }
  add.pvalue <- sum(add.emp.F >= add.F)/length(add.emp.F) # calculate additive p-value
  AB.pvalue <- sum(AB.emp.F >= AB.F)/length(AB.emp.F) # calculate interaction p-value

  # store results
  rf <- 3
  vartable <- data.frame("Source"=c(components[1],components[2],paste(components, collapse="+"),paste(components, collapse=":"),"Residuals"),"Df"=c(A.df,B.df,add.df,AB.df,resid.df),
                         "Sum Sq"=c(round(A.SS,rf),round(B.SS,rf),round(add.SS,rf),round(AB.SS,rf),round(resid.SS,rf)),
                         "Mean Sq"=c(round(A.msq,rf),round(B.msq,rf),round(add.msq,rf),round(AB.msq,rf),round(mse,rf)),
                         "F value"=c(round(A.F,rf),round(B.F,rf),round(add.F,rf),round(AB.F,rf),NA),
                         "p value"=c(round(A.pvalue,rf+1),round(B.pvalue,rf+1),round(add.pvalue,rf+1),round(AB.pvalue,rf+1),NA))

  return(vartable)

}

bootstrap <- function (data, group, g) {
  yb <- data
  for (gv in g){
    j <- (group == gv)
    yb[j] <- sample(data[j],length(data[j]),replace=T)
  }
  return(yb)
}


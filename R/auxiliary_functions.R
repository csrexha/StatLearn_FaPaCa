#' Predict with an mboost model, split by base-learner
#'
#' @description
#' Predicts the linear predictor of a fitted mboost model for `newdata`, in total and
#' split by base-learner, together with the predicted probability, and attaches
#' identifier columns.
#'
#' @details
#' For [mboost::Binomial()] the linear predictor `eta` is on the half-log-odds scale, so
#' the predicted probability is `plogis(2 * eta)`. `eta` is the sum of `offset` and the
#' contributions of the base-learners in `bl` only; it equals the full linear predictor
#' (and `response` equals `plogis(2 * eta)`) only if `bl` contains every selected
#' base-learner.
#'
#' @param fit A fitted mboost model, e.g. from [fit_glmboost()] or [fit_gamboost()].
#' @param data A data.table with the identifier columns `cols` and one row per row of
#'   `newdata`, in the same order (typically `newdata` itself).
#' @param newdata Data to predict for, passed to [stats::predict()].
#' @param bl Character vector of base-learner names whose contributions are returned
#'   (`which` in [mboost::predict.mboost()]).
#' @param link Currently not used, kept for backward compatibility. `response` is always
#'   the prediction on the response scale of the model.
#' @param cols Character vector of columns of `data` to keep in the output.
#'
#' @returns A data.table with one row per row of `newdata` and the columns `cols`,
#'   `offset`, one column per base-learner in `bl`, `eta` (`offset` plus the
#'   base-learner contributions) and `response` (predicted probability).
#' @export
predict_result <- function(fit, data, newdata, bl, link = "logit",
							cols = c("ID", "Status", "FamilyID", "Clinical_Findings")){
  # prediction of each base leaner
  prd <- predict(fit, type = "link", newdata = newdata, which = bl)

  rsp <- predict(fit, type = "response", newdata = newdata)

  # logit link function
  if(is.character(link)){
  	link_fun <- make.link(link)$linkinv
  } else{
  	link_fun <- link
  }

  return(data.table(data[, ..cols],
                    offset = attr(prd, "offset"),
                    prd,
                    eta = attr(prd, "offset") + rowSums(prd),
                    response = as.vector(rsp)))
}

#' Summarise a variable by groups
#'
#' @description
#' Summarises a numeric variable by groups: number of observations, mean, standard
#' deviation, standard error of the mean and confidence interval of the mean.
#'
#' @param data A data.frame.
#' @param measurevar Name of the numeric column to summarise.
#' @param groupvars Character vector of the grouping columns. Default `NULL`, no grouping
#'   variables.
#' @param na.rm Whether to ignore missing values. Default `FALSE`.
#' @param conf.interval Confidence level of the interval. Default 0.95.
#' @param .drop Whether to drop combinations of the grouping variables without
#'   observations, passed to [plyr::ddply()]. Default `TRUE`.
#'
#' @returns A data.frame with one row per group: the grouping columns, `N`, the mean
#'   (named after `measurevar`), `sd`, `se` (standard error of the mean) and `ci` (half
#'   width of the confidence interval, from the t distribution with `N - 1` degrees of
#'   freedom).
#' @export
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                             c(N    = length2(xx[[col]], na.rm=na.rm),
                               mean = mean   (xx[[col]], na.rm=na.rm),
                               sd   = sd     (xx[[col]], na.rm=na.rm)
                               )
                             },
                         measurevar
                         )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

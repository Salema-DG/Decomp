#' Logit Decomposition
#'
#' @param .data
#' @param y
#' @param x
#' @param independent
#' @param ...
#'
#' @return A tibble with the with the contribution of each covariate to the change of the structural coefficient.
#' @import stringr,
#' @import dplyr,
#' @import magrittr,
#' @import rlang,
#' @import purrr,
#' @import CovTools,
#' @import colorDF,
#' @import tidyr
#' @export
#'
#' @examples
logit_decomposition <- function(.data, y, x, independent = FALSE,...) {

  .data %<>% rename(y = rlang::enquos(y) %>% paste() %>% str_remove("~"), #this is the way to retrieve characters from fucntions
                    x = rlang::enquos(x) %>% paste() %>% str_remove("~"))
  #full regression
  formula <- paste("y ~ x + ", enquos(...) %>% paste() %>% str_remove("~") %>% paste(collapse = "+")) %>% as.formula()
  reg_full <- glm(formula, data=.data, family=binomial(link = "logit"))

  #base regression
  formula <- "y ~ x" %>% as.formula()
  reg_base <- glm(formula, data=.data, family=binomial(link = "logit"))

  #Generalize the 1st auxiliary regressions
  aux <- rlang::enquos(...) %>% paste() %>% stringr::str_remove("~") %>%
    purrr::map(~{
      paste(.x, " ~ x") %>% as.formula() %>%
        lm(data=.data)})

  lf <- aux %>% length() #number of auxiliary variables/regressions

  #input the residuals to the dataset with purrr:
  res_df <- aux %>% purrr::map(~{.x$residuals}) %>%
    as.data.frame() #retrive the residuals
  names(res_df) <- paste("res_", 1:lf, sep="") #change the names to res_
  .data %<>% bind_cols(res_df) #join


  #RE regression: I need to use the coeffs
  formula <- paste("y ~ x +", paste("res_", 1:lf, collapse = "+", sep = "")) %>% as.formula()
  reg_RE <- glm(formula, data=.data, family=binomial(link = "logit"))


  #build second set of auxiliary regressions
  aux_2 <- paste("res_", 1:lf, sep="") %>%
    purrr::map(~{
      paste(.x, " ~ x + y", sep = "") %>% as.formula() %>%
        lm(data=.data)
    })

  #here is the part of calculating the divisibility of the logit bias.

  #1st stage: warrant z independence
  aux_coeffs_1 <- aux_2 %>% purrr::map(~{.x$coefficients[2]}) %>% unlist()
  aux_coeffs_2 <- aux_2 %>% purrr::map(~{.x$coefficients[3]}) %>% unlist()

  # varaince-covaraince matrix:
  inverse <- as.matrix(.data %>% dplyr::select(starts_with("res_"))) %>% CovTools::PreEst.glasso(parallel = TRUE, method = list(type="BIC", param =  c(0.001, 0.01,0.1,1,10,100, 1000, 10000))) #select also works with "" on variables



  #thus the assumed diagonal of the inverse matrix is calculable:
  a <- c()
  for(i in (1:lf)){
    i_list <- Map("*", inverse$C[ ,i][-i], aux_coeffs_2[1:lf][-i])
    i_list %<>% Reduce("+", .) %>% unlist() #here is the summation
    a <- c(a, i_list)
  }

  diag <- (reg_RE$coefficients[3:(2+lf)] /aux_coeffs_2[1:lf])-(a/aux_coeffs_2)

  if (independent == TRUE) { #considerando os z's independents
    diag <- (reg_RE$coefficients[3:(2+lf)] /aux_coeffs_2[1:lf])
  }



  uncorr <-  (aux_coeffs_1*aux_coeffs_2*diag) + aux_coeffs_2 * (reg_RE$coefficients[3:(2+lf)]-(diag*aux_coeffs_2))


  #building the decomposition table with a dataset and then xtable, colorDF and crayon package
  #explaning the correlated part: full and aux regressions
  aux_coeffs <- aux %>% purrr::map(~{.x$coefficients[2]}) %>% unlist()
  decomp <- tibble(variable = rlang::enquos(...) %>% paste() %>% stringr::str_remove("~"),
                   corr_part = aux_coeffs*reg_full$coefficients[3:(2+lf)],
                   uncorr_part= uncorr)


  #The summations
  decomp %<>% bind_rows(decomp %>% summarise(corr_part = sum(corr_part), uncorr_part = sum(uncorr_part)) %>% mutate(variable= "sum")) %>% dplyr::select(variable, corr_part, uncorr_part)
  decomp %<>% mutate(sum= corr_part+uncorr_part)


  #the "goals to explain": both for correlated and uncorrelated part
  #total bias
  total_bias <- reg_base$coefficients[2]-reg_full$coefficients[2]
  #explained with the first parts
  correlated_bias <- reg_RE$coefficients[2]-reg_full$coefficients[2]
  #left to explain is
  uncorrealted_bias <- reg_base$coefficients[2]-reg_RE$coefficients[2]


  decomp %<>% bind_rows(c(correlated_bias, uncorrealted_bias, total_bias) %>% as_tibble() %>% dplyr::rename(Biases="value") %>%
                          tibble::rownames_to_column() %>%
                          pivot_longer(-rowname) %>%
                          pivot_wider(names_from=rowname, values_from=value) %>% setNames(names(decomp))
  )



  return(print(colorDF::colorDF(decomp, theme = "wb") ,
               cat(crayon::bgCyan("Decomposition of logit model"))) )
}

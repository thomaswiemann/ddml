---
title: "Constructing a new ML Wrapper"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Constructing a new ML Wrapper}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE #requireNamespace("keras")
)
```

# Introduction

Load ``keras`` package.

```{r setup}
library(ddml)
library(keras)
```


Construct keras wrapper. The prediction methods for keras already outputs a numerical vector. It is thus not necessary to also construct a wrapper for the prediction method.

The first part of the wrapper constructs a neural network architecture. Below is a simple
example using the relu activation function and potential ell_1 regularization in the hidden or output layers. This part of the wrapper should be adjusted for more complicated architectures.

```{r}
mdl_keras <- function(y, X,
                      units = 10, nhidden = 1, lambda1 = 0, lambda2 = 0,
                      optimizer_fun = "rmsprop",
                      loss = "mse",
                      epochs = 10,
                      batch_size = min(1000, length(y)),
                      validation_split = 0,
                      callbacks = NULL,
                      steps_per_epoch = NULL,
                      metrics = c("mae"),
                      verbose = 0) {

  # ============================================================================
  # ADJUST FOR DIFFERENT ARCHITECTURES =========================================

  # Construct neural network architecture
  nnet <- keras_model_sequential()
  for (k in 1:nhidden) {
    nnet <- nnet %>%
      layer_dense(units = units, use_bias = T,
                  activation = "relu",
                  kernel_regularizer = regularizer_l1(l = lambda1))
  }#FOR
  nnet <- nnet %>%
    layer_dense(units = 1, use_bias = T,
                kernel_regularizer = regularizer_l1(l = lambda2))

  # ============================================================================
  # ADJUST FOR DIFFERENT ARCHITECTURES =========================================

  # Compile model
  nnet %>% keras::compile(optimizer = optimizer_fun,
                           loss = loss,
                           metrics = metrics)

  # Fit neural net
  nnet %>% keras::fit(X, y,
                       epochs = epochs,
                       batch_size = batch_size,
                       validation_split = validation_split,
                       callbacks = callbacks,
                       steps_per_epoch = steps_per_epoch,
                       verbose = verbose)

  # Amend class
  class(nnet) <- c("mdl_keras", class(nnet))

  # Return fit
  return(nnet)
}#MDL_KERAS
```


```{r}
predict.mdl_keras <- function(obj, newdata = NULL){
  # Check for new data
  #if(is.null(newdata)) newdata <- obj$X
  # Predict data and output as matrix
  class(obj) <- class(obj)[-1] # Not a pretty solution...
  as.numeric(predict(obj, newdata))
}#PREDICT.MDL_KERAS
```


Let's test this on a simple empirical example. I consider the BLP_1995 data and, for simplicity, estimate a partially linear model (rather than a partially linear iv model).

```{r}
y <- log(BLP_1995$share) / log(BLP_1995$outshr)
D <- BLP_1995$price
X <- cbind(BLP_1995$space, BLP_1995$hpwt)
```

Next, we need to specify the particular neural network learner we want to use. In addition to the architecture, it's important to properly tune the optimization algorithm. callback_list below helps with specifying learning rate adjustments and defines the early stopping rule.

```{r}

callbacks_list <- list(callback_early_stopping(monitor = "val_loss",
                                                            patience = 15,
                                                            restore_best_weights = T),
                                    callback_reduce_lr_on_plateau(monitor = "val_loss",
                                                                  factor = 1/10,
                                                                  patience = 10,
                                                                  verbose = F),
                                    callback_learning_rate_scheduler(
                                      function(epoch, learning_rate){
                                        if(epoch == 0) learning_rate <- 0.1
                                        return(learning_rate)
                                      }))

learners = list(what = mdl_keras,
                      args = list(units = 10, nhidden = 1, lambda1 = 0, lambda2 = 0,
                                  epochs = 50,
                                  verbose = F,
                                  validation_split = 0.1,
                                  callbacks = callbacks_list))
```

We can now run the ddml estimator. Note that because we have set verbose to TRUE in the above list of arguments, live training output is automatically plotted. This is very helpful for assessing whether the defined training parameters are sensible. Since it substantially slows down the training process, however, verbose should be set to FALSE whenever you're not planning to actively monitor the output.


```{r}
plm_fit <- ddml_plm(y, D, X, learners = learners)
plm_fit$coef
```

```{r}

learners = list(list(fun = mdl_keras,
                      args = list(units = 10, nhidden = 1, lambda1 = 0, lambda2 = 0,
                                  epochs = 50,
                                  verbose = F,
                                  validation_split = 0.1,
                                  callbacks = callbacks_list)),
                list(fun = mdl_keras,
                      args = list(units = 5, nhidden = 2, lambda1 = 0, lambda2 = 0,
                                  epochs = 50,
                                  verbose = F,
                                  validation_split = 0.1,
                                  callbacks = callbacks_list)))

plm_fit <- ddml_plm(y, D, X, learners = learners, shortstack = TRUE)
plm_fit$coef
```
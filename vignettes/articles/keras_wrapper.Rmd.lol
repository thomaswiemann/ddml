---
title: "Constructing a User-Provided Base Learner"
description: "Tutorial on writing a simple wrapper for new user-provided base learner."
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Constructing a User-Provided Base Learner}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("keras")
)
```

# Introduction

This article illustrates how users can write their own base learners. The necessary requirements for a base learner to be compatible with ``ddml`` are minimal:

1. An estimation function that takes an outcome vector ``y`` and a feature matrix ``X`` and returns a fitted model object as an S3 class
2. A ``predict`` method for the above S3 class that takes a feature matrix ``newdata`` and returns a numerical vector of predictions

We outline the requirements on two concrete examples:

1. A simple wrapper for ``stats::glm()``
2. A more extensive wrapper for neural networks based on [``keras``](https://cran.r-project.org/web/packages/keras/index.html)

To test our wrappers, we apply them to the included random subsample of 5,000 observations from the data of Angrist & Evans (1998). See ``?AE98`` for details.

```{r}
# Load ddml and set seed
library(ddml)
set.seed(75039)

# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
```

# A Wrapper for ``stats::glm()``

We begin with writing an estimation function that takes named arguments ``y`` and ``X`` and returns a fitted model object. To allow for arguments to be passed from the wrapper to ``stats::glm()``, we also make use of the ``...`` argument.

```{r}
# Write glm wrapper
mdl_glm <- function(y, X, ...) {
  df <- data.frame(y, X) # transform data from matrices to data.frame
  glm_fit <- glm(y ~ ., data = df, ...) # fit glm
  class(glm_fit) <- c("mdl_glm", class(glm_fit)) # append class
  return(glm_fit) # return fitted glm object
}#MDL_GLM
```

Note that the above wrapper first transforms the inputs to a ``data.frame``. This is necessary because ``glm()`` will not take matrices as input but ``ddml`` supports only matrices (and sparse matrices, see ``vignette("sparse")``). As a consequence, we also need to wrap ``predict.glm()``:

```{r}
# Write predict.glm wrapper
predict.mdl_glm <- function(object, newdata) {
  df <- data.frame(newdata) # transform data from matrices to data.frame
  predict.glm(object, df, type = "response")
}#PREDICT.MDL_GLM
```

To test the wrapper, we estimate the local average treatment effect:
```{r}
# Use the new glm base learner
learner_glm <-  list(what = mdl_glm,
                     args = list(family = binomial(link = "logit")))

# Estimate ddml_late
late_fit <- ddml_late(y, D, Z, X,
                      learners = learner_glm,
                      sample_folds = 10,
                      silent = TRUE)
summary(late_fit)
```

All works well! See ``?stats::glm()`` and ``?stats::predict.glm()`` for additional details.


# A Wrapper for Neural Networks based on ``keras``

We now consider an example using the neural network implementation of the ``keras`` package. See this [guide to keras basics](https://tensorflow.rstudio.com/guides/keras/basics) for an excellent introduction to neural networks in R.

```{r}
# Load keras
library(keras)
```

The requirements for the wrapper to be compatible with ``ddml`` are the same as before, however, because estimating a neural network is more involved, we need to write a more extensive wrapper. In particular, our wrapper for the ``keras`` package has four components:

1. Building the neural network architecture
2. Compiling the neural network object
3. Estimating the neural network
4. Appending the class of the returned object (to write a custom ``predict`` method)

The wrapper below allows for custom specification of a feed-forward neural network architecture using the relu activation function. Of course, much more involved architectures are also supported (feel free to experiment!).

```{r}
mdl_keras <- function(y, X,
                      units = 10, nhidden = 1, 
                      optimizer_fun = "rmsprop",
                      loss = "mse",
                      epochs = 10,
                      batch_size = min(1000, length(y)),
                      validation_split = 0,
                      callbacks = NULL,
                      steps_per_epoch = NULL,
                      metrics = c("mae"),
                      verbose = 0) {

  # ===============================================
  # ADJUST THIS PART FOR DIFFERENT ARCHITECTURES ==

  # Construct neural network architecture
  nnet <- keras_model_sequential()
  for (k in 1:nhidden) {
    nnet <- nnet %>%
      layer_dense(units = units, use_bias = T,
                  activation = "relu")
  }#FOR
  nnet <- nnet %>%
    layer_dense(units = 1, use_bias = T)

  # ===============================================
  # ===============================================

  # Compile neural net
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

  # Append class
  class(nnet) <- c("mdl_keras", class(nnet))

  # Return fitted object
  return(nnet)
}#MDL_KERAS
```

The prediction method for our fitted object ensures the predictions are returned as a numerical vector:

```{r}
predict.mdl_keras <- function(object, newdata){
  class(object) <- class(object)[-1] # Remove custom class from object
  as.numeric(predict(object, newdata))
}#PREDICT.MDL_KERAS
```


To test the wrapper, we again estimate the local average treatment effect. In addition to the architecture, it is important to properly tune the optimization algorithm. ``callback_list`` helps with specifying learning rate adjustments and defines the early stopping rule.

```{r}
# Specify callbacks
callbacks_list <- list(callback_early_stopping(monitor = "val_loss",
                                               patience = 15,
                                               restore_best_weights = T),
                       callback_reduce_lr_on_plateau(monitor = "val_loss",
                                                     factor = 1/10,
                                                     patience = 10,
                                                     verbose = F))

# Use the neural network base learner
learner_keras = list(what = mdl_keras,
                     args = list(units = 10, nhidden = 1,
                                 epochs = 100,
                                 verbose = F,
                                 validation_split = 0.1,
                                 callbacks = callbacks_list))

# Estimate ddml_late
late_fit <- ddml_late(y, D, Z, X,
                      learners = learner_keras,
                      sample_folds = 10,
                      silent = T)
summary(late_fit)
```

All works well!

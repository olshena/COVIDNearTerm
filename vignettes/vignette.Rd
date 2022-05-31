---
title: "COVIDNearTerm Demonstration"
author: "Kristopher Kapphahn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{COVIDNearTerm Demonstration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format:

- Never uses retina figures
- Has a smaller default figure size
- Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style

## Vignette Info

Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output:
      rmarkdown::html_vignette:
        css: mystyles.css

## Figures

The figure sizes have been customised so that you can easily put two images side-by-side.

```{r, fig.show='hold'}
plot(1:10)
plot(10:1)
```

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[A footnote here.], and tables, e.g. using `knitr::kable()`.

```{r, echo=FALSE, results='asis'}
knitr::kable(head(mtcars, 10))
```

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))

## 1 Software purpose

## 2 Software methods

## 3 Example data description

The covid_testing data frame contains daily test hospitalization figures in Santa Clara County from 2020/05/01 to 2020/09/26.

## 4 simulateAR example

```{r, fig.show='hold'}

library(COVIDNearTerm)

simulate_ar_object <- simulateAR(y = covid_hospitalizations$covid_total_7davg,
                              wsize = 14,
                              method = "equal",
                              pdays = 28,
                              nsim = 100,
                              seed = 5,
                              output_type = "predictions",
                              lambda = seq(0,1,0.05))

```

## 5 run buildAR on predictAR output

```{r, fig.show='hold'}

library(COVIDNearTerm)

simulate_ar_object <- simulateAR(y = covid_hospitalizations$covid_total_7davg,
                              wsize = 14,
                              method = "equal",
                              pdays = 28,
                              nsim = 100,
                              seed = 5,
                              output_type = "predictions",
                              lambda = seq(0,1,0.05))

```


## 6 Plot predictions with median

## 7 Plot predictions vs. test set

## 8 adjust parameters, method, and wsize and generate another plot like 7

# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(nmdtx)

test_check("nmdtx")

context("render_gene_card")

test_that("render_gene_card returns a valid HTML character", {
  output <- nmdtx:::render_gene_card('ENSG00000160201')
  expect_is(output, "character",
            info = "The output of render_gene_card should be a character string")

})

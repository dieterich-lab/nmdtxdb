# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(nmdtxdb)

test_check("nmdtxdb")

context("render_gene_card")

test_that("render_gene_card returns a valid HTML character", {
  output <- nmdtx:::render_gene_card("ENSG00000160201")
  expect_is(output, "character",
    info = "The output of render_gene_card should be a character string"
  )
})

# Assuming you have a valid gene_id for testing
valid_gene_id <- "ENSG00000160201"

# Test 1: Function should not throw an error with valid input
test_that("Function does not error with valid input", {
  expect_silent(render_gene_card(valid_gene_id))
})

# Test 2: Function should return valid HTML
test_that("Function returns valid HTML", {
  output <- render_gene_card(valid_gene_id)
  expect_true(grepl("<div", output)) # Simple check for HTML tag
  # Additional checks can be added to ensure the structure of the HTML is as expected
})

# Test 3: Function should return a helpful error for invalid input
test_that("Function returns helpful error for invalid input", {
  invalid_gene_id <- "ABC123"
  expect_error(render_gene_card(invalid_gene_id),
    regexp = "Input should start with 'ENS'"
  )
})

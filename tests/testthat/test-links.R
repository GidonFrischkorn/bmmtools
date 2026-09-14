# Tests for inverse_link(), written against
# local/dev/spec-milestone-1-score-layer.md section 1.
#
# The closed forms below are written out independently of the
# implementation on purpose: a test that calls the function it is testing
# to compute its own expectation cannot fail. The vocabulary and the
# closed forms come from bmm:::link_transform()
# (bmm/R/helpers-parameters.R:81-121, bmm 1.4.1.9000, read 2026-09-07);
# the agreement test near the bottom is what keeps the two from drifting.

# One test per link ------------------------------------------------------

test_that("identity returns its input unchanged", {
  x <- c(-2.5, 0, 1, 7.25)
  expect_equal(inverse_link(x, "identity"), x)
})

test_that("log inverts to exp", {
  x <- c(-1, 0, 2)
  expect_equal(inverse_link(x, "log"), c(exp(-1), 1, exp(2)))
})

test_that("softplus inverts to log1p(exp(x))", {
  x <- c(-3, 0, 4)
  expect_equal(inverse_link(x, "softplus"), log(1 + exp(x)))
  # softplus is strictly positive everywhere, which is why bmm uses it
  # for scale parameters
  expect_true(all(inverse_link(c(-50, -1, 0, 10), "softplus") > 0))
})

test_that("log1p inverts to expm1", {
  x <- c(-0.5, 0, 1.5)
  expect_equal(inverse_link(x, "log1p"), exp(x) - 1)
})

test_that("logm1 inverts to exp(x) + 1", {
  # Measured 2026-09-07 on brms 2.23.0: brms::expp1(x) is exp(x) + 1 and
  # brms::logm1(x) is log(x - 1). NOT log(exp(x) + 1), which the name
  # suggests and which a first draft of the spec assumed.
  x <- c(-1, 0, 2)
  expect_equal(inverse_link(x, "logm1"), exp(x) + 1)
  # The natural scale of a logm1 link is bounded below by 1. The bound is
  # >= rather than > in double precision: exp(-100) is ~3.7e-44, which is
  # below the spacing of doubles at 1, so exp(-100) + 1 is exactly 1.
  expect_true(all(inverse_link(c(-100, 0, 5), "logm1") >= 1))
  expect_true(all(inverse_link(c(-2, 0, 5), "logm1") > 1))
})

test_that("inverse is its own inverse", {
  x <- c(-4, 0.5, 2)
  expect_equal(inverse_link(x, "inverse"), 1 / x)
})

test_that("sqrt inverts to squaring", {
  x <- c(0, 1.5, 3)
  expect_equal(inverse_link(x, "sqrt"), x^2)
})

test_that("logit inverts to plogis", {
  x <- c(-2, 0, 2)
  expect_equal(inverse_link(x, "logit"), 1 / (1 + exp(-x)))
  expect_equal(inverse_link(0, "logit"), 0.5)
})

test_that("probit inverts to pnorm", {
  x <- c(-1.96, 0, 1.96)
  expect_equal(inverse_link(x, "probit"), stats::pnorm(x))
  expect_equal(inverse_link(0, "probit"), 0.5)
})

test_that("tan_half inverts to 2 * atan", {
  x <- c(-1, 0, 1)
  expect_equal(inverse_link(x, "tan_half"), 2 * atan(x))
  # the natural scale is a circular location in (-pi, pi)
  wide <- inverse_link(c(-1e6, 1e6), "tan_half")
  expect_true(all(abs(wide) < pi))
})

test_that("loglog inverts to exp(-exp(x))", {
  x <- c(-1, 0, 1)
  expect_equal(inverse_link(x, "loglog"), exp(-exp(x)))
  out <- inverse_link(x, "loglog")
  expect_true(all(out > 0 & out < 1))
})

test_that("cloglog inverts to 1 - exp(-exp(x))", {
  x <- c(-1, 0, 1)
  expect_equal(inverse_link(x, "cloglog"), 1 - exp(-exp(x)))
  out <- inverse_link(x, "cloglog")
  expect_true(all(out > 0 & out < 1))
})

# Round trips ------------------------------------------------------------

test_that("each link round-trips against its forward transform", {
  # forward transforms, written out from link_transform()'s
  # inverse = FALSE branch, with a domain-appropriate value per link
  forward <- list(
    identity = list(f = function(v) v, v = c(-2, 0, 3)),
    log = list(f = log, v = c(0.5, 1, 8)),
    softplus = list(f = function(v) log(expm1(v)), v = c(0.2, 1, 5)),
    log1p = list(f = log1p, v = c(-0.5, 0, 4)),
    logm1 = list(f = function(v) log(v - 1), v = c(1.5, 2, 9)),
    inverse = list(f = function(v) 1 / v, v = c(0.25, 2, 10)),
    sqrt = list(f = sqrt, v = c(0.25, 1, 16)),
    logit = list(f = stats::qlogis, v = c(0.1, 0.5, 0.9)),
    probit = list(f = stats::qnorm, v = c(0.1, 0.5, 0.9)),
    tan_half = list(f = function(v) tan(v / 2), v = c(-2, 0, 2)),
    loglog = list(f = function(v) log(-log(v)), v = c(0.1, 0.5, 0.9)),
    cloglog = list(f = function(v) log(-log1p(-v)), v = c(0.1, 0.5, 0.9))
  )

  for (link in names(forward)) {
    v <- forward[[link]]$v
    expect_equal(
      inverse_link(forward[[link]]$f(v), link),
      v,
      tolerance = 1e-8,
      info = paste("round trip failed for link", link)
    )
  }
})

test_that("the vocabulary is the twelve links bmm's link_transform() covers", {
  expect_setequal(
    link_vocabulary(),
    c(
      "identity", "log", "softplus", "log1p", "logm1", "inverse",
      "sqrt", "logit", "probit", "tan_half", "loglog", "cloglog"
    )
  )
})

# Agreement with bmm -----------------------------------------------------

test_that("inverse_link() agrees with bmm's own link_transform()", {
  # This is the test that catches divergence between the two tables. It
  # reaches through ::: deliberately: link_transform() is @noRd, and this
  # is a test rather than package code. Request B5 asks bmm to export it,
  # after which the ::: goes away.
  skip_if_not_installed("bmm")
  x <- c(-1.5, -0.25, 0.5, 1.75)
  for (link in link_vocabulary()) {
    expect_equal(
      inverse_link(x, link),
      bmm:::link_transform(x, link, inverse = TRUE),
      info = paste("bmmtools and bmm disagree on link", link)
    )
  }
})

test_that("logm1 agrees with brms::expp1()", {
  skip_if_not_installed("brms")
  x <- c(-1, 0, 2.5)
  expect_equal(inverse_link(x, "logm1"), brms::expp1(x))
})

# Edge cases -------------------------------------------------------------

test_that("NA is passed through at its own position", {
  x <- c(1, NA, 3)
  out <- inverse_link(x, "log")
  expect_equal(out, c(exp(1), NA, exp(3)))
  expect_true(is.na(out[2]))
  expect_length(out, 3L)
})

test_that("a zero-length input returns a zero-length numeric", {
  for (link in link_vocabulary()) {
    out <- inverse_link(numeric(0), link)
    expect_equal(out, numeric(0), info = paste("link", link))
  }
})

test_that("NULL link is a synonym for identity", {
  # matches link_transform()'s own `if (is.null(link)) link <- "identity"`
  x <- c(-1, 0, 2)
  expect_equal(inverse_link(x, NULL), x)
})

test_that("integer input is accepted and returns double", {
  out <- inverse_link(0L:2L, "log")
  expect_type(out, "double")
  expect_equal(out, exp(0:2))
})

test_that("values outside a link's domain are passed through, not NA-ed", {
  # Deliberate: silently substituting NA would corrupt a recovery metric
  # more quietly than an Inf a user can see (spec section 1, Errors).
  expect_equal(inverse_link(0, "inverse"), Inf)
  expect_equal(inverse_link(-2, "sqrt"), 4)
  expect_false(is.na(inverse_link(-800, "log")))
})

# Errors -----------------------------------------------------------------

test_that("an unknown link errors, naming the value and the vocabulary", {
  expect_error(inverse_link(1, "cauchit"), "cauchit")
  expect_error(inverse_link(1, "cauchit"), "identity")
})

test_that("a link of length != 1 errors", {
  expect_error(inverse_link(1, c("log", "logit")), "single string")
  expect_error(inverse_link(1, character(0)), "single string")
})

test_that("a non-character link errors", {
  expect_error(inverse_link(1, 3), "single string")
})

test_that("a non-numeric x errors", {
  expect_error(inverse_link("a", "log"), "numeric")
  expect_error(inverse_link(list(1, 2), "log"), "numeric")
})

test_that("a logical x errors rather than coercing silently", {
  # TRUE/FALSE reaching a link transform means something upstream went
  # wrong; coercing to 1/0 would hide it.
  expect_error(inverse_link(c(TRUE, FALSE), "log"), "numeric")
})

# link_of (spec 5, naming D27) -------------------------------------------

test_that("link_of finds a term's link by name, then by prefix", {
  links <- c(kappa = "log", kappa2 = "softplus", thetat = "logit")

  expect_equal(link_of("kappa", links), "log")
  expect_equal(link_of("kappa_task1", links), "log")
  # the longest prefix wins, so kappa2's cells do not take kappa's link
  expect_equal(link_of("kappa2_task1", links), "softplus")
  expect_equal(link_of("G", links), "identity")
  # a prefix has to end at an underscore
  expect_equal(link_of("kappax", links), "identity")
  expect_equal(link_of("kappa", NULL), "identity")
  # a list link table, as a bmm model carries it, works as well
  expect_equal(link_of("thetat_task2", as.list(links)), "logit")
  expect_equal(
    link_of(c("kappa", "kappa2_x", "G"), links),
    c("log", "softplus", "identity")
  )
})

# Compare a fit's estimates with a reference

The question a model answers last: does it agree with what is already
known? `cross_check()` puts a fit's estimates beside a **reference** — a
closed form such as `bmm::sdt_d()`, another implementation, or published
values — and returns the comparison in the tibble shape every other
bmmtools scorer returns.

## Usage

``` r
cross_check(
  fit,
  reference,
  scale = c("natural", "link"),
  links = NULL,
  level = c("population", "subject"),
  group = "id",
  ci_level = 0.95,
  ...
)
```

## Arguments

- fit:

  A fit with an
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  method, so a `brmsfit` and therefore a `bmmfit`. One fit: a list of
  them is an error pointing at
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md),
  which scores replications against one truth and separates them with a
  `replication` column.

- reference:

  The right-hand side. A data frame with `term` and `estimate`,
  optionally `ci_low` and `ci_high` (both or neither), `source` and, at
  the subject level, `id`; its values are read on `scale`. Or another
  fit, which is extracted at the same level, transformed the same way,
  and marked `source = "fit"`.

- scale:

  `"natural"` compares on the scale a reader interprets, `"link"` on the
  scale the model was estimated on. The fit's estimates and interval
  bounds go through
  [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
  exactly as
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  transforms them, so a natural-scale `cross_check()` and a
  natural-scale
  [`recover()`](https://www.gfrischkorn.org/bmmtools/reference/recover.md)
  on the same fit agree to the last digit. A data-frame `reference` is
  taken to be on `scale` already and is not transformed.

- links:

  A named character vector mapping a term to one of the link names
  [`inverse_link()`](https://www.gfrischkorn.org/bmmtools/reference/inverse_link.md)
  understands. `NULL` reads the link table from a `bmmfit`; with no
  table available the comparison falls back to the link scale and says
  so.

- level:

  `"population"` or `"subject"`. `"subject"` needs an `id` on both
  sides. `"sd"` and `"cor"` are not compared in this version: a
  correlation reference would need `var1`, `var2` and a pair-ordering
  rule of its own.

- group:

  The grouping variable subject-level estimates come from, passed to
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).

- ci_level:

  The interval mass, passed to
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md).

- ...:

  Passed to
  [`extract_estimates()`](https://www.gfrischkorn.org/bmmtools/reference/extract_estimates.md)
  for both sides: `converged`, `ci_method`.

## Value

A `bmmtools_cross_check` object: a tibble subclass with the columns
`term`, `estimate`, `ci_low`, `ci_high`, `ci_method`, `ci_level`,
`rhat`, `ess_bulk`, `ess_tail`, `reference`, `ref_low`, `ref_high`,
`source`, `bias`, `covered`, `overlap`, `scale`, `level`, `id` and
`converged`. `ref_low`, `ref_high` and `overlap` are `NA` for a
reference without intervals. Call
[`summary()`](https://rdrr.io/r/base/summary.html) on it for the
per-parameter metrics and
[`plot_recovery()`](https://www.gfrischkorn.org/bmmtools/reference/plot_recovery.md)
for the picture.

## Details

**The reference is a comparison, not a truth.** `bias` is the signed
difference from it and `covered` says the fit's interval contains it;
neither asserts that the reference is correct. Where the two disagree,
the reference is as much a candidate for the error as the fit is, and a
reference derived under assumptions the fit does not make (an
equal-variance d′ against an unequal-variance fit, say) will disagree
for reasons that are nobody's bug.

A reference term the fit did not estimate is a warning listing what the
fit has, and is dropped; a fit term the reference does not mention is
dropped in silence, because a fit routinely estimates more than a
reference covers. No term matching at all is an error rather than an
empty result. At the subject level the same holds for ids.

**d′ and d_a.** `bmm::sdt_d()` returns the equal-variance d′. bmm's
`sdt_yn` model estimates d_a = √2·δ/√(1 + r²), which equals d′ only when
`sdratio` is fixed at 0. Against a fit with a free `sdratio` the
reference is converted by the user,
`d_a = d' * sqrt(2 / (1 + exp(sdratio)^2))`; `cross_check()` compares
what it is given and does not guess a conversion, because the right one
depends on which `sdratio` the fit used.

## Examples

``` r
if (FALSE) { # \dontrun{
# a signal-detection fit against the closed form, equal-variance case
fit <- bmm::bmm(
  bmm::bmf(d ~ 1, criterion ~ 1),
  data, bmm::sdt_yn(response = "resp", stimulus = "stim", n_trials = "n")
)
reference <- tibble::tibble(
  term = c("d", "criterion"),
  estimate = c(
    bmm::sdt_d(hit_rate = 0.8, fa_rate = 0.2),
    bmm::sdt_criterion(hit_rate = 0.8, fa_rate = 0.2)
  ),
  source = "closed form"
)
cross_check(fit, reference)

# with sdratio free the fit estimates d_a, so the reference converts
sdratio <- 0.2
reference$estimate[1] <- reference$estimate[1] *
  sqrt(2 / (1 + exp(sdratio)^2))
cross_check(fit, reference)
} # }
```

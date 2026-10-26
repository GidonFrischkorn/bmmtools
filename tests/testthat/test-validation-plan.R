# validation_defaults() and validation_plan() (Milestone 26, spec 26
# § Functions). A plan fits nothing, so none of these tests needs bmm,
# except the ones that build a model from its name.

#' A bmmodel-shaped stand-in of the class an adapter is looked up by
#' @noRd
plan_model <- function(name = "mixture2p",
                       links = list(kappa = "log", thetat = "logit")) {
  structure(list(links = links), class = c("bmmodel", "circular", name))
}

# validation_defaults() --------------------------------------------------

test_that("validation_defaults lists the mixture2p and ezdm entries", {
  entries <- validation_defaults()
  expect_named(entries, c("mixture2p", "ezdm"))
  fields <- c(
    "model_args", "pars", "sds", "cors", "grid", "reps", "sampler", "sbc",
    "prior_check", "reference", "owner", "version", "sections"
  )
  for (entry in entries) expect_true(all(fields %in% names(entry)))
})

test_that("an entry holds the standard grid of D81", {
  entry <- validation_defaults("mixture2p")
  expect_identical(entry$reps, 10L)
  expect_identical(nrow(entry$grid), 6L)
  expect_setequal(entry$grid$n_subjects, c(20L, 50L))
  expect_setequal(entry$grid$n_trials, c(30L, 100L, 300L))
  expect_identical(
    entry$sampler[c("chains", "iter", "warmup")],
    list(chains = 4L, iter = 2000L, warmup = 1000L)
  )
  expect_identical(entry$sbc$n_sims, 100L)
  expect_identical(c(entry$sbc$n_subjects, entry$sbc$n_trials), c(20L, 50L))
  expect_identical(
    entry$sbc$sampler[c("chains", "iter")],
    list(chains = 2L, iter = 500L)
  )
})

test_that("the mixture2p entry carries the truth and its SBC prior", {
  entry <- validation_defaults("mixture2p")
  expect_equal(entry$pars, c(kappa = log(8), thetat = stats::qlogis(0.75)))
  expect_equal(entry$sds, c(kappa = 0.3, thetat = 0.5))
  expect_null(entry$cors)
  # the article's narrowed prior, because bmm's default SD prior on log
  # kappa draws values rmixture2p() cannot simulate from
  expect_identical(
    entry$sbc$prior$prior,
    c("normal(2, 0.5)", "normal(0, 0.3)", "normal(0, 0.5)")
  )
  expect_identical(entry$sbc$prior$nlpar, c("kappa", "kappa", "thetat"))
  expect_match(entry$sbc$prior_reason, "stack")
  expect_identical(entry$owner, "Gidon Frischkorn")
})

test_that("the ezdm entry carries Milestone 31's generating values", {
  entry <- validation_defaults("ezdm")
  expect_equal(entry$pars, c(drift = 1.5, bound = 0, ndt = log(0.3)))
  expect_equal(entry$sds, c(drift = 0.3, bound = 0.2, ndt = 0.2))
  expect_identical(dimnames(entry$cors)[[1L]], c("drift", "bound", "ndt"))
  expect_equal(entry$cors["drift", "bound"], 0.4)
  expect_equal(entry$cors["bound", "drift"], 0.4)
  expect_equal(entry$cors["drift", "ndt"], 0)
  expect_identical(diag(entry$cors), c(drift = 1, bound = 1, ndt = 1))
  expect_null(entry$sbc$prior)
})

test_that("the ezdm prior check judges mean RT against a plausible range", {
  # bmm 1.4.1's posterior_predict() gives mean RT only, and prior_check()
  # knows no range for ezdm, so the entry gives one and says why
  entry <- validation_defaults("ezdm")
  expect_identical(entry$prior_check$range, c(0.15, 5))
  expect_match(entry$prior_check$range_reason, "plausib")
  expect_match(entry$prior_check$range_reason, "not a bound")
  expect_null(validation_defaults("mixture2p")$prior_check$range_reason)
  plan <- validation_plan("ezdm", seed = 1)
  expect_identical(plan$prior_check$range, c(0.15, 5))
  expect_identical(
    plan$prior_check$range_reason, entry$prior_check$range_reason
  )
})

test_that("validation_defaults takes a name or a model and is NULL without", {
  expect_identical(
    validation_defaults(plan_model()),
    validation_defaults("mixture2p")
  )
  expect_null(validation_defaults("sdm"))
  expect_null(validation_defaults(plan_model("toy")))
  expect_error(validation_defaults(1), "bmmodel")
})

# validation_plan(): resolution ------------------------------------------

test_that("a plan fills every field from the entry", {
  plan <- validation_plan(plan_model(), seed = 1)
  entry <- validation_defaults("mixture2p")
  expect_s3_class(plan, "bmmtools_validation_plan")
  expect_identical(plan$model_name, "mixture2p")
  expect_match(plan$defaults, "^bmmtools .* entry$")
  expect_equal(plan$pars, entry$pars)
  expect_equal(plan$sds, entry$sds)
  expect_identical(plan$grid, entry$grid)
  expect_identical(plan$reps, entry$reps)
  expect_identical(plan$sampler, entry$sampler)
  expect_identical(plan$sbc, entry$sbc)
  expect_identical(plan$prior_check, entry$prior_check)
  expect_identical(plan$owner, "Gidon Frischkorn")
  expect_identical(plan$version, "0.1.0")
  expect_identical(plan$criteria, default_validation_criteria())
  expect_null(plan$convergence)
  expect_null(plan$reference)
})

test_that("a given field replaces the entry's and a list merges into it", {
  plan <- validation_plan(
    plan_model(),
    pars = c(kappa = 2, thetat = 0),
    reps = 3,
    sampler = list(chains = 2L),
    sbc = list(n_sims = 20L, prior = NULL),
    seed = 1
  )
  expect_equal(plan$pars, c(kappa = 2, thetat = 0))
  expect_identical(plan$reps, 3L)
  # the named entry changes and the rest of the budget stays
  expect_identical(plan$sampler$chains, 2L)
  expect_identical(plan$sampler$iter, 2000L)
  expect_identical(plan$sbc$n_sims, 20L)
  # an explicit NULL is bmm's default prior, not a deleted entry
  expect_true("prior" %in% names(plan$sbc))
  expect_null(plan$sbc$prior)
  expect_identical(plan$sbc$n_subjects, 20L)
})

test_that("a model without an entry needs pars and gets the generic grid", {
  model <- plan_model("toy", links = list(a = "identity"))
  expect_error(
    validation_plan(model, seed = 1),
    "`pars` is required: \"toy\" has no validation defaults"
  )
  plan <- validation_plan(model, pars = c(a = 0), seed = 1)
  expect_identical(plan$defaults, "none")
  expect_identical(plan$grid, standard_validation_grid())
  expect_null(plan$sds)
  expect_true(is.na(plan$owner))
  expect_true(is.na(plan$version))
  expect_match(format(plan)[[1L]], "version unversioned")
})

test_that("a name is built through bmm with the entry's variable names", {
  skip_if_not_installed("bmm")
  plan <- validation_plan("ezdm", seed = 1)
  expect_s3_class(plan$model, "bmmodel")
  expect_s3_class(plan$model, "ezdm")
  expect_identical(plan$model_name, "ezdm")
  expect_error(
    validation_plan("toy", pars = c(a = 0), seed = 1),
    "has no validation defaults, so it cannot be built"
  )
})

# validation_plan(): checks ------------------------------------------------

test_that("seed is required and must be one number", {
  expect_error(validation_plan(plan_model()), "`seed` is required")
  expect_error(validation_plan(plan_model(), seed = "a"), "single finite")
  expect_error(validation_plan(plan_model(), seed = c(1, 2)), "single finite")
})

test_that("a model that is neither a bmmodel nor a name is an error", {
  expect_error(validation_plan(list(), seed = 1), "bmmodel")
})

test_that("grid needs both columns and positive whole numbers", {
  expect_error(
    validation_plan(plan_model(), grid = data.frame(n_subjects = 2), seed = 1),
    "n_subjects.*n_trials"
  )
  expect_error(
    validation_plan(plan_model(),
      grid = data.frame(n_subjects = 2, n_trials = 0.5), seed = 1
    ),
    "positive whole numbers"
  )
  plan <- validation_plan(plan_model(),
    grid = expand.grid(n_subjects = c(5, 10), n_trials = 20), seed = 1
  )
  # expand.grid's out.attrs never reach the plan or its hash
  expect_identical(
    plan$grid,
    data.frame(n_subjects = c(5L, 10L), n_trials = c(20L, 20L))
  )
})

test_that("pars and sds must be named numbers", {
  expect_error(validation_plan(plan_model(), pars = 1, seed = 1), "named")
  expect_error(
    validation_plan(plan_model(), pars = c(a = NA_real_), seed = 1),
    "missing"
  )
  expect_error(
    validation_plan(plan_model(), sds = c(kappa = "a"), seed = 1),
    "named numeric"
  )
})

test_that("a reference must carry its provenance", {
  reference <- tibble::tibble(term = "kappa", estimate = 8)
  expect_error(
    validation_plan(plan_model(), reference = reference, seed = 1),
    "lacks the columns"
  )
  full <- tibble::tibble(
    term = "kappa", estimate = 8, package = "mixtur", version = "1.2.3",
    call = "fit_mixtur()", data_hash = "abc", date = "2026-10-24",
    source = "script"
  )
  plan <- validation_plan(plan_model(), reference = full, seed = 1)
  expect_identical(plan$reference, full)
  expect_identical(
    plan$stages$status[plan$stages$stage == "agreement"], "planned"
  )
})

test_that("criteria merge by name and unknown names are errors", {
  plan <- validation_plan(plan_model(),
    criteria = list(coverage_level = 0.9, tolerance = c(rmse = 0.01)),
    seed = 1
  )
  expect_identical(plan$criteria$coverage_level, 0.9)
  expect_identical(plan$criteria$prior_bound_rate_max, 0.10)
  expect_identical(plan$criteria$tolerance[["rmse"]], 0.01)
  expect_identical(plan$criteria$tolerance[["default"]], 0)
  expect_error(
    validation_plan(plan_model(), criteria = list(nope = 1), seed = 1),
    "no criterion"
  )
  expect_error(
    validation_plan(plan_model(),
      criteria = list(coverage_level = "a"), seed = 1
    ),
    "single number"
  )
})

test_that("convergence is checked against check_convergence's arguments", {
  plan <- validation_plan(plan_model(),
    convergence = list(rhat_max = 1.01), seed = 1
  )
  expect_identical(plan$convergence, list(rhat_max = 1.01))
  expect_error(
    validation_plan(plan_model(), convergence = list(nope = 1), seed = 1),
    "no argument"
  )
})

test_that("a list argument must be a named list", {
  expect_error(
    validation_plan(plan_model(), sampler = 4, seed = 1),
    "named list"
  )
  expect_error(
    validation_plan(plan_model(), sbc = list(n_sims = 0), seed = 1),
    "sbc\\$n_sims"
  )
  expect_error(
    validation_plan(plan_model(), prior_check = list(n_trials = -1), seed = 1),
    "prior_check\\$n_trials"
  )
  expect_error(
    validation_plan(plan_model(), sbc = list(prior = "normal(0, 1)"), seed = 1),
    "brmsprior"
  )
})

# validation_plan(): stages ------------------------------------------------

test_that("the stage table lists the eleven stages with their status", {
  plan <- validation_plan(plan_model(), seed = 1)
  stages <- plan$stages
  expect_identical(stages$stage, validation_stage_names())
  expect_identical(stages$order, seq_len(11L))
  status <- stats::setNames(stages$status, stages$stage)
  expect_identical(
    unname(status[c("run_record", "prior_predictive", "sbc", "recovery")]),
    rep("planned", 4L)
  )
  # Milestone 22 has landed, so identifiability runs
  expect_identical(status[["identifiability"]], "planned")
  expect_identical(status[["individual_differences"]], "planned")
  # #17 has landed (Milestone 24), so data requirements runs
  expect_identical(status[["data_requirements"]], "planned")
  expect_identical(status[["model_recovery"]], "not built")
  expect_identical(status[["posterior_checks"]], "not built")
  expect_identical(status[["agreement"]], "not applicable")
  reason <- stats::setNames(stages$reason, stages$stage)
  expect_identical(reason[["data_requirements"]], NA_character_)
  expect_identical(reason[["model_recovery"]], "not run (#6)")
  expect_identical(reason[["posterior_checks"]], "not run (#20)")
  expect_identical(reason[["agreement"]], "no reference implementation")
  required <- stats::setNames(stages$required, stages$stage)
  expect_false(any(required[c("agreement", "model_recovery", "posterior_checks")]))
  expect_true(all(required[c("prior_predictive", "sbc", "recovery")]))
})

test_that("one subject-varying parameter makes individual differences n/a", {
  plan <- validation_plan(plan_model(), sds = c(kappa = 0.3), seed = 1)
  at <- plan$stages$stage == "individual_differences"
  expect_identical(plan$stages$status[at], "not applicable")
  expect_match(plan$stages$reason[at], "fewer than two")
  expect_identical(plan$re_cor, "none")
  expect_null(plan$correlations)

  two <- validation_plan(plan_model(), seed = 1)
  expect_identical(two$re_cor, "all")
  expect_identical(two$correlations, "model")
})

test_that("data requirements needs a varying parameter and trials to read", {
  status_of <- function(plan) {
    plan$stages$status[plan$stages$stage == "data_requirements"]
  }
  reason_of <- function(plan) {
    plan$stages$reason[plan$stages$stage == "data_requirements"]
  }
  # one subject-varying parameter has a person-level correlation to read
  one <- validation_plan(plan_model(), sds = c(kappa = 0.3), seed = 1)
  expect_identical(status_of(one), "planned")

  none <- validation_plan(plan_model(),
    sds = c(kappa = 0, thetat = 0), seed = 1
  )
  expect_identical(status_of(none), "not applicable")
  expect_identical(reason_of(none), "no subject-varying parameter")

  # D83 reads off at 50 subjects, which this grid does not have
  no_fifty <- validation_plan(plan_model(),
    grid = data.frame(n_subjects = 20L, n_trials = c(30L, 100L)), seed = 1
  )
  expect_identical(status_of(no_fifty), "not applicable")
  expect_identical(
    reason_of(no_fifty), "the grid has no design with 50 subjects"
  )

  # one number of trials at 50 subjects: there is nothing to read along
  flat <- validation_plan(plan_model(),
    grid = data.frame(
      n_subjects = c(20L, 50L, 20L), n_trials = c(30L, 30L, 100L)
    ),
    seed = 1
  )
  expect_identical(status_of(flat), "not applicable")
  expect_identical(
    reason_of(flat), "the grid has one number of trials at 50 subjects"
  )

  # the subjects it reads off at are a criterion
  moved <- validation_plan(plan_model(),
    grid = data.frame(n_subjects = 20L, n_trials = c(30L, 100L)),
    criteria = list(read_off_n_subjects = 20), seed = 1
  )
  expect_identical(status_of(moved), "planned")
  expect_identical(moved$criteria$read_off_n_subjects, 20)
  expect_identical(moved$criteria$read_off_r_low, 0.80)
  expect_identical(default_validation_criteria()$read_off_n_subjects, 50)
})

test_that("stages selects, always keeps the run record, and checks readers", {
  plan <- validation_plan(plan_model(),
    stages = c("prior_predictive", "sbc"), seed = 1
  )
  status <- stats::setNames(plan$stages$status, plan$stages$stage)
  expect_identical(status[["run_record"]], "planned")
  expect_identical(status[["prior_predictive"]], "planned")
  expect_identical(status[["recovery"]], "excluded")
  # a stage that does not exist stays "not built", not "excluded"
  expect_identical(status[["model_recovery"]], "not built")
  expect_error(
    validation_plan(plan_model(), stages = "nope", seed = 1),
    "names no stage"
  )
  expect_error(
    validation_plan(plan_model(), stages = "computation", seed = 1),
    "Add \"recovery\""
  )
})

# validation_plan(): seeds, hash, print, files -------------------------

test_that("every stage derives its seed from the master seed", {
  plan <- validation_plan(plan_model(), seed = 10)
  expect_identical(
    plan$seeds,
    c(recovery = 10, prior_predictive = 11, sbc = 12)
  )
})

test_that("the hash ignores seed, owner and version and nothing else", {
  plan <- validation_plan(plan_model(), seed = 1)
  same <- validation_plan(plan_model(),
    seed = 2, owner = "someone else", version = "9.9.9"
  )
  expect_identical(plan$hash, same$hash)
  expect_true(is.character(plan$hash) && nchar(plan$hash) > 8L)
  changed <- list(
    validation_plan(plan_model(), reps = 2, seed = 1),
    validation_plan(plan_model(), pars = c(kappa = 2, thetat = 1), seed = 1),
    validation_plan(plan_model(), sbc = list(n_sims = 50), seed = 1),
    validation_plan(plan_model(), criteria = list(z_mean_max = 3), seed = 1),
    validation_plan(plan_model(), stages = "sbc", seed = 1)
  )
  for (other in changed) expect_false(identical(other$hash, plan$hash))
})

test_that("a formula or a function enters the hash as its text", {
  generator <- function(pars, n_trials, model) NULL
  one <- validation_plan(plan_model(), generator = generator, seed = 1)
  # the same text, built in another environment
  local({
    generator <- function(pars, n_trials, model) NULL
    two <- validation_plan(plan_model(), generator = generator, seed = 1)
    expect_identical(two$hash, one$hash)
  })
})

test_that("print shows the model, the grid, the budget and the stages", {
  plan <- validation_plan(plan_model(), seed = 1)
  lines <- format(plan)
  expect_match(lines[[1L]], "mixture2p, version 0.1.0, hash ")
  expect_true(any(grepl("6 designs x 10 replications", lines)))
  expect_true(any(grepl("n_subjects 20, 50; n_trials 30, 100, 300", lines)))
  expect_true(any(grepl("SBC: 100 simulations at 20 subjects x 50", lines)))
  expect_true(any(grepl("the plan's own prior", lines)))
  expect_true(any(grepl("data_requirements +planned", lines)))
  expect_true(any(grepl("model_recovery +not built +not run \\(#6\\)", lines)))
  expect_output(print(plan), "<bmmtools validation plan>")
  expect_invisible(print(plan))
})

test_that("a plan round-trips through write_plan and read_plan", {
  plan <- validation_plan(plan_model(), seed = 1)
  file <- file.path(withr::local_tempdir(), "sub", "plan.rds")
  expect_identical(write_plan(plan, file), file)
  back <- read_plan(file)
  expect_identical(back, plan)
  expect_identical(plan_hash(back), plan$hash)

  other <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(a = 1), other)
  expect_error(read_plan(other), "not a validation plan")
  expect_error(read_plan(file.path(tempdir(), "none.rds")), "does not exist")
  expect_error(write_plan(list(), file), "must be a validation plan")
})

test_that("each malformed argument is an error naming it", {
  expect_error(
    validation_plan(plan_model(), reference = list(term = "a"), seed = 1),
    "`reference` must be a data frame"
  )
  expect_error(
    validation_plan(plan_model(), criteria = list(1), seed = 1),
    "`criteria` must be a named list"
  )
  expect_error(
    validation_plan(plan_model(),
      criteria = list(tolerance = "a"), seed = 1
    ),
    "criteria\\$tolerance"
  )
  expect_error(
    validation_plan(plan_model(), owner = c("a", "b"), seed = 1),
    "`owner` must be a single string"
  )
  expect_error(
    validation_plan(plan_model(), stages = 1, seed = 1),
    "character vector of stage names"
  )
  expect_error(write_plan(validation_plan(plan_model(), seed = 1), 1), "single path")
  expect_error(read_plan(c("a", "b")), "single path")
})

test_that("a formula enters the hash, and the print names a reference", {
  formula <- fake_bmmformula(kappa = kappa ~ 1, thetat = thetat ~ 1)
  plain <- validation_plan(plan_model(), seed = 1)
  given <- validation_plan(plan_model(), formula = formula, seed = 1)
  expect_false(identical(plain$hash, given$hash))
  expect_identical(plan_key(formula), formula_key(formula))

  reference <- tibble::tibble(
    term = "kappa", estimate = 8, package = "mixtur", version = "1.2.3",
    call = "fit()", data_hash = "abc", date = "2026-10-24", source = "script"
  )
  lines <- format(validation_plan(plan_model(),
    reference = reference,
    sampler = list(backend = NULL), seed = 1
  ))
  expect_true(any(lines == "Reference: mixtur"))
  # an explicit NULL entry prints as NULL, not as a missing value
  expect_true(any(grepl("backend = NULL", lines)))
  model <- plan_model("toy", links = list(a = "identity"))
  empty <- format(validation_plan(model, pars = c(a = 0), sampler = list(), seed = 1))
  expect_true(any(grepl("^Sampler: chains = 4", empty)))
})

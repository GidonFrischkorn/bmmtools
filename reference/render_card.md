# Render the validation card of a result

**\[experimental\]**

## Usage

``` r
render_card(result, output, format = "html", template = NULL, quiet = TRUE)
```

## Arguments

- result:

  A `bmmtools_validation`, its `result.rds`, or the directory it was run
  in.

- output:

  The HTML file to write. Its directory is created if missing, and an
  existing file is replaced.

- format:

  `"html"`, the only format so far.

- template:

  `NULL` for the card bmmtools ships
  (`system.file("card", "card.qmd", package = "bmmtools")`), or the path
  to a Quarto document of your own that takes the parameter `result`,
  the path to a `result.rds`.

- quiet:

  Passed to
  [`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html).

## Value

`output`, invisibly.

## Details

Renders the card of a
[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md)
result to one self-contained HTML file: the headline verdict, the run
record and the plan, the verdict of every stage beside its criterion, a
section per stage in the order of the workflow (results, or why the
stage did not run), the model's own sections, and the numbers that moved
since the earlier result the run was compared with.

The template prints no number of its own: every figure on the card is
read off the result (its `numbers`, `verdicts`, `stages`, `run_info` and
`plan`), so a card always describes the run it was rendered from.

## Requirements

The quarto R package and the Quarto CLI (<https://quarto.org>); both are
checked before anything is rendered.

## See also

[`run_validation()`](https://www.gfrischkorn.org/bmmtools/reference/run_validation.md),
[`verdict_table()`](https://www.gfrischkorn.org/bmmtools/reference/verdict_table.md),
[`card_numbers()`](https://www.gfrischkorn.org/bmmtools/reference/card_numbers.md),
[`compare_validations()`](https://www.gfrischkorn.org/bmmtools/reference/compare_validations.md).

## Examples

``` r
if (FALSE) { # \dontrun{
plan <- validation_plan("ezdm", seed = 2026)
result <- run_validation(plan, "validation/ezdm")
render_card(result, "validation/ezdm/card.html")
} # }
```

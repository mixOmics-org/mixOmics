<!--
Title: a plain sentence that says what changes, for example "Fix the legend order in plotIndiv".
Please read CONTRIBUTING.md before you open the pull request.
-->

## What changes and why

<!-- What the change does and the reason for it. Link the issue it closes, for example: Closes #123 -->

## How it was checked

- [ ] `devtools::test()` passes
- [ ] `devtools::check(args = "--no-manual", error_on = "warning")` passes
- [ ] Tests are added or updated for any change in behaviour

## Effect on results

<!-- Tick one. If results change, write under the boxes which functions, by how much, and what you checked the new results against. -->

- [ ] No numerical result and no plot changes
- [ ] Numerical results or plots change, as described below

## Checklist

- [ ] I have read CONTRIBUTING.md and the Code of Conduct
- [ ] A signed mixOmics contributor agreement covers this contribution (see CONTRIBUTING.md)
- [ ] Every commit is signed off (`git commit -s`), certifying the Developer Certificate of Origin
- [ ] The pull request covers one topic
- [ ] `NEWS.md` has a line for any change a user would notice
- [ ] Help pages are regenerated with `devtools::document()` if any roxygen comment, example file or template changed
- [ ] `README.md` is not edited by hand (edit `inst/README.Rmd` or `inst/README-WhatsNew.Rmd`)

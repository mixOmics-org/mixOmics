# Contributing to mixOmics

We welcome bug reports, feature suggestions and contributions to the mixOmics R package. You can report a bug or suggest an improvement in the repository’s Issues tab without signing an agreement. If you’d like to contribute a change to the package, follow the steps below.

## Contributor agreement

mixOmics grew out of Kim-Anh Lê Cao’s doctoral research. She leads and maintains the project. The contributor agreement gives her the rights to maintain and develop mixOmics and to decide how covered contributions are released and licensed over time.

Under the agreement, you assign the copyright you hold in your contribution to Kim-Anh Lê Cao personally, to the extent permitted by applicable law. Where and to the extent an assignment cannot take effect, you retain that copyright and grant her a worldwide, royalty-free licence for the full duration of the applicable rights, irrevocable to the extent permitted by law. The fallback licence allows her to use, copy, modify and distribute the contribution, and to sublicense it under open-source or proprietary terms. The signed agreement sets out the full terms, including applicable patent permissions and moral rights consents.

The `License` field in `DESCRIPTION` identifies the open-source licence for each version of the package. Recipients of earlier releases retain their rights under the terms of the licence that accompanied those releases. Contributors receive credit for their work.

## Before submitting a pull request

1. Email [contributors@mixomics.org](mailto:contributors@mixomics.org) with a brief description of your proposed contribution. We’ll send you the contributor agreement and signing instructions.
2. Tell us if the contribution was made in the course of employment, study or funded research, includes work by others, or contains third-party material. Identify any third-party material and its licence. If your employer, university or another organisation holds copyright in the work you are contributing, an authorised representative must sign the contributor agreement on its behalf for the rights it holds. If more than one person or organisation holds copyright in that work, each rights holder must sign for the rights they hold.
3. Submit your pull request once the signing requirements are settled. Include a short description of the change and how you tested it, where relevant.

A maintainer will verify that the signed agreement covers the contribution and that any additional permissions are in place before merging the pull request. If you’ve already signed, let us know if the rights holder or the circumstances in which you contribute have changed.

## Sign-off (Developer Certificate of Origin)

Every commit must be signed off, with a line at the end of its message:

```
Signed-off-by: Your Name <you@example.org>
```

With this line you certify the [Developer Certificate of Origin](DCO), version 1.1: that you wrote the change, or otherwise have the right to submit it under the open source licence indicated in the file. Use your own name and an email address you can be reached at; the sign-off is public and stays in the history of the repository. `git commit -s` adds the line for you, and `git commit --amend -s` adds it to your last commit if you forgot (if you had already pushed that commit, push again with `git push --force-with-lease`). If you edit a file on github.com, add the line to the commit message there. The sign-off says where your contribution comes from; the contributor agreement says what may be done with it. You need both.

## Licences

mixOmics is released under the GNU Affero General Public License, version 3 only ([LICENSE.md](LICENSE.md)). Code that came from other packages keeps its own licence.

Licences are named in two ways in this repository:

- **SPDX identifiers**, such as `AGPL-3.0-only` or `GPL-2.0-or-later`, come from the [SPDX License List](https://spdx.org/licenses/). SPDX is an open standard (ISO/IEC 5962:2021) that gives each licence one exact name. Every file in `R/` states its licence in a `SPDX-License-Identifier:` line.
- **R licence names**, such as `AGPL-3` or `GPL (>= 2)`, are the form R accepts in the `License` field of `DESCRIPTION` ([Writing R Extensions, Licensing](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Licensing)).

| SPDX identifier | R licence name |
|---|---|
| `AGPL-3.0-only` | `AGPL-3` |
| `GPL-2.0-or-later` | `GPL (>= 2)` |
| `GPL-2.0-only OR GPL-3.0-only` | `GPL-2 \| GPL-3` |
| `GPL-3.0-or-later` | `GPL (>= 3)` |

If you add code from another package, keep its licence and copyright notice, mark it with `SPDX-SnippetBegin` and `SPDX-SnippetEnd` as the existing files in `R/` do, and say so in the pull request.

## Questions about an analysis

Please ask on the discussion forum: <https://mixomics-users.discourse.group>. It is the place for questions about which method to use, how to set its parameters and how to read a result.

## Bugs and feature requests

Open an issue at <https://github.com/mixOmics-org/mixOmics/issues>. For a bug, please include:

- a small example that reproduces it, if possible with one of the datasets that come with the package;
- what you expected and what happened instead;
- the output of `sessionInfo()`.

If you think you have found a security problem, do not open an issue. Follow [SECURITY.md](SECURITY.md) instead.

## Proposing a change

1. For anything larger than a small fix, open an issue first, so that the approach can be agreed before you write the code.
2. Fork the repository and create a branch from `master`. If you have write access and create the branch in this repository, its name must start with `feat/`, `fix/`, `chore/`, `docs/`, `refactor/`, `test/` or `release/`. Other names are refused, apart from `master`, `main`, `gh-pages` and the names that GitHub and automated tools such as Dependabot create.
3. Make the change, and add or update tests for any change in behaviour.
4. Run the tests and the package check on your machine (see below).
5. Make sure a signed contributor agreement covers the change (see Before submitting a pull request), sign off every commit (see Sign-off, above), then open a pull request against `master`, and say what changed and why.

### Setting up

The README describes how to [set up a development environment](README.md#bug-reports-and-pull-requests).

### Code and documentation

- Follow [Bioconductor's guidance on R code](https://contributions.bioconductor.org/r-code.html).
- The help pages in `man/` and the `NAMESPACE` file are generated by roxygen2 from the comments in `R/`, the example files in `examples/` and the templates in `man-roxygen/`. Edit those sources and run `devtools::document()`; do not edit the generated files.
- `README.md` is rendered from `inst/README.Rmd` and `inst/README-WhatsNew.Rmd`. Edit those two files, not `README.md`. When a pull request into `master` changes them, or changes `DESCRIPTION`, the Render README workflow renders `README.md` and, if it has changed, commits it to your branch, so pull before you push again. From a fork the workflow cannot commit: if the rendered file differs, it fails, and you download the rendered `README.md` from the `rendered-readme` artifact of that run and commit it yourself.
- Add a line to `NEWS.md` for any change that a user would notice.

### Tests and checks

- The tests are in `tests/testthat` and use testthat. Plots are compared with the snapshots stored in `tests/testthat/_snaps`, using vdiffr.
- Run `devtools::test()` for the tests and `devtools::check(args = "--no-manual", error_on = "warning")` for the full package check, as CI runs it.
- On each push, and for each pull request against `master`, GitHub Actions runs R CMD check on Linux, macOS and Windows, and treats a warning as a failure. It uses the development version of Bioconductor, except on a branch named after a Bioconductor release, such as the existing `RELEASE_3_15`, where it uses the current release. It then runs BiocCheck. On a push, a commit message containing `skip-ci` or `worksave` skips this R CMD check run, and one containing `skip-check` skips its check steps. The pull request run is skipped when an earlier run on the same commit ended as a success. A separate workflow measures test coverage on pushes to `master` and on pull requests against it.

### Commit messages

Use the form `type(scope): summary` (the scope may be left out), for example `fix(citation): correct the author names and add the DOI`. The usual types are `feat`, `fix`, `docs`, `test`, `chore` and `refactor`.

## Review and merging

A pull request into `master` needs an approving review from a maintainer, given after the last push to it by someone other than the person who made that push, and every review conversation must be resolved. It is merged with a merge commit. The reviewers are researchers, so a review can take some weeks.

## Releases

mixOmics is released through [Bioconductor](https://bioconductor.org/packages/mixOmics). The maintainers send the changes on `master` to Bioconductor; contributors do not need to do anything for a release.

## Code of Conduct

Everyone taking part is expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

# Contributing to brutifus

If you have a question  about brutifus, [jump here](#asking-a-question).
If you want to report a bug with brutifus, [jump here](#reporting-a-bug).

If you are still reading this, you may actually be considering contributing to the development of
brutifus. :heart_eyes: :tada:

There are many ways that you can do so, including by:
- [reporting a bug](#reporting-a-bug)
- fixing a [known issue](https://github.com/brutifus/brutifus/issues?q=is%3Aissue+),
- implementing a new functionality, and/or
- improving the documentation:
  * in the code, with better docstrings
  * in this repository (for example this very file !)
  * in the website, via the docs `.rst` files

All these contributions are welcome, and what follows should help you get started. Note that contributing to brutifus does *not* necessarily require an advanced knowledge of Python and/or Github. Helping us fix typos in the docs, for example, could be an excellent first contribution. Plus, :anger: typos :anger: are the worst !

## Table of contents

- [Code of conduct](#code-of-conduct)
- [Reporting a bug](#reporting-a-bug)
- [Essential things to know about brutifus](#essential-things-to-know-about-brutifus)
- [Styles](#styles)
- [Step-by-step guide to contributing](#step-by-step-guide-to-contributing)

## Code of conduct

This project and everyone participating in it is governed by the [brutifus Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [frederic.vogt@alumni.anu.edu.au](mailto:frederic.vogt@alumni.anu.edu.au).

## Asking a question

The [brutifus Github Discussions page](https://github.com/brutifus/brutifus/discussions) is the place to ask your questions about brutifus, share your ideas for the code, show something cool you've done with brutifus, and generally engage with other brutifus users.

## Reporting a bug

If you find something odd/wrong/broken with brutifus, first check if it is a [known issue](https://github.com/brutifus/brutifus/issues?q=is%3Aissue+). If not, please create a new [Github Issue](https://github.com/brutifus/brutifus/issues). This is the best way for everyone to keep track of new problems and past solutions.

## Essential things to know about brutifus

brutifus is primarly a Python package. But in practice, brutifus also includes a series of parameter and
utilitarian files related to its Github repository, and a dedicated documentation hosted using Github pages.

For the sake of clarity, and to facilitate the code maintenance, we list here (succinctly) a series of key facts about the brutifus code, its repository, and the associated organization:

1. **Source code:**
   * brutifus is distributed under the terms of the GNU General Public License v3.0 or later.
   * brutifus adheres to a versioning system specifying the year, month, and release id (during that month),
     i.e. yyyy.mm.v
   * The adopted coding styles are described [here](#styles).
   * brutifus *operational* dependencies are specified in `setup.py`.
   * There is a human-readable [Changelog](CHANGELOG).

2. **Github repository:**
   * Contributions to brutifus get typically merged into the `develop` branch. Pull requests to the
     `master` branch should only originate from the `develop` branch.
   * Any successful pull request to the `master` branch should trigger a new code release.
   * A series of Github Actions are implemented for CI purposes. These include the execution of
     the brutifus tests on Windows, macOS and Linux, a linting of the code, a validation
     of the docs, and a check of the `CHANGELOG`. Upon any push to the `master` branch, the docs
     will also be automatically compiled and published onto the `brutifus.github.io` repository.
   * The `.pylintrc` file refines the behavior of pylint for brutifus.
   * The `master` branch of the `brutifus` repository mirrors the latest release of the code.
     `develop` is the default branch of the repository: this is where the latest not-yet-released
     live.

3. **Documentation:**
   * The brutifus documentation is generated using Sphinx, with the Read-the-docs theme. The compiled
     documentation is hosted on the `master` branch of the `brutifus.github.io` repository of the
     [brutifus Github Organization](https://github.com/brutifus/brutifus.github.io).

## Styles

- **linting:**
  * The following [pylint](https://www.pylint.org/) error codes are forbidden in brutifus: ``E, C0303, C0304, C0112, C0114, C0115, C0116, C0411, W0611, W0612.`` Any pull request will be automatically linted, and these will be flagged accordingly.
  * We encourage contributors to follow PEP8 as closely as possible/reasonable. You should check
    often how well you are doing using the command `pylint some_modified_file.py`.

- **doctrings:**
    Google Style. Please try to stick to the following MWE:
    ```
    """ A brief one-liner description, that finishes with a dot.

    Use some
    multi-line space for
    more detailed info.

    Args:
        x (float | int): variable x could be of 2 types ... note the use of `|` to say that !

           - *float*: x could be a float
           - *int*: x could also be an int

        y (list[str] | str, optional): variable y info

    Returns:
        bool: some grand Truth about the World.

    Raises:
        Exception: if blah and blah occurs.

    Example:
        If needed, you can specify chunks of code using code blocks::

            def some_function():
                print('hurray!')

    Note:
        `Source <https://github.com/sphinx-doc/sphinx/issues/3921>`__
        Please note the double _ _ after the link !

    Caution:
        Something to be careful about.

    """
    ```
    You should of course feel free to use more of the tools offered by [sphinx](https://www.sphinx-doc.org/en/master/), [napoleon](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html), and
    [Google Doc Strings](https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html#example-google). But if you do, **please make sure there are no errors upon generating the docs !**

## Step-by-step guide to contributing

0. Make sure you have git installed. Check that the setup is correct:

       git config --list

   If `user.name` and `user.email` are missing or do not match those of your Github account, [change them](https://docs.github.com/en/free-pro-team@latest/github/setting-up-and-managing-your-github-user-account/setting-your-commit-email-address):

       git config --global user.name "ID+username"
       git config --global user.email "ID+username@users.noreply.github.com"

   **:closed_lock_with_key: Optional but recommended:** use a GPG key to sign your commits. Quoting from [the instructions](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/about-commit-signature-verification): "*Github will verify these signatures so other people will know that your commits come from a trusted source.*" See also [this SO post](https://superuser.com/questions/1512137/which-email-to-sign-commits-with-for-github-and-retain-privacy) and the reply by *Toby* if you use the *@users.noreply...* email for your commits. Do not forget
   to then actually enable git to auto-sign all commits:

       git config --global commit.gpgsign true

0. If you are not a [brutifus collaborator](https://github.com/orgs/brutifus/teams), you will
   have to fork the brutifus repository. If you'd like to become a collaborator instead, you should
   contact the brutifus [dev team](https://github.com/orgs/brutifus/teams).

1. Clone the `develop` branch locally:

        git clone -b develop git@github.com:brutifus/brutifus.git brutifus

   :eyes: If you forked the repo, you should clone your copy of it instead.
   :warning: `develop` is the default branch for brutifus that contains all the latest *approved* changes. **Unless you have a good reason to do otherwise**, this is the branch you want to clone and branch-off from.

2. Actually create your new branch locally:

       cd brutifus
       git checkout -b your_branch_name

3. Check that it all went as expected:

       git branch -a
       git config --list
       git status

4. Install brutifus from this local repo:

       pip install -e ./

   :warning: if you had brutifus already installed, you should remove it first: `pip uninstall brutifus`

5. Modify the code locally. This could be the source code, or the docs `.rst` source files.

   :warning: Please read carefully (and adhere to!) the brutifus [style conventions](#styles) above.

6. Commit changes regularly, trying to bundle them in a meaningful manner.

       git add a_modified_file (OR possibly: git rm a_file_to_delete)
       git commit -m "Some useful, clear, and concise message. Use present tense."

   You can/should also push your branch to the remote repository every now and then:

       git push origin your_branch_name

7. Lint your contributions using the command `pylint some_modified_file.py`. If you want to run the
   checks that will be executed automatically at the pull request stage, you can run the following
   commands from the brutifus repository:

       python ./.github/workflows/pylinter.py --restrict E C0303 C0304 C0112 C0114 C0115 C0116 C0411 W0611 W0612
       python ./.github/workflows/pylinter.py --min_score 8

    Note that this may pick-up linting problems outside of your contribution as well.

8. If warranted, make sure that the docs still compile without errors/warnings:

       cd docs
       sh build_docs.sh

9. Run pytest to check that all is fine with your changes. From the core brutifus folder, type:

       pytest

   :warning: Needless to say, your code tweaks will *evidently* come with dedicated tests. Right ?

10. Once ready with all your modifications, we'll ask that you do a rebase of your branch to
    incorporate any modification that may have occurred on the original `develop` branch in the meantime:

        git fetch origin develop
        git pull --rebase origin develop

    If this leads to conflicts that you are unsure about, get in touch with @brutifus-devs.
    :warning: If you forked the brutifus repository, you'll first need to [sync your fork](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/syncing-a-fork).

11. You can now push your branch to the remote repository. If warranted (it most likely will be!),
    remember to update the `CHANGELOG` and add your name to the `AUTHORS` before doing so:

        git push -f origin your_branch_name

    Note the `-f` flag, required because of the `--rebase` to update the commit history of the
    branch stored on Github.

12. Next, go to `your_branch_name` on the brutifus Github repository, and draft a new pull request. By
    default, the pull request should go from `your_branch_name` to the `develop` branch. Do not forget to link the pull request to a specific issue if warranted. Once the pull request is issued, automated checks will be run (pytest, pylint, changelog, ...). These must all succeed before the changes can be merged (if they do not, there might be something wrong with your changes).

    The @brutifus-devs will then come to formally review the pull request (some reviewers might also be
    set automatically based on the `.github/CODEOWNERS` info).

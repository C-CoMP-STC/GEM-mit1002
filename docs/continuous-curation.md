# Continuous Curation: Best Practices for Model Curation Inspired by Software Development

GEMs are, at their heart, a software product, and we took lessons from software development and applied them to the model curation process. We term this constant testing and iterative model improvement strategy “continuous curation”, inspired by continuous integration/continuous delivery (CI/CD) for traditional software. This included tracking all changes using version control (i.e., Git), having multiple curators collaborate and propose changes by working on branches and opening pull requests (i.e., trunk-based development), testing changes for the model with defined pass/fail tests (i.e., unit tests), and automatically generating artifacts for curator inspection.

None of this is really new, the field is crystallizing

## Motivation
### Why do we need this?

### What have people been doing before this?
#### standard-GEM

### What's new here

## Set-Up
### What is GitHub?
Version control keeps a historical record of changes made to tracked files in a specialized database called a repository (or “repo”).  Git is the software tool that enables version control, and GitHub is one popular cloud-based platform to host Git repositories, that also offers other functionalities such as issue tracking and wiki hosting. While we used GitHub, and will use GitHub-focused terminology (e.g., pull requests, actions) it not the only option for hosting Git repositories, other popular options include GitLab, Bitbucket and Azure DevOps, each of which have analogous tools to those we describe here and could similarly be used for a continuous curation pipeline.

Version control is critical for model curation because it tracks the “who”, “what”, and “when” of all changes made to the model. Who edited the model file, when did they make the edits, and what exactly was changed. It also maintains the historical versions of the model file, so at any time you can revert changes and return to an older version of the model.

Local vs Remote, basic terminology (commit, push, pull, etc.)

![An introduction to basic Git terminology: Local and remote repositories, commit, push, and pull](./figures/png/github-intro.png)

#### Diffs
* How to read a diff

![](./figures/png/diff-sbml.png)

* Not all diffs are equally helpful
    * File types to avoid
        * Microsoft office files: e.g., xlsx

![](./figures/png/diff-excel-vs-csv.png)

For a more in depth coverage, see...

### What's in a Name?: Choosing your Model and Repository Name
* What does BiGG/Palsson do
* What does standard-GEM do

### Branches
Branching is a key feature of Git- it allows developers to isolate their changes so that the main version of the repository is not affected. This allows multiple developers to work simultaneously, and allows developers to test out changes where they will not affect anyone else. We chose to use a branching strategy based on the popular GitFlow strategy. We had two long-lived branches, “main”, the main branch, which had the official releases of the model, and “dev”, the development branch, where all accepted changes to the model were integrate before an official release. All changes made the model were made on feature branches that branched off of and were merged back into the dev branch. This ensured that any new feature development did not disturb the main model. Early on in development, branching is critical to XXX, and later branching became increasingly important to differentiate the version of model from users vs from developers.

![](./figures/png/branches.png)

### Repository Strucure
#### Model
##### What file type to use?
#### Data

### GitHub Actions
We used automation through GitHub actions to run tests and scripts upon the opening of a pull request.
#### Defining an Action with a YML file


## The Continuous Curation Loop

**PUT CONTINUOUS CURATION LOOP HERE**

The Continuous Curation loop consists of 6 steps:
1. Curate
2. Test
3. Report
4. Release
5. Run
6. Monitor

### Step 1) Curate
*NOTE: We do not discuss here how to make curation decisions, but rather how to implement them*
* How big is one curation task?
One critical component of the history of changes to the model is the “why”- why was a change to the model made (e.g., was a reaction found to have genomic evidence, was there a mistake in the biochemistry database, etc.). There are text fields in the model file itself where this information can be stored, and there have been cases in the past of defined “codes” used to represent different types of evidence that support each reaction (CITE EXAMPLES) however we have found that these are not well used, lack standardization across the community, and are often not comprehensive enough to fully explain the reasoning behind each change. We instead elected to documented these in issues and pull requests on the repository. Issues can be used as a sort of electronic lab notebook. To ensure that all curators (present and future) are reminded to document their reasoning, a pull request template was used.
* Open a pull request, that starts the cycle

### Step 2) Test
#### What is a Unit Test?
#### What makes a Good Unit Test?
#### Examples of Unit Tests for Model Curation
Tests were written using the unittest framework and gave a Boolean result- pass or fail. We added additional tests as we went on, and the ones we present here are by no means an exhaustive list of everything that could or should be tested. We tested that the biomass metabolite (cpd11416_c0) added up to 1 g, this is important for dFBA simulations. We tested that there were no erroneous energy generating cycles capable of regenerating ATP without an input carbon source. We checked that there were no dead-end transporters (i.e. external metabolites without an exchange reaction). We checked that the model was not capable of growth without a carbon source in the medium. And that the model recapitulated all known experimental growth phenotypes (depending on exactly how you implement this, it might “fail” for the majority of time of curation. We tested that the SBML file was valid- important as COBRApy may fail to load a model with a malformed fail, and KBase created such files. We tested for isolated genes and metabolites. We tested that all reactions were mass and charge balances. Many of these tests use previously published tools (e.g. MEMOTE), but we found that buy implementing them with unittests on a GitHub action it was easier to track model performance over time and recognize errors introduced into the model quickly.
#### Examples of Unit Tests for Enforcement

### Step 3) Report
#### Scripts vs Tests
#### Examples of Scripts Used for Model Curation
Also ran through GitHub actions, were a set of “scripts”. These differ from tests because they cannot “pass” or “fail”, and instead generate artifacts (e.g., plots) for a human curator to look at. Using the same underlying code as the growth test, we generated a plot of which experimentally known growth phenotypes the model matched or not. The heatmap visualization was useful for sharing results. While this graph is useful to grasp model performance at a glance, it was not the most instructive when gapfilling the model (just knowing that the model does not grow does not help you find a gap). Instead we checked the model’s ability to produce each individual biomass component (i.e., we added a demand/sink reaction for each biomass component that take the metabolite and removes it from the system (similar to an exchange reaction), and looped through the list of biomass components and set each as the objective, maximizing the flux through that sink reaction. A positive value indicated that the model was capable of producing that biomass component. This helped narrow down searches for gaps (e.g., could say that a subset of amino acids was not producible, therefore there must be a gap in that pathway). When using an objective other than the biomass, we considered if there should be free transport/exchange/sinks for all biomass components simultaneously or if only the one being maximized should have a sink. Theoretically, there could be components whose production is tied and without flux through the biomass reaction, dead ends could appear that block flux.
##### Biomass Component Prodcuibility
##### Growth Reports

### Step 4) Release
* What counts as a new version
* Chores upon release
    * different file types
    * MEMOTE
    * MACAW
* GitHub release
* Zenodo release

### Step 5) Run
* This is the fun part, where you, or others, actually try to use the model

### Step 6) Monitor
* Open an issue
* What belongs in an issue

## Why can't you just make it for me?
* Can't you just make an installable tool that makes all of this for me?

## What's the cost?

## How to think like a Programmer
### Questions to Ask Yourself:
1. 

## Glossary
* **Artifact**:
  * *In software engineering*:
  * *In continuous curation*:
* **Branch**: A branch is a parallel version of a repository. It is contained within the repository, but does not affect the primary or main branch allowing you to work freely without disrupting the "live" version. [^gh-glossary].
* **Commit**: A commit, or "revision", is an individual change to a file (or set of files). When you make a commit to save your work, Git creates a unique ID (a.k.a. the "SHA" or "hash") that allows you to keep record of the specific changes committed along with who made them and when. Commits usually contain a commit message which is a brief description of what changes were made. [^gh-glossary].
* **Commit Message**: Short, descriptive text that accompanies a commit and communicates the change the commit is introducing [^gh-glossary].
* **Conflict**: A situation where two branches have changes in a file that Git cannot automatically merge, requiring manual resolution [^harvard].
* **Continuous Integration (CI)**:A development practice where team members frequently integrate their code into a shared repository, often multiple times a day. Each integration is verified by automated builds and tests to detect errors early [^agile].
* **Continuous Delivery**: A software development practice where teams keep their product in a deployable state at all times, but deployment still requires a manual decision [^agile].
* **Continuous Deployment**: A software development practice where every change that passes all automated tests is released to production automatically [^agile].
* **Curation**:
  * *In continuous curation*:
* **Diff**: A diff is the difference in changes between two commits, or saved changes. The diff will visually describe what was added or removed from a file since its last commit [^gh-glossary].
* **Feature Branch**: A branch used to experiment with a new feature or fix an issue that is not in production. Also called a topic branch [^gh-glossary].
* **Git**: Git is an open source program for tracking changes in text files [^gh-glossary].
* **GitFlow**:
* **GitHub**: A web-based platform that facilitates Git's use for collaboration between individuals. Other web-based platforms include GitLab  and BitBucket [^harvard].
* **GitHub Actions**:
* **Issue**: Issues are suggested improvements, tasks or questions related to the repository. Issues can be created by anyone (for public repositories), and are moderated by repository collaborators. Each issue contains its own discussion thread. You can also categorize an issue with labels and assign it to someone [^gh-glossary].
* **JSON**:
* **Local**: 
* **Merge**: Merging takes the changes from one branch (in the same repository or from a fork), and applies them into another. This often happens as a "pull request" (which can be thought of as a request to merge), or via the command line [^harvard].
* **Monitor**:
  * *In software engineering*:
  * *In continuous curation*:
* **Pull**: The process of integrating changes from one version of a repository to another (e.g. from a fork back to the original repo, or from a branch back to the main branch). There are two general use cases: 1) When the owner of a repository makes changes to it, you pull those changes into your local copy. 2) When you make changes to a forked repository or a branch of a repository and want to incorporate the changes back to the original repo or branch, you initiate a pull request, and then whoever is in charge of the original repository can pull those changes in  [^harvard].
* **Pull Request**: Pull requests are proposed changes to a repository submitted by a user and accepted or rejected by a repository's collaborators [^gh-glossary].
* **Push**: To push means to send your committed changes to a remote repository on GitHub.com. For instance, if you change something locally, you can push those changes so that others may access them [^gh-glossary].
* **Release**:
  * *In software engineering*:
  * *In continuous curation*:
* **Remote Repository**: This is the version of a repository or branch that is hosted on a server, most likely GitHub.com [^gh-glossary].
* **Report**:
  * *In software engineering*:
  * *In continuous curation*:
* **Repository/Repo**: A repository is the most basic element of GitHub. They're easiest to imagine as a project's folder. A repository contains all of the project files (including documentation), and stores each file's revision history. Repositories can have multiple collaborators and can be either public or private [^gh-glossary].
* **Run**:
  * *In software engineering*:
  * *In continuous curation*:
* **SBML**:
* **Script**:
  * *In software engineering*:
  * *In continuous curation*:
* **Semantic Versioning**:
* **Test**:
  * *In software engineering*:
  * *In continuous curation*:
* **Trunk-Based Development**:
* **Version Control**: The process of tracking changes to files over time, allowing you to recall specific versions later [^harvard].
* **XML**:

[^harvard]: https://informatics.fas.harvard.edu/resources/glossary/#git-terms
[^agile]: https://www.agile-academy.com/en/agile-dictionary/
[^gh-glossary]: https://docs.github.com/en/get-started/learning-about-github/github-glossary

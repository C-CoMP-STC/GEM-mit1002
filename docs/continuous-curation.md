# Continuous Curation: Best Practices for Model Curation Inspired by Software Development

GEMs are, at their heart, a software product, and we took lessons from software development and applied them to the model curation process. We term this constant testing and iterative model improvement strategy “continuous curation”, inspired by continuous integration/continuous delivery (CI/CD) for traditional software. This included tracking all changes using version control (i.e., Git), having multiple curators collaborate and propose changes by working on branches and opening pull requests (i.e., trunk-based development), testing changes for the model with defined pass/fail tests (i.e., unit tests), and automatically generating artifacts for curator inspection.

## Motivation
### Why do we need this?

### What have people been doing before this?

### What's new here

## Set-Up
### What is GitHub?
Version control keeps a historical record of changes made to tracked files in a specialized database called a repository (or “repo”).  Git is the software tool that enables version control, and GitHub is one popular cloud-based platform to host Git repositories, that also offers other functionalities such as issue tracking and wiki hosting. While we used GitHub, and will use GitHub-focused terminology (e.g., pull requests, actions) it not the only option for hosting Git repositories, other popular options include GitLab, Bitbucket and Azure DevOps, each of which have analogous tools to those we describe here and could similarly be used for a continuous curation pipeline. Version control is critical for model curation because it tracks the “who”, “what”, and “when” of all changes made to the model. Who edited the model file, when did they make the edits, and what exactly was changed. It also maintains the historical versions of the model file, so at any time you can revert changes and return to an older version of the model.

GitHub is an implementation of Git- and so much more. A repositry also incldued issues, GitHub Actions, pull requests

#### Diffs
* How to read a diff
* Not all diffs are equally helpful
    * File types to avoid
        * Microsoft office files: e.g., xlsx

For a more in depth coverage, see...

### What's in a Name?: Choosing your Model and Repository Name

### Branches
Branching is a key feature of Git- it allows developers to isolate their changes so that the main version of the repository is not affected. This allows multiple developers to work simultaneously, and allows developers to test out changes where they will not affect anyone else. We chose to use a branching strategy based on the popular GitFlow strategy. We had two long-lived branches, “main”, the main branch, which had the official releases of the model, and “dev”, the development branch, where all accepted changes to the model were integrate before an official release. All changes made the model were made on feature branches that branched off of and were merged back into the dev branch. This ensured that any new feature development did not disturb the main model. Early on in development, branching is critical to XXX, and later branching became increasingly important to differentiate the version of model from users vs from developers.

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
  * *In biology*:
  * *In continuous curation*:
* **Branch**: A copy of a repository from a certain point within that repository's history. Typically a repository has a "main" branch and other branches are created off of it. Changes on the main branch are not reflected in the split branch unless explicitly synced and vice versa  [^1].
* **Commit**: The process of saving changes to the repository. This is done after adding files to the staging area  [^1].
* **Conflict**: A situation where two branches have changes in a file that Git cannot automatically merge, requiring manual resolution  [^1].
* **Continuous Integration**:
* **Continuous Delivery**:
* **Continuous Deployment**:
* **Curation**:
  * *In ???*:
  * *In continuous curation*:
* **Diff**:
* **Fetch**:
* **Git**: Software for version control, which keeps track of changes to files in a given directory [^1].
* **GitFlow**:
* **GitHub**: A web-based platform that facilitates Git's use for collaboration between individuals. Other web-based platforms include GitLab  and BitBucket [^1].
* **GitHub Actions**:
* **Issue**:
* **JSON**:
* **Local**: 
* **Merge**: The process of combining changes from one branch into another branch, typically done as part of a pull request [^1].
* **Monitor**:
  * *In software engineering*:
  * *In continuous curation*:
* **Pull**: The process of integrating changes from one version of a repository to another (e.g. from a fork back to the original repo, or from a branch back to the main branch). There are two general use cases: 1) When the owner of a repository makes changes to it, you pull those changes into your local copy. 2) When you make changes to a forked repository or a branch of a repository and want to incorporate the changes back to the original repo or branch, you initiate a pull request, and then whoever is in charge of the original repository can pull those changes in  [^1].
* **Pull Request**: When someone has made changes to a fork or a branch that they wish the owner's or the original repository to incorporate, they initiate a pull request so the owner can review and potentially pull the changes [^1].
* **Push**: The process of uploading committed changes from a local repository to a remote repository to a remote platform (e.g. Github) [^1].
* **Release**:
  * *In software engineering*:
  * *In continuous curation*:
* **Remote Repository**: A repository that is hosted on a server, typically on a web-based platform like Github [^1].
* **Report**:
  * *In software engineering*:
  * *In continuous curation*:
* **Repository/Repo**: A directory of files that has been initialized by Git for syncing, possibly including code, documentation, or data [^1].
* **Run**:
  * *In software engineering*:
  * *In continuous curation*:
* **SBML**:
* **Script**:
  * *In software engineering*:
  * *In continuous curation*:
* **Test**:
  * *In software engineering*:
  * *In continuous curation*:
* **Trunk-Based Development**:
* **Version Control**: The process of tracking changes to files over time, allowing you to recall specific versions later [^1].
* **XML**:

[^1]: https://informatics.fas.harvard.edu/resources/glossary/#git-terms

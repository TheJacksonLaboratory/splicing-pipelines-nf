# GitHub basic commands

The following is a list of useful GitHub commands. If your still confused after reading it don't hesitate to ping [@PhilPalmer](https://github.com/PhilPalmer) or [@adeslatt](https://github.com/adeslatt). You can also read this [this guide](https://guides.github.com/introduction/git-handbook/) (10mins read) for a more comprehensive introduction to GitHub.

## Download the repository from GitHub.com to your machine

As the repository is currently private you wil need to sign in to your GitHub account. (You can also [chache your password or generate SSH keys](https://help.github.com/en/github/getting-started-with-github/set-up-git#next-steps-authenticating-with-github-from-git) to prevent needing to enter your username and password everytime you push or pull from GitHub).

```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git
```

## Change into the \`repo\` directory

```bash
cd splicing-pipelines-nf
```
 
## To ensure you are up-to-date with the `master` branch or whatever branch you may be on

This command updates the local branch (line of development) with updates from its remote counterpart. Developers use this command if a teammate has made commits to a branch on a remote, and they would like to reflect those changes in their local environment.

```bash
git pull
```

## To checkout on your own branch 

Remember to not stay on a branch long, so your changes do not get too far out of sync with the master. A common convention for branch names is `yourname-feature` eg `phil-docs-update`

```bash
git checkout -b [name of your branch]
```

Once you are on your own branch you can then make changes, for example, edit `README.md` using the text editor

## To show the status of changes

Changes are shown as either untracked, modified, or staged. It will also show you what branch you're on

```bash
git status
```

## To add, commit and push your changes.

### add
Stage the changed files eg `git add README.md`
```bash
git add [whatever you have added]
```

### commit
This is like taking a snapshot of all of the files that have been added
```bash
git commit -m “a useful message regarding what you have changed”
```
 
### push
Push changes to github
```bash
git push
```

![git_cmds](https://raw.githubusercontent.com/lifebit-ai/images/master/jax_splicing/git_cmds.gif)

## To push your local branch onto master

```bash
git push https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git [name of your branch]
```
 
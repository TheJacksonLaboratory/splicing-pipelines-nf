# github basic commands

## git clone your repository

```bash
git clone https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git
```
 
## To ensure you are up-to-date with the `master` branch or whatever branch you may be on:

```bash
git pull
```

## To add, commit and push your changes.

### add
```bash
git add [whatever you have added]
```

### commit

```bash
git commit -m “a useful message regarding what you have changed”
```
 
### push

```bash
git push
```

## To checkout on your own branch 

Remember to not stay on a branch long, so your changes do not get too far out of sync with the master:

```bash
git checkout -b [name of your branch]
```

## To push your local branch onto master


```bash
git push https://github.com/TheJacksonLaboratory/splicing-pipelines-nf.git [name of your branch]
```
 
# PHLAWD

There are several versions of phlawd floating around:

* The official version of phlawd at https://github.com/blackrim/phlawd

* The "chinchliff" fork of phlawd was last updated in Mar 2016 and contains a number of bug fixes and features, including the ability to use whole genome sequences and "shrinking" to the target gene.

This is an unofficial fork, which updates the code for the GenBank ID change (https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/) and downloads GenBank daily updates in addition to the regular releases.

## Requirements

phlawd requires the GCC compiler, as well as `wget`, `mafft`, `muscle`, `quicktree`, and `sqlite` to be in the user's path. These can be installed via Homebrew (https://brew.sh) or your system package manager.

## Installing

### via Homebrew

Install Homebrew by following the instructions at <https://brew.sh>, then run:

```
brew install brewsci/bio/phlawd
```

### Building from source

```bash
git clone https://github.com/jonchang/phlawd.git
cd phlawd/
cd src/
./configure
make
install PHLAWD /usr/local/bin # or wherever...
```

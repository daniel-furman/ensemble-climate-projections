### Requirements

---

Python dependencies are listed in a `requirements-py.txt` file, including the library version numbers. You can replicate the environment your codebase needs by using virtualenv:

```
# This creates the virtual environment
cd $PROJECT-PATH
virtualenv ensemble-climate-projections

# Then install the dependencies by referring to the requirements-py.txt:

# This installs the modules
pip install -r requirements-py.txt

# This activates the virtual environment
source ensemble-climate-projections/bin/activate
```
R dependencies are listed in a `requirements-R.txt` file. You can replicate the environment your codebase needs by using install.packages():

```
reqs <- read.table('PROJECT-PATH/requirements-R.txt', col.names = c('package', 'version'))
for (i in 1:length(reqs$package)){
  install.packages(reqs$package[i]) # devtools::install_version if versions desired
  }
```

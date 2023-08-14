# Regression Tests

The regression tests are organized to the `pytest` framework. Please ensure you 
have the `pytest` package installed in your Python environment before doing 
regression tests. More information for `pytest` can be found 
[here](https://docs.pytest.org/en/latest/getting-started.html)

### Run test

1. Download the data files that are needed for the regression tests

```
wget --no-check-certificate https://software.rcc.uchicago.edu/mscg/downloads/data.tar.gz
tar xzvf data.tar.gz
mv data ${path_to_openmscg_package}/tests/
```


2. Go into the root folder of the OpenMSCG package (where the file`pytest.ini`
exists)

3. Run `pytest`

```
pytest
```


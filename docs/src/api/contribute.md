# Contributing

While not strictly adhering to the [SciML Style Guide](https://github.com/SciML/SciMLStyle), we recommend consulting the guide as a baseline for contributions to this package. Refer to the [SciML Contributors Guide](https://github.com/SciML/ColPrac/blob/master/README.md) as well. Ideally, your contribution features:

 - tests for new or changed functionality
 - docstrings for relevant functions
 - documentation and examples

We would also love to feature your work with this package as part of the **Examples** section.

## Running the tests

To run the full test suite, use `] test` in the package environment. To run a single test
file standalone, run `julia --project=test -e 'using Pkg; Pkg.instantiate()'` once, then
`julia --project=test test/<path>.jl`. Tests live in a Pkg workspace subproject
(`test/Project.toml`); new test-only dependencies are added there, not to the root
`Project.toml`.

# Contributors

The following is a list of people who have made significant contributions to the development of BMO:

- Hugo Uittenbosch
- Oliver Kliebisch
- Aurelius Manny
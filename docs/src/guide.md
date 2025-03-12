# Developers guide

## Contributing

While not strictly adhering to the [SciML Style Guide](https://github.com/SciML/SciMLStyle), we recommend consulting the guide as a baseline for contributions to this package. Refer to the [SciML Contributors Guide](https://github.com/SciML/ColPrac/blob/master/README.md) as well. Ideally, your contribution features:

 - tests for new or changed functionality
 - docstrings for relevant functions
 - documentation and examples

We would also love to feature your work with this package as part of the **Examples** section.

## Documentation development

If you want to edit the package documentation locally, follow these steps:

1. Create your local dev. repository via `] dev BeamletOptics`
2. Switch into the `docs` environment, e.g. `] activate .` inside of the `docs` folder
    1. Inside of [VS Code](https://code.visualstudio.com/) you can activate the local environment by right-clicking the `make.jl` file
    2. If you have the Julia plugin installed, you will be able to select `Julia: Activate This Environment`
3. Inside of the `docs` environment switch the dependency onto your local `BeamletOptics` dev folder via `] dev BeamletOptics`
    1. This step is **important**, otherwise an incompatible version of `BeamletOptics` might be used to generate the docs
4. Run the `make.jl` file

Changes you have made will then be saved into the `build` folder. You can host the website locally by opening the `index.html` starting page.
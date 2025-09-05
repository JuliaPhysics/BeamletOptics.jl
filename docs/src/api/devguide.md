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

### Section titles

When creating a custom section in the documentation, you should avoid naming the section the same way as your type, e.g. for `MyCustomType` you should not create a section that is called `# MyCustomType`. The reason for this is that the `@ref` macro will confuse the docstring of your type with the section header, leading to undefined behavior for any links pointing to the embedded docstring via `[`MyCustomType`](@ref)`.

### Creating figures

In general, you can generate and include figures into your documentation section any way you see fit. We strongly urge you to use the existing `CairoMakie` or `GLMakie` backend. However, with the increasing amount of plots and corresponding scripts the build time for the docs in a local environment has become unsustainable. Therefore, for the BMO docs we recommend that you adhere to the following design pattern:

1. Place your code in a standalone `.jl`-script within the `docs\src\assets` folder. 
    - refer to existing scripts for a rough guideline
    - do not forget to activate the appropriate backend via e.g. `GLMakie.activate!()`
2. Make sure that each figure in your script is saved via a `save("my_fig.png", my_fig, ...)` statement
    - files will be saved with respect to the calling environment
3. In the markdown file that contains your documentation and should load your images, do the following:
    1. create a `@setup` code block
    2. include the `conditional_include` via e.g. `include(joinpath(@__DIR__, "..", "assets", "cond_save.jl"))`
        - ensure that the relative file path is correct
    3. load your script during build via `conditional_include`
        - more info on this function is provided in the `cond_save.jl` file
4. Load the image within your markdown file via `![My figure](my_fig.png)`

Examples for this pattern can be found at the top of most .md files of the documentation, e.g. `beamsplitters.md`.
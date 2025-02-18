# Developers guide

## Documentation development

If you want to edit the package documentation locally, follow these steps:

1. Go to your local repository via `] dev SCDI`
2. Switch into the `docs` environment, e.g. `] activate .` inside of the `docs` folder
    1. Inside of [VS Code](https://code.visualstudio.com/) you can activate the local environment by right-clicking the `make.jl` file
    2. If you have the Julia plugin installed, you will be able to select `Julia: Activate This Environment`
3. Inside of the `docs` environment switch the dependency onto your local `SCDI` dev folder via `] dev SCDI`
    1. This step is **important**, otherwise an incompatible version of `SCDI` might be used to generate the docs
4. Run the `make.jl` file

Changes you have made will then be saved into the `build` folder. You can host the website locally by opening the `index.html` starting page.
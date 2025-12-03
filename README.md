# Ponca
Ponca is a header only C++/CUDA library for point cloud analysis and surface reconstruction.

[![codecov](https://codecov.io/github/poncateam/ponca/branch/master/graph/badge.svg?token=NWSHQWK6NO)](https://codecov.io/github/poncateam/ponca)

## Development status
Ponca is currently under refactoring as we are working on the next release v2.0, which introduces breaking changes in the API. We do not change everything, the most important concepts will remain, but the library will offer more versatile code, will be easier to understand, and its structure will be more consistent with the involved mathematical models.

Hence, depending on you needs, we recommend to:

 - Stick to tag `v1.4`: this is the last release that follows to the current API. Mandatory if you have production code already working and are not rushing for updated. Be aware that bug fixes are planned to be introduced directly to v2.0, so an update will be necessary in the end.
```bash
git checkout v1.4
```
 - Use tag`v2.0.alphaX`: the latest release introducing the new changes in the API. Each alpha release might introduce new breaking change overtime depending on our progresses. Each alpha release is however tested and its behavior should be reliable: the `alpha` state is related more to instabilities in the API, than in the data-structures or algorithms. Recommended for newcomers, or for developers who want to anticipate the upcoming API changes.
```bash
git checkout v2.0.alphaX
```
 - Use head of `master` branch, to monitor our progresses or to contribute.
```bash
git checkout master
```


## Documentation
The official documentation is available at: https://poncateam.github.io/ponca/index.html, and includes:
 - [Starting guide](https://poncateam.github.io/ponca/ponca_getting_started_page.html)
 - [User manual](https://poncateam.github.io/ponca/user_manual_page.html)
 - [Examples](https://poncateam.github.io/ponca/example_page.html)
 - [C++ programming reference](https://poncateam.github.io/ponca/annotated.html)


Ponca is developed by the STORM Research group (https://www.irit.fr/STORM/site/) at the Institut de Recherche en Informatique de Toulouse (https://www.irit.fr).
Prior development on Patate (previous fancy name) have been carried out at Inria Bordeaux (Manao Research group).


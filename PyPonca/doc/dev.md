# Ponca C++ Binding Dev doc

## API Design choices

### Core API design resonning

The Python API from the user perspective is much different from the C++ side. The typical workflow in C++ is:

```c++
Fit f;
f.setNeighborFilter({loc, radius});
f.compute(pointcloud);

// Do whatever with f, knowing its capabilities
```

In python, the API looks like:

```python
f = Fit()
f.setNeighborFilter(loc, radius) # Or an array of such

# Describe what to do with f once fitted
f.addComputation(id, params)
f.addComputation(id, params)

# Get results of computation
result = f.compute(pointcloud)
```

This API arises from the fact that we don’t want the user to iterate over the point cloud to keep as much performance as possible. For this reason, the user may supply the list of filter locations and radii as a whole instead of iterating over them one by one and performing the fit, which would be both slow and cumbersome to parallelize in Python.

Copying the C++ pattern in this case is possible (and was tried), but two problems arise. From a performance perspective, the best solution would be to keep a copy of all the fitted objects and apply the result extraction separately to each. This has the disadvantage of implying a high memory cost and a loss of parallelization advantages, as the mechanism of creating threads (or launching CUDA kernels) can be costly in terms of computational overhead. On the other hand, not keeping every object in memory destroys performance, as each call to an extraction function would require re-fitting the object.

The choice was made for the user to specify the list of computations before the call to compute. This allows to reconcile performances in both time and memory but at the expense of lacks of user controls. The user must remain within the bounds of the predefined operations. 


## Name mangling

Ponca C++ uses templates to provide generic yet performant code. However, templates are compile-time types and are therefore not directly compatible with Python. Instead, the bindings must instantiate and expose every template specialization that should be available to the user.

In Ponca, this is particularly challenging because the core type, `Point`, is itself a template parameter. Consequently, every class templated on `Point` would need to be bound under a different name. From the user's perspective, this would require using different classes depending on the scalar type or the dimension of the points, which is inconvenient.

To address this, we instead rely on name mangling and a Python-generated dispatcher. Within the C++ bindings, each class is registered under a unique mangled name. This name is computed from both the algorithm name (e.g., `APSS`) and the template type information (e.g., `Point<float, 3>` becomes `3f`).

In `__init__.py`, we define a Python class that, once the required type information is known (i.e., when `compute()` is called), computes the mangled name of the corresponding C++ class and forwards all method calls to the appropriate implementation.


### Name mangling rules

We list here the rules for mangling. 

#### Array mangling rules

We encourage using the `MangleArray` function on the C++ side and `_pyponca.internal.mangleArray` on the Python side to ensure consistency between both implementations.

The first dimension of an array is never encoded in the mangled type. Subsequent dimensions are encoded as a list of their sizes, separated by underscores (`_`). For example, an array of shape `(12, 8, 4)` is encoded as `"8_4"`.

The array's scalar type is encoded according to the following mapping (using the first letter of the corresponding C++ type whenever possible):

* `float` → `"f"`
* `double` → `"d"`
* any other type → `"unknown"`

The shape and type encodings are then concatenated directly. For example, a `float` point cloud of shape `(N, 3)` has the mangled name `"3f"`.

#### Filter mangling

Filter mangling is specified manually within Fitting/Filters.h. It commonly uses a shorthand in order for names to remain short. 

#### Compute Object Mangling

Compute object names are constructed by concatenating the method name, the array mangling, the point type mangling then the filter one. For instance:

* "MongePatchRestrictedQuadratic3fPNCW" is built by the concatenation of "MongePatchRestrictedQuadratic" + "3f" + "PN" + "CW", hence is the specialization of MongePatchRestrictedQuadratic method for 3d float point cloud with normals and the constant weight filter. 

## Test suites


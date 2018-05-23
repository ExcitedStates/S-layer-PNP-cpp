# S-layer-PNP-cpp
multiscale electro-diffusion reaction transport model

## How to compile?

### Mac

```
clang++ -Ofast -o pnp3d cpp_pnp3d.cpp
```

### Linux (with Intel compilers)

```
icc -Ofast -o pnp3d cpp_pnp3d.cpp
```

### Matlab (MEX)

```
mex CXXFLAGS='-Ofast' mex_pnp3d.cpp
```

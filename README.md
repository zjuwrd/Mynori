# Nori with debiasing methods

Based on `Nori`, this repository is aimed for the implementation of unbiased photon mapping and unbiased photon mapping.

## 1. Building
Before building the project, please use the following command to initialize the submodules.
```bash
git submodule update --init --recursive
```

The project is constructed and tested in `WSL2-ubuntu22.04` and adopts `cmake` to manage project building. To build the project in linux you can use the following commands.
```bash
mkdir build && cd build
cmake ..
make
```
To build the project in `windows`, you will need `cmake` together with `Visual Studio 2022`.


## 2. testing
The codes for photon-mapping-related algorithm locate in `./src/integrator/`. 

The scene resources for debiasing-method tests locate in scenes/extra/water-caustic. To test debiased photon mapping, run the following command:

```bash
./build/nori ./scenes/extra/water-caustic-XXX.xml
```
where `XXX` is replaced with specific rendering method. `debiasedPM` for debiased photon mapping, `debiasedPPM` for debiased progressive photon mapping, and `SPPM` for stochastic photon mapping.

Parameters for rendering integration are passed by the `xml` description files. For example, the parameters passed towards debiased PM are as follows:
``` xml
<integrator type = "debiasedpm">
    <float name="radius" value="0.04"/>
    <integer name="photonsunit" value="131072" />
    <integer name="k" value="10"/>
    <integer name="iterations" value="200" />
    <integer name="record" value="1" />
</integrator>	
```
+ `radius`: initial radius 
+ `phtonsunit`: number of photons to be collected per iteration
+  `k` : the initial k for the discrete distribution $p(j)$,
+ `iterations` : total iterations to make
+ `record`: determines whether to record results after each iteration

You can modify them to see what difference would be made by the changes.

To create new scenes and test them, please look into `Nori.md` for further `xml` file format constraints.

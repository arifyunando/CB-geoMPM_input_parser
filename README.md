# ASCII Mesh and Particles File Generator for CB-GEO MPM
> [CB-Geo Computational Geomechanics Research Group](https://www.cb-geo.com)

## Documentation
Main code are all written in mp_parser.py. To use it, import the module with Jupyter Notebook and call the class. File will be automatically generated in the result folder

> examples
- Generating Mesh
```python
import mp_parser

#Mesh
a = mp_parser.mesh(1)   # Instantiate Object
a.setMesh(20, 20)       # Set Mesh Size
a.printFile(filename)   # print file
```

- Generating Particles
```python
import mp_parser
import matplotlib.pyplot as plt

points = [
    (0,0),
    (5,0),
    (7,3),
    (4,5),
    (0,5)
]

# Generating Files
b = mp_parser.particle(points, delta = 0.1)
b.printFile(digits=4)

# Optionally, you can also preview the shape with Matplotlib
filled = b.nodes
x = [i[0] for i in filled]
y = [i[1] for i in filled]
plt.plot(x,y, '.')
plt.axis('scaled')
plt.show()
```

For CB-Geo Documentations, please refer to [CB-Geo MPM Documentation](https://mpm.cb-geo.com/). They also include more spesific ASCII Criterion for the input files. 

## Author
* Arif Yunando - [more information..](https://arifyunando.github.io)

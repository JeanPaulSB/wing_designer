# Wing Designer ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)

Wing Designer is a project that help you to analyze the aerodynamical performance of bodies, such as wings, airfoils and so on.
It is built in the top of Python. The project is still in development, If you'd like to contribute, let me a message
through my email.

# Features
  - Airfoils 
    - Simple Airfoil
      - Naca(4,5) digits
      - Custom airfoil
    - Flapped Airfoil (plain flap)
      - Naca(4,5) digits
      - Custom airfoil
  - Wings
    - any geometry
  - Physical Models
    - Thin Airfoil Theory
    - Lifting Line Theory

## Getting Started

### Overall architecture
The tools inside the project are organized as python packages, the availables until now are
- Geometry (pretty self explanatory, it contains python object definitions for objects such as airfoils, wings, and surfaces)
- Models (it packages all the solvers that implements the mathematical and physical background required to estimate the aerodynamical performance of the specified geometry)


Geometry currently implements definitions for:
- Airfoil
- Custom Airfoil
- Naca 

Models currently implemented:
- Thin Airfoil Theory:
  Mathematical model that replaces the mean camberline of a given airfoil for a vortex sheet, under that abstraction it is possible to estimate the aerodynamical perfomance of an airfoil assuming that the flow is inviscid (still being a reasonable aproximation and it states good results)
- Lifting Line Theory:
This models deals with the wing finite representation in contrast to the Thin Airfoil Theory, with treats the Airfoil as an infinite wing, it is important to analyze the differences between those two behaviors:
A finite wing inmersed in a flow vs An infinite wing in the some condition. It uses a collection of horseshoe vortex along the span of the wing and it is able to compute the aerodynamical performance of an actual wing with the considerations of its vorticity and the corresponding induced drag due to that.
> [!NOTE]  
> Infinite wing actually means that the span of the airfoil is 1.

> [!NOTE]  
> In the future, more geometry definitions and models will be added.


### How to run it

The workaround is pretty straigthforward. The operations are stated to be intuitive and setted in a high level so anyone regardless its coding experience
could use it to work with it.

Reading classes and docs definition may ligth your perspective of the scope of each object,geometry or model.

### Case 1: Let's suppose that you'd like to estimate the Cl for a NACA 2412 airfoil.
The strategy would be as simple as, create a Naca object, passing parameters such as NACA digits, chord and number of points for its definition.
Then, you just have to call the model ThinAirfoilTheory and apply it to your airfoil.
What matters to this result is the aL_0 of your airfoil, since TaT predicts that regardless the airfoil shape, the lift slopte would be the same.
Now you can just extract the aL_0 from your airfoil and use it considering that $Cl = {2\pi *(alpha-aL_0)}.

Let's traduce that to code.

```python
from geometry import Naca
from models import ThinAirfoilTheory
import numpy as np

airfoil = Naca("2412",15,50) # NACA 2412, 15m of chord and 50 points for its definition.
ThinAirfoilTheory.solve(airfoil)
# cl @ 5deg
cl = 2*np.pi*(np.deg2rad(5) - airfoil.aL_0)

```


### Sample outputs

<img src="https://i.imgur.com/cuNHSXT.png" width="500"><br>

<img src="https://i.imgur.com/eBbF9up.png" width="500"><br>
<img src="https://i.imgur.com/khXyXJK.png" width="200">

### Prerequisites

![Static Badge](https://img.shields.io/badge/Python-blue?logo=python&logoColor=white&link=https%3A%2F%2Fwww.python.org%2F)

install all the required libraries by using

```

pip install -r requirements.txt

```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```

Give the example

```

And repeat

```

until finished

```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```

Give an example

```

### And coding style tests

Explain what these tests test and why

```

Give an example

```

## Built With

- [Python](http://www.dropwizard.io/1.0.2/docs/)

## What's next for the project

- Implementing a GUI.
- Making the project capable of computing lift,drag and moment coefficient using the correspoding theoretical model in Aerodynamics.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Author

- **Jean Paul Sierra**

```

```

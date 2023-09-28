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
- Geometry (pretty self explanatory, it contains python object definitions for object such as airfoils, wings, and surfaces)
- Models (it packages all the solvers that implements the mathematical and physical background required to estimate the aerodynamical performance of the specified geometry)


Geometry currently implements definitions for:
- Airfoil
- Custom Airfoil
- Naca 
> [!NOTE]  
> In the future, more geometry definitions will be added.


### How to run it

```
As the software architecture suggests, start by importing the required tools from its corresponding folder, i.e
if you'd like to i
```

**a few more methods will be added in future versions.**

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

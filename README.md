# BAC2

Binding Affinity Calculator 2, a ground up rewrite of the BAC tool intended to unify and 
simplify the maintenance and use of the tool, and provide a platform on 
which new functionality can be built.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine 
for development and testing purposes.

### Prerequisites

You will need Python 3.6 or higher. Using a conda environment is the easiest way to go.


### Installing

You need to clone this repository and run 

```
pip install BAC2/
```

## Authors

* **Kristof Farkas-Pall**
* **Dave W. Wright**
* **Stefan Zasada**
* **Fouad Husseini**

See also the list of [contributors](https://github.com/UCL-CCS/BAC2/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details


## Structure of this project

- bac - code
    - build - build layer API top level directory
    - present - present layer API top level directory
    - simulate - simulate layer API top level directory
    - remote - execution of simulation on remote servers
- client - client code
- docs - documentation

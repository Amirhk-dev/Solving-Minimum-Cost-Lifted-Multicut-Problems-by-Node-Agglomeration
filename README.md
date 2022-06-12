# This is a repository related to the "Solving-Minimum-Cost-Lifted-Multicut-Problems-by-Node-Agglomeration" paper

[Project Page](https://web.informatik.uni-mannheim.de/akardoos/)

## In this paper two heuristic solvers are proposed for the minimum cost (lifted) multicut problem:
* Balanced Edge Contraction (BEC): Creates balanced clusters
* BEC-cut: Creates clusters along the object boundaries

## Setup
* The code is based on the [graph and graph algorithms library](https://github.com/bjoern-andres/graph)

* To create the instance of the Lifted Multicut Problem (LMP) we use [this folder.](https://www.mpi-inf.mpg.de/fileadmin/inf/d2/levinkov/iccv-2015/code.tar.gz)

* After downloading the folder add the two solvers to the directory:
"code\include\andres\graph\multicut-lifted\"

* Compile the library

## References
````
@inproceedings{Kardoost2018,
title     = {Solving Minimum Cost Lifted Multicut Problems by Node Agglomeration},
author    = {A Kardoost and M. Keuper},
year      = {2018},
booktitle = {ACCV 2018, 14th Asian Conference on Computer Vision}
}
````
````
@article{Keuper2015,
title    = {Efficient Decomposition of Image and Mesh Graphs by Lifted Multicuts},
author   = {M. Keuper and E. Levinkov and N. Bonneel and G. Lavou{\'e} and T. Brox and B. Andres},
journal  = {ICCV},
year     = {2015}
}
````

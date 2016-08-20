Algorithm
-------
group_lapl; sparse group Laplacian shrinkage.


Maintainer
-------
Jin Liu   <jin.liu@duke-nus.edu.sg>


Publication
-------
Liu, J., Huang, J., & Ma, S. (2013). Incorporating network structure in integrative analysis of cancer prognosis data. Genetic epidemiology, 37(2), 173-183.


Description
-------
SGLS implements penalization method for integrative analysis of multiple high-throughput cancer prognosis studies incorporating network structures.  This method is based on a combination of the group MCP penalty and a Laplacian penalty. The group MCP is adopted for gene selection and Laplacian penalty is applied to smooth the differences between regression coefficients of tightly-connected genes.


Usage
-------
source("sp.r")
- gmcp_lapl() solves for a fixed tuning.
- sp() produces the solution path.

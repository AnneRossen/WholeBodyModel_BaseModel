# WholeBodyModel_BaseModel
Whole Body Model recreation of “A whole-body multi-scale mathematical model for dynamic simulation of the metabolism in man", Carstensen et al. 2022

## Mathematical approach
The general model is described as a system of metabolites, reactions and flow in and out of organs. 
The general dynamics of a single compartment is defined by:

$$V \frac{dC}{dt} = M (Q_{in}C_{in} - Q_{out}C) + RV , $$

where $V$ is the volume, $C$ is a vector containing the concentration of the metabolites, $M$ is the external and internal component ordering, 
$Q_{in}$ is the flow rate of what goes in, $C_{in}$ is a vector containing the concentration of the metabolites that flow in,
$Q_{out}$ is the flow rate of what goes out and $R$ is the production rates. 
The compartments are coupled through concentration gradients in the blood vessels that connect the compartments. 
$M$ is a square matrix containing only ones and zeros in the diagonal corresponding to the metabolites distributed through the blood vessels (circulating metabolites). 
For instance, the circulating metabolite, $C_i$, corresponds to $M_{i,i}=1$. The production rate $R$ is incorporated as a vector defined by
$$R = (T S)' T r $$
where $T$ is a matrix of reactions that occur, $S$ is a stoichiometric matrix containing all reactions and $r$ is a vector with the kinetics for the reactions. 
$T$ contains ones and zeros corresponding to which reactions from the stoichiometric matrix, that are present in the compartment. 
For instance, a compartment which involves reaction $1,3$ and $5$ from the stoichiometric matrix have $T_{1,1} = 1, T_{2,3} = 1, T_{3,5} = 1$, and zeros elsewhere. 
The reaction rate vector, $r$, is a function of the concentration of each metabolite:
$$r = r(C) $$
To utilize this, in a whole-body model, it must be formulated for each compartment.

The model also consist of a SIMO sub-model that accounts for metabolism of the three macronutrients, carbohydrates, proteins and lipids and a Hormonal submodel is also implemented and accounts for insulin and glucagon effects.

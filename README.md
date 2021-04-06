# Robust Dual Dynamic Programming (RDDP)
## Usage
```julia
pkg> add RDDP
```
## Example
```julia
# inventory management under uncertain demand
using RDDP,Gurobi,JuMP
N = 10
msro = RDDP.buildMultiStageRobustModel(
    N_stage = N, # number of stages
    optimizer = Gurobi.Optimizer, # choose an optimizer
    MaxIteration = 100, # maximum epoches
    MaxTime = 60, # maximum trainning time
    Gap = 0.01, # optimal gap
    use_maxmin_solver = false # using a vertex enumeration by default
) do ro::JuMP.Model,t::Int
    # defining operator's strategy
    @variable(ro,0 <= x <= 5>,RDDP.State,initial_value = 0) # state variable (stock level)
    @variable(ro,0 <= u <= 1>) # control variable (supply)
    @variable(ro,0 <= ξ[i = 1:2] <= 1,RDDP.Uncertain) # uncertain parameter (demand)
    @constraint(ro,sum(ξ) <= 1>) # constraint abount uncertain parameter
    @constraint(ro, x.out == x.in - u + sum(ξ)) # state transition (temporal dependence of stock level)
    @objective(ro,Min,x) # minimize the stock
end
RDDP.train(msro)
```
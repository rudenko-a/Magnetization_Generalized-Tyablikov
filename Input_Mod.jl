module Input
using TOML

struct Parameters
    J1::Float64
    J2::Float64
    J3::Float64
    J4::Float64    
    D::Float64
    E::Float64
end

input = TOML.parsefile("parameters.toml")
param = Parameters(input["J1"], input["J2"], input["J3"], input["J4"], input["D"], input["E"])
const J = [param.J1 ; param.J2 ; param.J3 ; param.J4]
const D = param.D
const E = param.E
# D > 0: y-easy; E > 0: z-hard
# D > 0: y-easy; E < 0: x-hard

end

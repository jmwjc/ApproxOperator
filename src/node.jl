ParametricCoordinates = Union{Float64,NTuple{2,Float64},NTuple{3,Float64}}
## PhysicalNode
@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])
# ------------- Node -------------
struct Node <: PhysicalNode
    coordinates::NTuple{3,Float64}
end
Node(x::Float64,y::Float64,z::Float64) = Node((x,y,z))

## ParametricNode
# --------------- Gauss integration point ----------------
struct GaussPoint <: ParametricNode
    coordinates::ParametricCoordinates
    w::Float64
end
GaussPoint(ξ₁::Float64,ξ₂::Float64,w::Float64) = GaussPoint((ξ₁,ξ₂),w)
GaussPoint(ξ₁::Float64,ξ₂::Float64,ξ₃::Float64,w::Float64) = GaussPoint((ξ₁,ξ₂,ξ₃),w)

const QuadratureRule = Dict(
:PoiGI1 => [GaussPoint(-1.0,1.0)],
:SegGI1 => [GaussPoint(0.0,2.0)],
:SegGI2 =>
[
    GaussPoint(-0.5773502691896257645091487805,1.0),
    GaussPoint( 0.5773502691896257645091487805,1.0)
],
:SegGI3 =>
[
    GaussPoint(-0.774596669241483377035853079957,0.555555555555555555555555555556),
    GaussPoint( 0.0,0.88888888888888888888888888889),
    GaussPoint( 0.774596669241483377035853079957,0.555555555555555555555555555556)
],
:SegGI4 =>
[
    GaussPoint(-0.861136311594052575223946488893,0.347854845137453857373063949222),
    GaussPoint(-0.339981043584856264802665759103,0.652145154862546142626936050778),
    GaussPoint( 0.339981043584856264802665759103,0.652145154862546142626936050778),
    GaussPoint( 0.861136311594052575223946488893,0.347854845137453857373063949222)
],
:SegGI5 =>
[
    GaussPoint(-0.906179845938663992797626878299,0.23692688505618908751426404072),
    GaussPoint(-0.5384693101056830910363144207,0.47862867049936646804129151484),
    GaussPoint( 0.0,0.568888888888888888888888888889),
    GaussPoint( 0.5384693101056830910363144207,0.47862867049936646804129151484),
    GaussPoint( 0.906179845938663992797626878299,0.23692688505618908751426404072)
],
:SegGI6 =>
[
    GaussPoint(-0.932469514203152027812301554494,0.171324492379170345040296142173),
    GaussPoint(-0.6612093864662645136613995950,0.360761573048138607569833513838),
    GaussPoint(-0.238619186083196908630501721681,0.46791393457269104738987034399),
    GaussPoint( 0.238619186083196908630501721681,0.46791393457269104738987034399),
    GaussPoint( 0.66120938646626451366139959502,0.360761573048138607569833513838),
    GaussPoint( 0.932469514203152027812301554494,0.171324492379170345040296142173)
],
:SegGI7 =>
[
    GaussPoint(-0.949107912342758524526189684048,0.129484966168869693270611432679),
    GaussPoint(-0.741531185599394439863864773281,0.27970539148927666790146777142),
    GaussPoint(-0.405845151377397166906606412077,0.38183005050511894495036977549),
    GaussPoint( 0.0,0.417959183673469387755102040816),
    GaussPoint( 0.405845151377397166906606412077,0.381830050505118944950369775489),
    GaussPoint( 0.741531185599394439863864773281,0.279705391489276667901467771424),
    GaussPoint( 0.949107912342758524526189684048,0.129484966168869693270611432679)
],
:SegGI8 =>
[
    GaussPoint(-0.96028985649753623168356086857,0.10122853629037625915253135431),
    GaussPoint(-0.796666477413626739591553936476,0.22238103445337447054435599443),
    GaussPoint(-0.525532409916328985817739049189,0.313706645877887287337962201987),
    GaussPoint(-0.18343464249564980493947614236,0.36268378337836198296515044928),
    GaussPoint( 0.18343464249564980493947614236,0.362683783378361982965150449277),
    GaussPoint( 0.525532409916328985817739049189,0.31370664587788728733796220199),
    GaussPoint( 0.796666477413626739591553936476,0.222381034453374470544355994426),
    GaussPoint( 0.96028985649753623168356086857,0.10122853629037625915253135431)
],
:SegGI9 =>
[
    GaussPoint(-0.968160239507626089835576202904,0.0812743883615744119718921581105),
    GaussPoint(-0.83603110732663579429942978807,0.180648160694857404058472031243),
    GaussPoint(-0.613371432700590397308702039342,0.260610696402935462318742869419),
    GaussPoint(-0.32425342340380892903853801464,0.31234707704000284006863040658),
    GaussPoint( 0.0,0.330239355001259763164525069287),
    GaussPoint( 0.32425342340380892903853801464,0.31234707704000284006863040658),
    GaussPoint( 0.613371432700590397308702039342,0.260610696402935462318742869419),
    GaussPoint( 0.83603110732663579429942978807,0.180648160694857404058472031243),
    GaussPoint( 0.968160239507626089835576202904,0.081274388361574411971892158111)
],
:SegGI10 =>
[
    GaussPoint(-0.973906528517171720077964012085,0.066671344308688137593568809893),
    GaussPoint(-0.865063366688984510732096688424,0.149451349150580593145776339658),
    GaussPoint(-0.679409568299024406234327365115,0.219086362515982043995534934228),
    GaussPoint(-0.433395394129247190799265943166,0.26926671930999635509122692157),
    GaussPoint(-0.14887433898163121088482600113,0.295524224714752870173892994651),
    GaussPoint( 0.14887433898163121088482600113,0.295524224714752870173892994651),
    GaussPoint( 0.433395394129247190799265943166,0.26926671930999635509122692157),
    GaussPoint( 0.679409568299024406234327365115,0.219086362515982043995534934228),
    GaussPoint( 0.865063366688984510732096688424,0.149451349150580593145776339658),
    GaussPoint( 0.973906528517171720077964012085,0.066671344308688137593568809893)
],
:TriGI1 =>
[
    GaussPoint(0.0,0.0,1.0)
],
:TriGI3 =>
[
    GaussPoint(2/3,1/6,1/3),
    GaussPoint(1/6,2/3,1/3),
    GaussPoint(1/6,1/6,1/3)
],
:TriGI4 =>
[
    GaussPoint(0.333333333333333,0.333333333333333,-0.562500000000000),
    GaussPoint(0.600000000000000,0.200000000000000, 0.520833333333333),
    GaussPoint(0.200000000000000,0.600000000000000, 0.520833333333333),
    GaussPoint(0.200000000000000,0.200000000000000, 0.520833333333333)
],
:TriGI6 =>
[
    GaussPoint(0.108103018168070,0.445948490915965,0.223381589678011),
    GaussPoint(0.445948490915965,0.108103018168070,0.223381589678011),
    GaussPoint(0.445948490915965,0.445948490915965,0.223381589678011),
    GaussPoint(0.816847572980459,0.091576213509771,0.109951743655322),
    GaussPoint(0.091576213509771,0.816847572980459,0.109951743655322),
    GaussPoint(0.091576213509771,0.091576213509771,0.109951743655322)
],
:TriGI7 =>
[
    GaussPoint(0.101286507323500,0.101286507323500,0.125939180544800),
    GaussPoint(0.797426985353100,0.101286507323500,0.125939180544800),
    GaussPoint(0.101286507323500,0.797426985353100,0.125939180544800),
    GaussPoint(0.470142064105100,0.059715871789800,0.132394152788500),
    GaussPoint(0.470142064105100,0.470142064105100,0.132394152788500),
    GaussPoint(0.059715871789800,0.470142064105100,0.132394152788500),
    GaussPoint(0.333333333333300,0.333333333333300,0.225000000000000)
],
:TriGI12 =>
[
    GaussPoint(0.501426509658179,0.249286745170910,0.116786275726379),
    GaussPoint(0.249286745170910,0.501426509658179,0.116786275726379),
    GaussPoint(0.249286745170910,0.249286745170910,0.116786275726379),
    GaussPoint(0.873821971016996,0.063089014491502,0.050844906370207),
    GaussPoint(0.063089014491502,0.873821971016996,0.050844906370207),
    GaussPoint(0.063089014491502,0.063089014491502,0.050844906370207),
    GaussPoint(0.053145049844817,0.310352451033784,0.082851075618374),
    GaussPoint(0.053145049844817,0.636502499121399,0.082851075618374),
    GaussPoint(0.310352451033784,0.053145049844817,0.082851075618374),
    GaussPoint(0.310352451033784,0.636502499121399,0.082851075618374),
    GaussPoint(0.636502499121399,0.053145049844817,0.082851075618374),
    GaussPoint(0.636502499121399,0.310352451033784,0.082851075618374)
],
:TriGI13 =>
[
    GaussPoint(0.065130102902200,0.065130102902200,0.053347235608800),
    GaussPoint(0.869739794195600,0.065130102902200,0.053347235608800),
    GaussPoint(0.065130102902200,0.869739794195600,0.053347235608800),
    GaussPoint(0.312865496004900,0.048690315425300,0.077113760890300),
    GaussPoint(0.638444188569800,0.312865496004900,0.077113760890300),
    GaussPoint(0.048690315425300,0.638444188569800,0.077113760890300),
    GaussPoint(0.638444188569800,0.048690315425300,0.077113760890300),
    GaussPoint(0.312865496004900,0.638444188569800,0.077113760890300),
    GaussPoint(0.048690315425300,0.312865496004900,0.077113760890300),
    GaussPoint(0.260345966079000,0.260345966079000,0.175615257433200),
    GaussPoint(0.479308067841900,0.260345966079000,0.175615257433200),
    GaussPoint(0.260345966079000,0.479308067841900,0.175615257433200),
    GaussPoint(0.333333333333300,0.333333333333300,-0.14957004446770)
],
:TriGI16 =>
[
    GaussPoint(0.333333333333333,0.333333333333333,0.144315607677787),
    GaussPoint(0.081414823414554,0.459292588292723,0.095091634267285),
    GaussPoint(0.459292588292723,0.081414823414554,0.095091634267285),
    GaussPoint(0.459292588292723,0.459292588292723,0.095091634267285),
    GaussPoint(0.658861384496480,0.170569307751760,0.103217370534718),
    GaussPoint(0.170569307751760,0.658861384496480,0.103217370534718),
    GaussPoint(0.170569307751760,0.170569307751760,0.103217370534718),
    GaussPoint(0.898905543365938,0.050547228317031,0.032458497623198),
    GaussPoint(0.050547228317031,0.898905543365938,0.032458497623198),
    GaussPoint(0.050547228317031,0.050547228317031,0.032458497623198),
    GaussPoint(0.008394777409958,0.263112829634638,0.027230314174435),
    GaussPoint(0.008394777409958,0.728492392955404,0.027230314174435),
    GaussPoint(0.263112829634638,0.008394777409958,0.027230314174435),
    GaussPoint(0.263112829634638,0.728492392955404,0.027230314174435),
    GaussPoint(0.728492392955404,0.008394777409958,0.027230314174435),
    GaussPoint(0.728492392955404,0.263112829634638,0.027230314174435)
],
:QuadGI1 =>
[
    GaussPoint(0.0,0.0,2.0)
],
:QuadGI2 =>
[
    GaussPoint(-0.5773502691896258,-0.5773502691896258,1.0),
    GaussPoint( 0.5773502691896258,-0.5773502691896258,1.0),
    GaussPoint( 0.5773502691896258, 0.5773502691896258,1.0),
    GaussPoint(-0.5773502691896258, 0.5773502691896258,1.0)
]
)

## Actions
getindex(x::AbstractNode,i::Int) = x.coordinates[i]

# set_integration_rule!
function set_integration_rule!(ap::Approximator,𝓖::Symbol)
    ap.𝓖 = QuadratureRule[𝓖]
end
function set_integration_rule!(aps::Vector{Approximator},𝓖::Symbol)
    for ap in aps
        set_integration_rule!(ap,𝓖)
    end
end

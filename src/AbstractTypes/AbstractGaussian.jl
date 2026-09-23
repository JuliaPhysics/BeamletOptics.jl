"""From Wilhelm (2001)"""
isparaxial(w0, λ; margin::Int=100) = w0 / λ > margin ? true : false
